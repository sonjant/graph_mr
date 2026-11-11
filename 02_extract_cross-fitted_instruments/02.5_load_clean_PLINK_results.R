# Purpose: Load and clean PLINK GWAS results.

# Clean the environment
rm(list=ls())

# Paths
if (Sys.getenv("GRAPH_MR_ROOT") == "") stop("Set GRAPH_MR_ROOT to the pipeline root directory (see README).")
setwd(Sys.getenv("GRAPH_MR_ROOT"))
source("config/config.R")

# Load packages
library("tidyverse")
library("vroom")
library("data.table")
library("parallel")

# Set the number of cores you require
ncores <- 2

# Start array indexing
task_id_string <- Sys.getenv("PBS_ARRAY_INDEX")
k <- as.numeric(task_id_string)

# Specify a base_dir based on the command line inputs
args <- commandArgs(trailingOnly = TRUE)
method_type <- args[1]
input_type <- args[2]
base_dir <- paste0(method_type, "_", input_type)

# Check if METHOD_TYPES environment variable is set (for writing to multiple method folders)
method_types_str <- Sys.getenv("METHOD_TYPES", unset = "")
if (method_types_str != "") {
  method_types_list <- strsplit(method_types_str, ",")[[1]]
  method_types_list <- trimws(method_types_list)  # Remove any whitespace
} else {
  method_types_list <- c(method_type)  # Use single method_type if not specified
}

# List input (use primary method_type for reading)
input_list <- basename(list.dirs(paste0(data_dir, base_dir, "/gwas_input/"), full.names = TRUE))[-1]

# Dedicate one PBS array element to each element in the input_list 
input_i <- input_list[k]

# Set the input_files path to the "output" path from script 02.3 and 02.4
input_files <- paste0(ephemeral_dir, "results/", base_dir, "/gwas_results/")

# Set output file paths (will write to all method_types if METHOD_TYPES is set)
# For reading, use the primary method_type
output <- paste0(results_dir, base_dir, "/instruments/")

# Create output directories for all method types
output_dirs <- sapply(method_types_list, function(mt) {
  paste0(results_dir, mt, "_", input_type, "/instruments/")
})

# Protein cis window, matching script 02.6: hg19 gene body ± 500 kb.
# Non-protein exposures keep the genome-wide clumped set.
olink_csv <- paste0(annotation_dir, "OlinkProteomics/olink-explore-3072-assay-list-2024-03-19.csv")
protein_names <- if (file.exists(olink_csv)) {
  vroom(olink_csv) %>%
    dplyr::rename(UniProt = "UniProt ID", protein_name = "Protein name") %>%
    dplyr::select(UniProt, protein_name)
} else {
  data.frame(UniProt = character(), protein_name = character())
}
protein38_path <- paste0(annotation_dir, "OlinkProteomics/proteinID_hg38_coordinates.txt")
protein38coords <- if (file.exists(protein38_path)) {
  vroom(protein38_path, col_names = c("UniProtID", "hg38chr", "hg38start", "hg38end")) %>%
    distinct(.)
} else {
  data.frame(UniProtID = character(), hg38chr = integer(), hg38start = integer(), hg38end = integer())
}
bed_path <- paste0(annotation_dir, "OlinkProteomics/protein_coords_hg19_hg38.bed")
protein_19_38_coords <- if (file.exists(bed_path)) {
  vroom(bed_path, col_names = c("hg19chr", "hg19start", "hg19end", "hg38chr_start_end", "rm")) %>%
    dplyr::select("hg19chr", "hg19start", "hg19end", "hg38chr_start_end") %>%
    dplyr::mutate(hg19chr = str_remove(hg19chr, "^chr")) %>%
    tidyr::separate(hg38chr_start_end, into = c("hg38chr", "pos"), sep = ":") %>%
    tidyr::separate(pos, into = c("hg38start", "hg38end"), sep = "-") %>%
    dplyr::mutate(hg38chr = str_remove(hg38chr, "^chr"),
                  hg38chr = as.integer(hg38chr),
                  hg38start = as.integer(hg38start),
                  hg38end = as.integer(hg38end),
                  hg19chr = as.integer(hg19chr)) %>%
    distinct(.)
} else {
  data.frame(hg19chr = integer(), hg19start = integer(), hg19end = integer(), hg38chr = integer(), hg38start = integer(), hg38end = integer())
}

merged_protein_coords <- dplyr::inner_join(protein_19_38_coords, protein38coords, by = c("hg38chr", "hg38start", "hg38end")) %>%
  dplyr::mutate(hg19_length = hg19end - hg19start) %>%
  dplyr::group_by(UniProtID) %>%
  dplyr::slice_max(order_by = hg19_length, n = 1, with_ties = F) %>%
  dplyr::ungroup()

is_protein <- input_i %in% protein_names$UniProt
cis_window <- merged_protein_coords %>%
  dplyr::filter(UniProtID == input_i) %>%
  dplyr::mutate(hg19start = hg19start - 500000,
                hg19end = hg19end + 500000)
has_cis_coords <- is_protein && nrow(cis_window) == 1

empty_instruments <- function() {
  data.frame(ID = character(),
             CHROM = character(),
             POS = integer(),
             BETA = numeric(),
             SE = numeric(),
             P = numeric(),
             REF = character(),
             ALT = character(),
             stringsAsFactors = FALSE)
}

# Linear results are *.glm.linear. Logistic results are *.glm.logistic,
# *.glm.logistic.hybrid (PLINK default, Firth fallback), or *.glm.firth.
# One chromosome can have a leftover linear file after a logistic rerun;
# keep the logistic/Firth file when both exist.
glm_result_files <- function(path) {
  candidates <- list.files(path = path,
                           pattern = "PHENO1\\.glm\\.(linear|logistic|firth)",
                           full.names = TRUE,
                           recursive = TRUE)
  if (length(candidates) == 0) {
    return(candidates)
  }
  preference <- function(candidate) {
    name <- basename(candidate)
    if (grepl("glm\\.logistic\\.hybrid", name)) {
      1
    } else if (grepl("glm\\.logistic", name)) {
      2
    } else if (grepl("glm\\.firth", name)) {
      3
    } else {
      4
    }
  }
  chromosome <- sub("\\.PHENO1\\.glm\\..*$", "", basename(candidates))
  keep <- tapply(seq_along(candidates), chromosome, function(idx) {
    idx[[which.min(vapply(candidates[idx], preference, numeric(1)))]]
  })
  candidates[sort(unlist(keep, use.names = FALSE))]
}

# Logistic output reports the log-odds as BETA when cols=+beta is set, and
# names the standard error LOG(OR)_SE when an odds ratio is also reported.
# Scores downstream still use BETA, which is the log-odds coefficient.
read_glm_results <- function(path) {
  header <- names(fread(path, nrows = 0))
  chrom_col <- if ("#CHROM" %in% header) "#CHROM" else "CHROM"
  se_col <- if ("SE" %in% header) "SE" else "LOG(OR)_SE"
  beta_from_or <- !("BETA" %in% header)
  beta_col <- if (beta_from_or) "OR" else "BETA"
  needed <- c("ID", chrom_col, "POS", beta_col, se_col, "P", "REF", "ALT")
  missing_cols <- setdiff(needed, header)
  if (length(missing_cols) > 0) {
    stop("GLM file is missing columns: ", paste(missing_cols, collapse = ", "), " (", path, ")")
  }
  result <- fread(path,
                  select = needed,
                  col.names = c("ID", "CHROM", "POS", "BETA", "SE", "P", "REF", "ALT"))
  if (beta_from_or) {
    result[, BETA := log(as.numeric(BETA))]
  }
  result
}

# Keep clumped SNPs inside the protein cis window.
filter_to_cis <- function(clumped_snps, cis_window) {
  clumped_snps %>%
    dplyr::mutate(CHROM = as.integer(str_remove(as.character(CHROM), "^chr")),
                  POS = as.integer(POS)) %>%
    dplyr::filter(CHROM == cis_window$hg19chr,
                  POS >= cis_window$hg19start,
                  POS <= cis_window$hg19end)
}



# Start script ---- ---- ---- ---- ---- ----
writeLines(c(" ", " ", "Started script 02.5 - PLINK data loading and cleaning", " ", " "))

# Bind and save the GWAS result files for each exposure fold
save_gwas <- function(input_i, nfolds = 10){
  # For each exposure fold, bind and save the GWAS result files
  mclapply(1:nfolds, mc.cores = ncores, FUN = function(fold_i){
    # Create result directories for the current exposure in all method_type folders
    for (output_dir in output_dirs) {
      dir.create(paste0(output_dir, input_i), showWarnings = F, recursive = T)
    }
    # List, read, and bind the 22 GWAS result files for the current exposure fold.
    # A missing file (for example a skipped constant binary fold) leaves an empty table.
    gwas_res_filelist <- glm_result_files(paste0(input_files, input_i, "/fold", fold_i))
    if (length(gwas_res_filelist) == 0) {
      gwas_res_fold_i <- empty_instruments()
    } else {
      gwas_res_fold_i <- rbindlist(l = lapply(gwas_res_filelist, read_glm_results),
                                   use.names = T, fill = T)
      # Filter out rows where the column "ID" contains the value "ID"
      gwas_res_fold_i <- gwas_res_fold_i[gwas_res_fold_i$ID != "ID", ]
    }
    n_gwas_snps <- nrow(gwas_res_fold_i)

    # If GWAS result is empty, keep an empty instrument table.
    # Otherwise join the clumped SNP IDs back to the GWAS rows.
    if (n_gwas_snps == 0) {
      n_clumped_snps <- 0
      instruments_fold_i <- empty_instruments()
      n_clumped_cis_snps <- if (is_protein) 0L else NA_integer_
    } else {
      clumped_res_filelist <- list.files(path = paste0(input_files, input_i, "/fold", fold_i), pattern = "*.clumped.clumps", full.names = T, recursive = T)
      clumped_res_fold_i <- rbindlist(l = lapply(clumped_res_filelist, fread, select = c("ID")), use.name = T, fill = T)
      clumped_res_fold_i <- clumped_res_fold_i %>%
        inner_join(y = gwas_res_fold_i, by = "ID")
      n_clumped_snps <- nrow(clumped_res_fold_i)

      # Proteins: save cis SNPs only. A protein with no coordinates has no cis window, so its instrument file is empty.
      # Other exposures: save the genome-wide clumped SNPs.
      if (has_cis_coords) {
        instruments_fold_i <- filter_to_cis(clumped_res_fold_i, cis_window)
        n_clumped_cis_snps <- nrow(instruments_fold_i)
      } else if (is_protein) {
        instruments_fold_i <- empty_instruments()
        n_clumped_cis_snps <- 0L
      } else {
        instruments_fold_i <- clumped_res_fold_i
        n_clumped_cis_snps <- NA_integer_
      }
    }

    # Write the instrument table and the SNP counts to all method_type folders
    gwas_snps_dim <- data.frame(Exposure = input_i,
                                Fold = paste0("fold", fold_i),
                                n_gwas_snps = n_gwas_snps,
                                n_clumped_snps = n_clumped_snps,
                                n_clumped_cis_snps = n_clumped_cis_snps,
                                check.names = F)
    for (output_dir in output_dirs) {
      vroom_write(instruments_fold_i, col_names = T, file = paste0(output_dir, input_i, "/fold", fold_i, "_instruments.tsv"))
      vroom_write(gwas_snps_dim, col_names = T, file = paste0(output_dir, input_i, "/fold", fold_i, "_instrument_dimensions.csv"))
    }
    # Print current exposure and fold
    writeLines(c(" ", " ", "Joined PLINK results for:", input_i, fold_i, " ", " "))
    }
  )
}

saved_gwas_res <- save_gwas(input_i)

writeLines(c(" ", " ", "Finished script 02.5 - finished PLINK data loading and cleaning", " ", " "))
