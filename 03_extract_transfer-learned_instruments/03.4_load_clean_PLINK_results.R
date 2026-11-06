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

# Specify a base_dir based on the command line inputs
args <- commandArgs(trailingOnly = TRUE)
method_type <- args[1]
input_type <- args[2]
base_dir <- paste0(method_type, "_", input_type)

# Set the number of cores you require
ncores <- 6

# Set the input_files path to the "output" path from script 03.3
input_list <- list.dirs(path = paste0(ephemeral_dir, "results/", base_dir, "/gwas_results"))
input_list <- grep("complete_cases$", input_list, value = T)

# Set an output file path for the selected instrumental variables
output <- paste0(results_dir, base_dir, "/instruments/")

# Protein cis window, matching scripts 02.5 and 03.5: hg19 gene body ± 500 kb.
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
writeLines(c(" ", " ", "Started script 03.4 - PLINK data loading and cleaning", " ", " "))

# Bind and save the GWAS result files for each exposure
save_gwas <- function(input_list){
  # For each exposure, bind and save the GWAS result files
  mclapply(input_list, mc.cores = ncores, FUN = function(input_i){
    # Extract the name of the input exposure
    input_name <- basename(dirname(input_i))
    is_protein <- input_name %in% protein_names$UniProt
    cis_window <- merged_protein_coords %>%
      dplyr::filter(UniProtID == input_name) %>%
      dplyr::mutate(hg19start = hg19start - 500000,
                    hg19end = hg19end + 500000)
    has_cis_coords <- is_protein && nrow(cis_window) == 1
    # Create output directory
    dir.create(paste0(output, input_name), showWarnings = F, recursive = T)
    # List, read, and bind the 22 GWAS result files for the current exposure
    gwas_res_filelist <- list.files(path = input_i, pattern = "*glm.linear", full.names = T, recursive = T)
    gwas_res <- rbindlist(l = lapply(gwas_res_filelist, fread, 
                                     select = c("ID", "#CHROM", "POS", "BETA", "SE", "P", "REF", "ALT"),
                                     col.names = c("ID", "CHROM", "POS", "BETA", "SE", "P", "REF", "ALT")), 
                          use.names = T, fill = T) 
    # Filter out rows where the column "ID" contains the value "ID"
    gwas_res <- gwas_res[gwas_res$ID != "ID", ]
    n_gwas_snps <- nrow(gwas_res)

    # If GWAS result is empty, keep an empty instrument table.
    # Otherwise join the clumped SNP IDs back to the GWAS rows.
    if (n_gwas_snps == 0) {
      n_clumped_snps <- 0
      instruments <- empty_instruments()
      n_clumped_cis_snps <- if (is_protein) 0L else NA_integer_
    } else {
      clumped_res_filelist <- list.files(path = input_i, pattern = "*.clumped.clumps", full.names = T, recursive = T)
      clumped_res <- rbindlist(l = lapply(clumped_res_filelist, fread, select = c("ID")), use.name = T, fill = T)
      clumped_res <- clumped_res %>%
        inner_join(y = gwas_res, by = "ID")
      n_clumped_snps <- nrow(clumped_res)

      # Proteins: save cis SNPs only. A protein with no coordinates has no cis window, so its instrument file is empty.
      # Other exposures: save the genome-wide clumped SNPs.
      if (has_cis_coords) {
        instruments <- filter_to_cis(clumped_res, cis_window)
        n_clumped_cis_snps <- nrow(instruments)
      } else if (is_protein) {
        instruments <- empty_instruments()
        n_clumped_cis_snps <- 0L
      } else {
        instruments <- clumped_res
        n_clumped_cis_snps <- NA_integer_
      }
    }

    vroom_write(instruments, col_names = T, file = paste0(output, "/", input_name, "/complete_cases_instruments.tsv"))
    gwas_snps_dim <- data.frame(Exposure = input_name,
                                n_gwas_snps = n_gwas_snps,
                                n_clumped_snps = n_clumped_snps,
                                n_clumped_cis_snps = n_clumped_cis_snps,
                                check.names = F)
    vroom_write(gwas_snps_dim, col_names = T, file = paste0(output, "/", input_name, "/complete_cases_instruments_dims.csv"))
    # Print the current exposure
    writeLines(c(" ", " ", "Joined PLINK results for:", input_i, " ", " "))
    }
  )
}

saved_gwas_res <- save_gwas(input_list)

writeLines(c(" ", " ", "Finished script 03.4 - finished PLINK data loading and cleaning", " ", " "))
