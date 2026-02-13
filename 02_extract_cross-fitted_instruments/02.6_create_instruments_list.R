# Purpose: Create a list of instrumental variables (IVs) for the dosage extraction.

# Clean the environment
rm(list=ls())

# Paths
if (Sys.getenv("GRAPH_MR_ROOT") == "") stop("Set GRAPH_MR_ROOT to the pipeline root directory (see README).")
setwd(Sys.getenv("GRAPH_MR_ROOT"))
source("config/config.R")

# Load packages
library("tidyverse")
library("vroom")
library("parallel")
library("gtools")
library("ensembldb")
library("EnsDb.Hsapiens.v86")
library("AnnotationFilter")
library("GenomicRanges")

# Set the number of cores you require
ncores <- 8

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

# Optional: protein annotation for cis-region filtering (proteomics; set GRAPH_MR_ANNOTATION_DIR or use data/annotation/OlinkProteomics/)
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
  vroom(protein38_path, col_names = c("UniProtID", "hg38chr", "hg38start", "hg38end")) %>% distinct(.)
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

# Merge protein coordinates and for any duplicates, pick the longest coordinate duplicate
merged_protein_coords <- dplyr::inner_join(protein_19_38_coords, protein38coords, by = c("hg38chr", "hg38start", "hg38end")) %>% 
  dplyr::mutate(hg19_length = hg19end - hg19start) %>%      
  dplyr::group_by(UniProtID) %>%
  dplyr::slice_max(order_by = hg19_length, n = 1, with_ties = F) %>%  
  dplyr::ungroup()

# Set the output file paths (will write to all method_types if METHOD_TYPES is set)
# For reading, use the primary method_type
output <- paste0(results_dir, base_dir, "/instruments/")

# Create output directories for all method types
output_dirs <- sapply(method_types_list, function(mt) {
  paste0(results_dir, mt, "_", input_type, "/instruments/")
})
output_summary_dirs <- sapply(method_types_list, function(mt) {
  dir_path <- paste0(results_dir, mt, "_", input_type, "/summary_stats/")
  dir.create(dir_path, recursive = T, showWarnings = F)
  return(dir_path)
})



# Start script ---- ---- ---- ---- ---- ----
writeLines(c(" ", " ", "Started script 02.6 - Make instruments list.", " ", " "))

# Initialize empty data frame
nemptyfolds <- data.frame(Exposure = character(), N.empty_rows = integer())

# Create SNP lists for the exposures with <5 empty folds
create_iv_lists <- function(input_list) {
  mclapply(input_list, mc.cores = ncores, FUN = function(exposure_i) {
    
    # List the 10 SNP result files by exposure 
    fold_files <- list.files(path = paste0(output, exposure_i), full.names = TRUE, recursive = TRUE, pattern = "fold[0-9]+_instruments\\.tsv$")
    fold_files <- fold_files[mixedorder(fold_files)]
    
    # Read the 10 folds as a list of lists
    read_data <- function(fold_files) {
      vroom(fold_files, delim = "\t", 
            col_types = c("ID" = "c", "CHROM" = "i", "POS" = "i", 
                          "BETA" = "d", "SE" = "d", "P" = "d", 
                          "REF" = "c", "ALT" = "c"))
    }
    
    fold_data_list <- lapply(fold_files, read_data)
    empty_folds_count <- sum(sapply(fold_data_list, nrow) == 0)
    
    fold_data_mutated <- list() 
    
    if (empty_folds_count < 5) {
      fold_data_mutated <- lapply(fold_data_list, function(fold_i) {
        if (nrow(fold_i) > 0) {
          # Remove variants that contain >1 base
          fold_i <- fold_i %>%
            dplyr::filter(nchar(REF) == 1 & nchar(ALT) == 1)
          
          if (nrow(fold_i) > 0) {
            # Case 1: Protein with coordinates available — filter for cis region
            if (exposure_i %in% protein_names$UniProt && exposure_i %in% merged_protein_coords$UniProtID) {
              coords <- merged_protein_coords %>%
                dplyr::filter(UniProtID == exposure_i) %>%
                dplyr::mutate(hg19start = hg19start - 500000,
                              hg19end = hg19end + 500000)
              
              mutated_data <- fold_i %>%
                dplyr::filter(CHROM == coords$hg19chr,
                              POS >= coords$hg19start,
                              POS <= coords$hg19end) %>%
                dplyr::mutate(list = paste0(CHROM, "\t", POS)) %>%
                dplyr::select(list)
              
              if (nrow(mutated_data) == 0) return(NULL)
              return(mutated_data)
              
              # Case 2: Protein but no coordinates — skip this fold
            } else if (exposure_i %in% protein_names$UniProt) {
              return(NULL)
              
              # Case 3: Not a protein — use full variant list
            } else {
              mutated_data <- fold_i %>%
                dplyr::mutate(list = paste0(CHROM, "\t", POS)) %>%
                dplyr::select(list)
              
              return(mutated_data)
            }
          } else {
            return(NULL)
          }
        } else {
          return(NULL)
        }
      })
      
      # Save the lists as .txt files (write to all method_type folders)
      for (i in seq_along(fold_data_mutated)) {
        if (!is.null(fold_data_mutated[[i]])) {
          for (output_dir in output_dirs) {
            dir.create(paste0(output_dir, exposure_i), recursive = T, showWarnings = F)
            write.table(fold_data_mutated[[i]],
                        paste0(output_dir, exposure_i, "/fold", i, "_instruments_list.txt"),
                        col.names = FALSE, row.names = FALSE, quote = FALSE)
          }
        }
      }
    }
    
    # Print message and return metadata
    cat(c(exposure_i, "-", "N. empty folds:", empty_folds_count, "\n"))
    new_row <- data.frame(Exposure = exposure_i, N.empty_rows = empty_folds_count)
    return(new_row)
  })
}


iv_lists <- create_iv_lists(input_list)

nemptyfolds <- base::do.call(rbind, iv_lists)
# Write summary to all method_type folders
for (output_summary_dir in output_summary_dirs) {
  vroom_write(nemptyfolds, paste0(output_summary_dir, "n_empty_folds.csv"), delim = ",", append = F)
}

writeLines(c(" ", " ", "Finished script 02.6 - Made instruments list.", " ", " "))
