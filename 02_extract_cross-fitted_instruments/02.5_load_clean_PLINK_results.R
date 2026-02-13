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
    # List, read, and bind the 22 GWAS result files for the current exposure fold
    gwas_res_filelist <- list.files(path = paste0(input_files, input_i, "/fold", fold_i), pattern = "*glm.linear", full.names = T, recursive = T)
    gwas_res_fold_i <- rbindlist(l = lapply(gwas_res_filelist, fread, 
                                            select = c("ID", "#CHROM", "POS", "BETA", "SE", "P", "REF", "ALT"),
                                            col.names = c("ID", "CHROM", "POS", "BETA", "SE", "P", "REF", "ALT")), 
                                 use.names = T, fill = T) 
    # Filter out rows where the column "ID" contains the value "ID"
    gwas_res_fold_i <- gwas_res_fold_i[gwas_res_fold_i$ID != "ID", ]
    # If GWAS result is empty, create a dataframe, if not, proceed with binding and saving the results
    if (nrow(gwas_res_fold_i) == 0) {
      # Create an empty dataframe with specified columns and save it as a file so that the number of files remains n_phenotypes * k_folds)
      gwas_res_fold_i <- data.frame(ID = character(),
                                    CHROM = character(),
                                    POS = integer(),
                                    BETA = numeric(),
                                    SE = numeric(),
                                    P = numeric(),
                                    REF = character(),
                                    ALT = character(),
                                    stringsAsFactors = FALSE)
      # Write to all method_type folders
      for (output_dir in output_dirs) {
        vroom_write(gwas_res_fold_i, col_names = T, file = paste0(output_dir, input_i, "/fold", fold_i, "_instruments.tsv"))
      }
      # Also create a clumped_res_fold_i object with 0 rows
      clumped_res_fold_i <- data.frame(SNP = character())
      } else {
      # List, read, bind, and save the clumped results for the current exposure fold
      clumped_res_filelist <- list.files(path = paste0(input_files, input_i, "/fold", fold_i), pattern = "*.clumped.clumps", full.names = T, recursive = T)
      clumped_res_fold_i <- rbindlist(l = lapply(clumped_res_filelist, fread, select = c("ID")), use.name = T, fill = T)
      clumped_res_fold_i <- clumped_res_fold_i %>%
        inner_join(y = gwas_res_fold_i, by = "ID") 
      # Write to all method_type folders
      for (output_dir in output_dirs) {
        vroom_write(clumped_res_fold_i, col_names = T, file = paste0(output_dir, input_i, "/fold", fold_i, "_instruments.tsv"))
      }
      }
    # Save the number of SNPs selected by the GWAS and after clumping in a dataframe
    gwas_snps_dim <- data.frame(Exposure = input_i,
                                Fold = paste0("fold", fold_i),
                                "n_gwas_snps" = nrow(gwas_res_fold_i),
                                "n_clumped_snps" = nrow(clumped_res_fold_i),
                                check.names = F)
    # Write to all method_type folders
    for (output_dir in output_dirs) {
      vroom_write(gwas_snps_dim, col_names = T, file = paste0(output_dir, input_i, "/fold", fold_i, "_instrument_dimensions.csv"))
    }
    # Print current exposure and fold
    writeLines(c(" ", " ", "Joined PLINK results for:", input_i, fold_i, " ", " "))
    }
  )
}

saved_gwas_res <- save_gwas(input_i)

writeLines(c(" ", " ", "Finished script 02.5 - finished PLINK data loading and cleaning", " ", " "))
