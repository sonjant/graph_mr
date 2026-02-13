# Purpose: Split each exposure (X) into complete and missing cases. 

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

# Set the number of cores, typically one per fold (k)
ncores <- 8

# Arguments: method_type, dataset
args <- commandArgs(trailingOnly = TRUE)
method_type <- args[1]
dataset     <- if (length(args) >= 2L) args[2] else Sys.getenv("GRAPH_MR_DATASET")
if (!nzchar(dataset)) stop("Provide dataset name or set GRAPH_MR_DATASET.")
base_dir <- paste0(method_type, "_", dataset)

exposures_path <- paste0(data_dir, base_dir, "/preprocessed/exposures.tsv")
if (!file.exists(exposures_path)) stop("Preprocessed exposures not found: ", exposures_path)
input_df <- vroom(exposures_path)
output <- paste0(data_dir, base_dir, "/gwas_input/")


# Start script ---- ---- ---- ---- ---- ----
writeLines(c(" ", " ", "Started script 03.1 - Create missing and complete cases lists", " ", " "))

# For each exposure, filter missing cases and create a text file for running a GWAS
create_exposure_dfs <- function(input_df) {
  exposure_names <- colnames(input_df)[-which(colnames(input_df) == "barcode")]
  create_exposure_i_df <- function(exposure_name) {
    # Extract and format each exposure for running a GWAS  
    exposure_data <- input_df %>%
      mutate(fullbarcode = barcode) %>% 
      dplyr::select(barcode, fullbarcode, exposure_name) 
    # Remove rows with NAs
    exposure_data_no_na <- exposure_data %>% 
      filter(!is.na(.[[exposure_name]]))
    # Standardise exposure data
    exposure_data_no_na[, 3] <- scale(exposure_data_no_na[, 3])
    # Create list of participants with missing data for the current exposure
    exposure_data_only_na <- exposure_data %>% 
      filter(is.na(.[[exposure_name]])) %>% 
      mutate(eid = barcode) %>% 
      dplyr::select(eid)
    # Save the complete cases as a text file
    gwas_participants_file_name <- paste0(output, exposure_name, "/complete_cases_gwas_participants.txt")
    write.table(exposure_data_no_na, file = gwas_participants_file_name, row.names = F, col.names = F, quote = F, sep = "\t")
    # Save missing cases as a text file
    na_participants_file_name <- paste0(output, exposure_name, "/participants_missing_exposure.txt")
    write.table(exposure_data_only_na, file = na_participants_file_name, row.names = F, col.names = F, quote = F, sep = "\t")
    cat("Finished:", exposure_name, "\n")
  }
  mclapply(exposure_names, create_exposure_i_df, mc.cores = ncores)
}

# Process each exposure
exposure_dfs <- create_exposure_dfs(input_df)



writeLines(c(" ", " ", "Finished script 03.1 - Create missing and complete cases lists.", " ", " "))

