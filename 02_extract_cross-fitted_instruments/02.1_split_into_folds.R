# Purpose: Split exposure (X) data into folds containing train and test sets. Also save the participants with missing data for each exposure.

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

# Arguments: method_type, dataset (dataset = any label, e.g. proteomics, mydata)
args <- commandArgs(trailingOnly = TRUE)
method_type <- args[1]
dataset     <- if (length(args) >= 2L) args[2] else Sys.getenv("GRAPH_MR_DATASET")
if (!nzchar(dataset)) stop("Provide dataset name as second argument or set GRAPH_MR_DATASET (e.g. proteomics, mydata).")
base_dir <- paste0(method_type, "_", dataset)

# Check if METHOD_TYPES environment variable is set (for writing to multiple method folders)
method_types_str <- Sys.getenv("METHOD_TYPES", unset = "")
if (method_types_str != "") {
  method_types_list <- strsplit(method_types_str, ",")[[1]]
  method_types_list <- trimws(method_types_list)
} else {
  method_types_list <- c(method_type)
}

# Read preprocessed exposures (generic: one file per dataset)
exposures_path <- paste0(data_dir, base_dir, "/preprocessed/exposures.tsv")
if (!file.exists(exposures_path)) {
  stop("Preprocessed exposures not found: ", exposures_path, ". Run generic preprocessing or create layout (see docs/INPUT_SPEC.md).")
}
input_df <- vroom(exposures_path)

# Output directories for all method types
output_dirs <- sapply(method_types_list, function(mt) {
  paste0(data_dir, mt, "_", dataset, "/gwas_input/")
})



# Start script ---- ---- ---- ---- ---- ----
writeLines(c(" ", " ", "Started script 02.1 - Split exposures into folds", " ", " "))

# Split each exposure into folds, process folds, and save as text files for running GWASes
split_exposures <- function(input_df, k = 10, unique_threshold = 4) {
  # Extract exposure names
  exposure_names <- colnames(input_df)[-which(colnames(input_df) == "barcode")]
  # Split each exposure into k folds
  split_exposure_i <- function(exposure_name) {
    print(exposure_name)
    # Extract and format each exposure for running a GWAS  
    exposure_data <- input_df %>%
      mutate(fullbarcode = barcode) %>% 
      dplyr::select(barcode, fullbarcode, all_of(exposure_name))
    # Remove rows with NAs for splitting into folds
    exposure_data_no_na <- exposure_data %>% 
      filter(!is.na(.[[exposure_name]]))
    # Create 10-fold cross-validation indices
    create_folds <- function(data) {
      set.seed(123)
      n <- nrow(data)
      sample_indices <- sample(seq_len(n))
      split(sample_indices, cut(seq_len(n), breaks = k, labels = FALSE))
    }
    folds <- create_folds(exposure_data_no_na)
    # Create train and test sets within each exposure fold
    process_fold <- function(fold_indices, fold_number) {
      train_indices <- unlist(folds[-fold_number])
      test_indices <- fold_indices
      train_set <- exposure_data_no_na[train_indices, ]
      test_set <- exposure_data_no_na[test_indices, ]
      # If the train set is continuous, standardise it
      train_set <- train_set %>%
        mutate(across(3, ~ {
          unique_vals <- unique(.)
          if (all(unique_vals %in% c(0, 1)) && length(unique_vals) == 2) {
            .
          } else if (length(unique_vals) < unique_threshold && all(unique_vals == floor(unique_vals))) {
            .
          } else {
            scale(.)[, 1]
          }
        }))
      # Fill the test set barcode column with "-9" so they're skipped by PLINK in the GWAS
      test_set[, 1] <- -9 
      # Combine train set, test set, and NAs
      combined_set <- rbind(train_set, test_set)
      # Save exposure fold data as text files to all method_type folders
      for (output_dir in output_dirs) {
        dir.create(paste0(output_dir, exposure_name), recursive = T, showWarnings = F)
        gwas_participants_file_name <- paste0(output_dir, exposure_name, "/fold", fold_number, "_gwas_participants.txt")
        write.table(combined_set, file = gwas_participants_file_name, row.names = F, col.names = F, quote = F, sep = "\t", append = F)
      }
      print(paste0("Finished:", exposure_name))
    }
    # Process each fold using parallel processing
    mclapply(seq_along(folds), function(fold_number) process_fold(folds[[fold_number]], fold_number), mc.cores = ncores)
  }
  sapply(exposure_names, split_exposure_i)
}

# Split each exposure
split_exposures(input_df)

writeLines(c(" ", " ", "Finished script 02.1 - Exposures have been split into folds.", " ", " "))

