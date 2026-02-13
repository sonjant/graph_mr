# Purpose: Create test-set-specific genetic scores (GSes) for each exposure (X).

# Clean the environment
rm(list=ls())

# Get the task ID from the environment variable
task_id_string <- Sys.getenv("PBS_ARRAY_INDEX")
k <- as.numeric(task_id_string)

# Paths
if (Sys.getenv("GRAPH_MR_ROOT") == "") stop("Set GRAPH_MR_ROOT to the pipeline root directory (see README).")
setwd(Sys.getenv("GRAPH_MR_ROOT"))
source("config/config.R")

# Load packages
library("tidyverse")
library("vroom")
library("parallel")
library("gtools")
library("data.table")

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

# Set the first input to the instrumental variable files, which is the output path from script 02.5 (use primary method_type for reading)
input_ivbetas <- list.files(path = paste0(results_dir, base_dir, "/instruments"), full.names = T, recursive = T, pattern = "fold[0-9]+_instruments\\.tsv$")
# Extract fold numbers and sort numerically to avoid mixedorder warnings with full paths
fold_numbers_ivbetas <- as.numeric(gsub(".*fold([0-9]+)_.*", "\\1", basename(input_ivbetas)))
input_ivbetas <- input_ivbetas[order(fold_numbers_ivbetas)]

# Set the second input to the GWAS participant list files, which is the output path from script 02.1 (use primary method_type for reading)
input_test_ids <- list.files(path = paste0(data_dir, base_dir, "/gwas_input"), full.names = T, recursive = T, pattern = "fold[0-9]+_gwas_participants\\.txt$")
# Extract fold numbers and sort numerically to avoid mixedorder warnings with full paths
fold_numbers_test_ids <- as.numeric(gsub(".*fold([0-9]+)_.*", "\\1", basename(input_test_ids)))
input_test_ids <- input_test_ids[order(fold_numbers_test_ids)]

# Set the third input to the instrument dosages, which is the output from 02.8 (use primary method_type for reading)
input_dosages <- vroom(paste0(ephemeral_dir, "results/", base_dir, "/instrument_dosages/cross_fitted_instruments/all_instrument_dosages.txt")) %>%
  as.data.frame() %>%
  distinct()

# Set output paths for all method types
output_dirs <- sapply(method_types_list, function(mt) {
  dir_path <- paste0(results_dir, mt, "_", input_type, "/genetic_scores/")
  dir.create(dir_path, recursive = T, showWarnings = F)
  return(dir_path)
})



# Start script ---- ---- ---- ---- ---- ----
writeLines(c(" ", " ", "Started script 02.9 - Create cross-fitted genetic scores", " ", " "))

# Modify the column names of input_dosages, excluding the first two columns
clean_name <- function(name) {
  # Remove brackets and their contents
  name <- gsub("\\[.*?\\]", "", name)
  # Remove the last three characters ":DS"
  name <- substr(name, 1, nchar(name) - 3)
  return(name)
}

colnames(input_dosages)[-c(1:4)] <- sapply(colnames(input_dosages)[-c(1:4)], clean_name)

input_dosages <- input_dosages %>% 
  dplyr::rename(CHROM = 1, POS = 2, REF = 3, ALT = 4) %>% 
  tidyr::unite(variant, c(CHROM, POS, REF, ALT), sep = "_", remove = T)

# Create the GS vector for each exposure test fold
create_gs <- function(ivbetas_x_i, test_ids_x_i, input_dosages, standardise = TRUE) {
  # Extract the strings of the exposure name and split number from the file extension
  exposure_name <- strsplit(ivbetas_x_i, "/")[[1]][11]
  fold_i <- strsplit(ivbetas_x_i, "/")[[1]][12] 
  fold_i <- strsplit(fold_i[1], "_")[[1]][1]
  
  # Initialize status message
  status_msg <- ""

  # Read IV betas 
  ivbetas_x_i <- vroom(ivbetas_x_i, delim = "\t", col_types = cols(REF = col_character(), ALT = col_character())) %>%
    dplyr::filter(nchar(REF) == 1 & nchar(ALT) == 1)
  
  # If IV betas is empty, then create a GS vector filled with NAs
  if (nrow(ivbetas_x_i) == 0) {
    status_msg <- "nrow(ivbetas_x_i) == 0 ... "
    test_ids_x_i <- vroom(test_ids_x_i, col_select = 1:2, col_names = c("barcode", "fullbarcode"), col_types = c("cc")) 
    gs <- test_ids_x_i %>% 
      dplyr::filter(barcode == -9) %>% 
      dplyr::select(fullbarcode) %>% 
      dplyr::rename(barcode = fullbarcode) %>% 
      add_column(gs = NA)
  } else {
    ivbetas_x_i <- ivbetas_x_i %>% 
      tidyr::unite(variant, c(CHROM, POS, REF, ALT), sep = "_", remove = TRUE) %>% 
      dplyr::select(variant, BETA) %>% 
      as.data.frame()
    
    # Find overlap between ivbetas_x_i and input_dosages
    overlap_indices <- which(input_dosages$variant %in% ivbetas_x_i$variant)
    
    # If no overlapping variants are found, fallback to NA-filled gs vector
    if (length(overlap_indices) == 0) {
      status_msg <- "No overlapping variants found between ivbetas_x_i and input_dosages, creating NA-filled gs ... "
      test_ids_x_i <- vroom(test_ids_x_i, col_select = 1:2, col_names = c("barcode", "fullbarcode"), col_types = c("cc")) 
      gs <- test_ids_x_i %>% 
        dplyr::filter(barcode == -9) %>% 
        dplyr::select(fullbarcode) %>% 
        dplyr::rename(barcode = fullbarcode) %>% 
        add_column(gs = NA)
    } else {
      test_dosages <- input_dosages[overlap_indices, ]
      
      # Find overlap between test-set participants (test_ids_x_i) and test_dosages participants
      test_ids_x_i_full <- vroom(test_ids_x_i, col_select = 1:2, col_names = c("barcode", "fullbarcode"), col_types = c("cc")) 
      test_ids_x_i <- test_ids_x_i_full %>% 
        dplyr::filter(barcode == -9) %>% 
        dplyr::select(fullbarcode) %>% 
        dplyr::rename(barcode = fullbarcode) %>% 
        pull(barcode)
      
      common_barcodes <- intersect(test_ids_x_i, colnames(test_dosages)) 
      common_barcodes <- c("variant", common_barcodes)
      test_dosages <- test_dosages %>% 
        dplyr::select(all_of(common_barcodes))
      
      # Find overlap between ivbetas_x_i with test_dosages
      test_dosages_for_ivs <- inner_join(ivbetas_x_i, test_dosages, by = "variant")
      dosages <- test_dosages_for_ivs %>% 
        dplyr::select(c(-variant, -BETA)) %>% 
        t() 
      dosages <- apply(dosages, 2, as.numeric)
      dosages <- as.matrix(dosages)
      ivbetas <- as.matrix(as.numeric(test_dosages_for_ivs$BETA))
      
      # If all entries in gs$gs are 0, then replace them with NA
      # Else if the number of columns in the matrix and the number of rows in the ivbetas_x_i vector match, create the gs vector
      if (length(ivbetas) == 0) {
        status_msg <- "All entries in gs vector are 0, creating NA-filled gs ... "
        gs <- test_ids_x_i_full %>% 
          add_column(gs = NA) %>% 
          dplyr::select(barcode, gs)
      } else if (ncol(dosages) == length(ivbetas)) { 
        status_msg <- "Dosages * Betas ... "
        gs <- dosages %*% ivbetas
        barcode <- colnames(test_dosages_for_ivs)[-(1:2)]
        gs <- data.frame(barcode = barcode, gs = gs) 
        if (standardise) {
          gs$gs <- scale(gs$gs)
        }
      } else {
        # If dimensions don't match, create NA-filled gs to ensure output is always written
        status_msg <- "Dimensions don't match, creating NA-filled gs ... "
        gs <- test_ids_x_i_full %>% 
          add_column(gs = NA) %>% 
          dplyr::select(barcode, gs)
      }
    }
  }
  
  # Create result directories for this exposure's genetic risk scores in all method_type folders
  for (output_dir in output_dirs) {
    dir.create(paste0(output_dir, exposure_name), showWarnings = FALSE, recursive = TRUE)
  }

  # Save the gs vector to all method_type folders
  for (output_dir in output_dirs) {
    vroom_write(gs, file = paste0(output_dir, exposure_name, "/", fold_i, "_scores.csv"), delim = ",", append = FALSE)
  }
  cat(status_msg, "Genetic score for", exposure_name, fold_i, " created!\n")
}

    

# Set the chunk size to your k (the number of folds) 
chunk_size <- 10

# Set the number and size of chunks
num_chunks_per_subjob <- 10  
total_chunks <- ceiling(length(input_ivbetas) / chunk_size)
num_subjobs <- ceiling(total_chunks / num_chunks_per_subjob)

# Determine the start and end chunk index for the current subjob
start_chunk_index <- (k - 1) * num_chunks_per_subjob + 1
end_chunk_index <- min(k * num_chunks_per_subjob, total_chunks)

# Loop over the chunks for the current subjob
for (chunk_index in start_chunk_index:end_chunk_index) {
  start_index <- (chunk_index - 1) * chunk_size + 1
  end_index <- min(chunk_index * chunk_size, length(input_ivbetas))
  ivbetas_x_i <- input_ivbetas[start_index:end_index]
  test_ids_x_i <- input_test_ids[start_index:end_index]
  gs <- mapply(create_gs, ivbetas_x_i, test_ids_x_i, MoreArgs = list(input_dosages = input_dosages))
}

cat("\n","Finished script 02.9 - Cross-fitted genetic scores created.")
