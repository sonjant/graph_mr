# Purpose: for each exposure (X), create transfer-learned scores for the missing cases 

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

# Set the first input to the instrumental variable files, which is the "output" path of script 03.4
input_ivbetas <- list.files(path = paste0(results_dir, base_dir, "/instruments"), full.names = T, recursive = T, pattern = "complete_cases_instruments.tsv")
input_ivbetas <- input_ivbetas[mixedorder(input_ivbetas)]

# Set the second input to the GWAS participant list files, which is the output path from script 03.1
input_missing_x_ids <- list.files(path = paste0(data_dir, base_dir, "/gwas_input"), full.names = T, recursive = T, pattern = "participants_missing_exposure.txt")
input_missing_x_ids <- input_missing_x_ids[mixedorder(input_missing_x_ids)]

# Set the third input to the IV dosages, which is the output of the script 03.7
input_dosages <- vroom(paste0(ephemeral_dir, "results/", base_dir, "/instrument_dosages/transfer_learned_instruments/all_instrument_dosages.txt")) %>% 
  as.data.frame() %>% 
  distinct()

# Set the output path for the genetic scores
output <- paste0(results_dir, base_dir, "/genetic_scores/")



# Start script ---- ---- ---- ---- ---- ----
writeLines(c(" ", " ", "Started script 03.8 - Create missing-case genetic scores (GSes)", " ", " "))

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

# Create the GS vector for each participant with missing exposure cases
create_gs <- function(ivbetas_x_i, missing_x_i_ids, input_dosages, standardise = TRUE) {
  # Extract the strings of the exposure name and split number from the file extension
  exposure_name <- strsplit(ivbetas_x_i, "/")[[1]][11]
  
  # Initialize status message
  status_msg <- ""
  
  # Read IV betas 
  ivbetas_x_i <- vroom(ivbetas_x_i, delim = "\t", col_types = cols(REF = col_character(), ALT = col_character())) %>%
    dplyr::filter(nchar(REF) == 1 & nchar(ALT) == 1) 
  
  # Read IDs of participants missing exposure - if file is empty, create empty output file with NA, else continue
  missing_x_i_ids <- suppressMessages(suppressWarnings(vroom(missing_x_i_ids, col_names = "barcode", delim = "\n", col_types = "c"))) %>% 
    as.data.frame()
  
  if (nrow(missing_x_i_ids) == 0) {
    # Create empty file with NA to ensure output is always written
    gs <- data.frame(barcode = character(0), gs = numeric(0))
  } else {
    # If IV betas is empty, then create a GS vector filled with NAs
    if (nrow(ivbetas_x_i) == 0) {
      gs <- missing_x_i_ids %>%
        add_column(gs = NA)
    } else {
      ivbetas_x_i <- ivbetas_x_i %>% 
        tidyr::unite(variant, c(CHROM, POS, REF, ALT), sep = "_", remove = TRUE) %>% 
        dplyr::select(variant, BETA) %>% 
        as.data.frame()
      
      # Find overlap between ivbetas_x_i and input_dosages
      overlap_indices <- which(input_dosages$variant %in% ivbetas_x_i$variant)
      
      # If no overlapping variants, return NA-filled GSs
      if (length(overlap_indices) == 0) {
        gs <- missing_x_i_ids %>%
          add_column(gs = NA)
      } else {
        overlap_dosages <- input_dosages[overlap_indices, ]
        
        # Find overlap between missing-case participants and overlap_dosages participants
        common_barcodes <- intersect(missing_x_i_ids$barcode, colnames(overlap_dosages)) 
        common_barcodes <- c("variant", common_barcodes)
        missing_x_i_dosages <- overlap_dosages %>% 
          dplyr::select(all_of(common_barcodes))
        
        # If none found, create NA-filled GSs
        if (nrow(overlap_dosages) == 0) {
          gs <- missing_x_i_ids %>%
            add_column(gs = NA)
        } else {
          # Find overlap between ivbetas_x_i with test_input_dosages
          missing_x_i_dosages_for_ivs <- inner_join(ivbetas_x_i, missing_x_i_dosages, by = "variant")
          
          dosages <- missing_x_i_dosages_for_ivs %>% 
            dplyr::select(c(-variant, -BETA)) %>% 
            t() 
          dosages <- apply(dosages, 2, as.numeric)
          dosages <- as.matrix(dosages)
          ivbetas <- as.matrix(as.numeric(missing_x_i_dosages_for_ivs$BETA))
          
          # If all entries in gs$gs are 0, then replace them with NA
          # Else if the number of columns in the matrix and the number of rows in the ivbetas_x_i vector match, create the GS vector
          if (length(ivbetas) == 0) {
            status_msg <- "All entries in GS vector are 0, creating NA-filled gs ... "
            gs <- missing_x_i_ids %>% 
              add_column(gs = NA) %>% 
              dplyr::select(barcode, gs)
          } else if (ncol(dosages) == length(ivbetas)) { 
            status_msg <- "Multi-variant case ... "
            gs <- dosages %*% ivbetas
            barcode <- colnames(missing_x_i_dosages_for_ivs)[-(1:2)]
            gs <- data.frame(barcode = barcode, gs = gs) 
            if (standardise) {
              gs$gs <- scale(gs$gs)
            } 
          } else if (length(dosages) == length(ivbetas)) {
            status_msg <- "Single-variant case ... "
            gs <- sum(dosages * ivbetas)
            barcode <- colnames(missing_x_i_dosages_for_ivs)[-(1:2)]
            gs <- data.frame(barcode = barcode, gs = gs) 
          } else {
            # If dimensions don't match, create NA-filled gs to ensure output is always written
            status_msg <- "Dimensions don't match, creating NA-filled gs ... "
            gs <- missing_x_i_ids %>% 
              add_column(gs = NA) %>% 
              dplyr::select(barcode, gs)
          }
        }
      }
    }
  }
  
  # Create result directory for this exposure's genetic risk scores
  dir.create(paste0(output, exposure_name), showWarnings = FALSE, recursive = TRUE)
  
  # Save the gs vector
  vroom_write(gs, file = paste0(output, exposure_name, "/missing-case_scores.csv"), delim = ",", append = FALSE)
  cat(status_msg, "Missing-case genetic scores for", exposure_name, " created!\n")
}


# Set the number of elements each subjob should process
elements_per_subjob <- 50

# Calculate the total number of subjobs needed
total_jobs <- ceiling(length(input_ivbetas) / elements_per_subjob)

# Determine the start and end index for the current job
start_index <- (k - 1) * elements_per_subjob + 1
end_index <- min(k * elements_per_subjob, length(input_ivbetas))

# Create the chunks for the current job
ivbetas_x_i_chunk <- input_ivbetas[start_index:end_index]
missing_x_i_ids_chunk <- input_missing_x_ids[start_index:end_index]

# Process each element pair individually
for (i in 1:length(ivbetas_x_i_chunk)) {
  ivbetas_x_i <- ivbetas_x_i_chunk[i]
  missing_x_i_ids <- missing_x_i_ids_chunk[i]
  create_gs(ivbetas_x_i, missing_x_i_ids, input_dosages)
}

cat("\n","Finished script 03.8 - Missing-case genetic scores created.")
