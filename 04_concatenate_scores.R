# Purpose: Concatenate test-set and missing-case genetic scores (GSes) into one overall GS per exposure (X)

# Clean environment
rm(list=ls())

# Paths
if (Sys.getenv("GRAPH_MR_ROOT") == "") stop("Set GRAPH_MR_ROOT to the pipeline root directory (see README).")
setwd(Sys.getenv("GRAPH_MR_ROOT"))
source("config/config.R")

# How many cores do you require?
ncores <- 7

# Load the packages required for reading and manipulating the data.
library("tidyverse")
library("vroom")
library("parallel")
library("data.table")

# Specify a base_dir based on the command line inputs
args <- commandArgs(trailingOnly = TRUE)
method_type <- args[1]
input_type <- args[2]
base_dir <- paste0(method_type, "_", input_type)

# Set the input to a list of exposure names
input_exposure_list <- basename(list.dirs(paste0(results_dir, base_dir, "/genetic_scores")))[-1]

# Set the input path
input_path <- paste0(results_dir, base_dir, "/genetic_scores")

# Set the output to the folder containing the genetic scores
output <- paste0(results_dir, base_dir, "/genetic_scores")



# Start script ----
writeLines(c(" ", " ", "Started script 04 - Concatenate scores", " ", " "))


# Create a concatenated PRS vector
join_prses <- function(input_exposure_list) {
  mclapply(input_exposure_list, mc.cores = ncores, FUN = function(exposure_i) {
    cat("Concatenate GSes for:", exposure_i, "\n")
    
    # choose input‐file pattern based on method_type
    file_pattern <- switch(
      method_type,
      cfmr   = "^fold([1-9]|10)_scores\\.csv$",       # fold1_scores.csv … fold10_scores.csv
      cftlmr = "^(fold([1-9]|10)_scores\\.csv|missing-case_scores\\.csv)$",
      stop("Unknown method_type: ", method_type)
    )
    
    # choose output file name based on method_type
    out_fname <- switch(
      method_type,
      cfmr   = "cfmr_scores.csv",
      cftlmr = "cftlmr_scores.csv",
      stop("Unknown method_type: ", method_type)
    )
    
    # List the PRS result files by exposure
    files <- list.files(
      path       = file.path(input_path, exposure_i),
      pattern    = file_pattern,
      full.names = TRUE,
      recursive  = TRUE
    )
    
    # Read and bind the PRS result files by exposure
    joined_gses <- rbindlist(
      l = lapply(files, vroom, delim = ",", show_col_types = FALSE),
      use.names = TRUE,
      fill      = TRUE
    )
    
    # Ensure output directory exists
    dir.create(file.path(output, exposure_i), recursive = TRUE, showWarnings = FALSE)
    
    # Save the bound result file by exposure
    vroom_write(
      joined_gses,
      file  = file.path(output, exposure_i, out_fname),
      delim = ",",
      append = FALSE
    )
    
    cat("Finished concatenating GSes for:", exposure_i, "\n")
  })
}


concatenated_gses <- join_prses(input_exposure_list)

writeLines(c(" ", " ", "Finished script 04 - Concatenated scores.", " ", " "))
