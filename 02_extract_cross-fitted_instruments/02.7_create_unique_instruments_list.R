# Purpose: List all unique instruments in a text file.

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

# Set the input directory (use primary method_type for reading)
input <- list.files(path = paste0(results_dir, base_dir, "/instruments"), full.names = T, recursive = T, pattern = "fold[0-9]+_instruments_list\\.txt")

# Set output directories for all method types
output_dirs <- sapply(method_types_list, function(mt) {
  dir_path <- paste0(results_dir, mt, "_", input_type, "/instrument_dosages/cross_fitted_instruments/")
  dir.create(file.path(dir_path), showWarnings = F, recursive = T)
  return(dir_path)
})



# Start script---- ---- ---- ---- ---- ----
writeLines(c(" ", " ", "Started script 02.7 - Create list of unique instruments.", " ", " "))

snps <- vroom(input, col_names = F, delim = "\n") %>% 
  distinct()

# Write to all method_type folders
for (output_dir in output_dirs) {
  write.table(snps, file = paste0(output_dir, "all_unique_instruments_list.txt"), row.names = F, col.names = F, quote = F, append = F)
}


writeLines(c(" ", " ", "Finished script 02.7 - Created list of unique instruments.", " ", " "))
