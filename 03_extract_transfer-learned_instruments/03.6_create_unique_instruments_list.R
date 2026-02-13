# Purpose: List all unique IVs in a text file.

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

# Set the input directory
input <- list.files(path = paste0(results_dir, base_dir, "/instruments"), full.names = T, recursive = T, pattern = "complete_cases_instruments_list.txt")

# Set the output directory for the list of IVs
output <- paste0(results_dir, base_dir, "/instrument_dosages/transfer_learned_instruments/")
dir.create(file.path(output), showWarnings = F, recursive = T)



# Start script---- ---- ---- ---- ---- ----
writeLines(c(" ", " ", "Started script 03.6 - Create unique instruments list.", " ", " "))

snps <- vroom(input, col_names = F, delim = "\n") %>% 
  distinct()

write.table(snps, file = paste0(output, "all_unique_instruments_list.txt"), row.names = F, col.names = F, quote = F, append = F)


writeLines(c(" ", " ", "Finished script 03.6 - Created unique instruments list.", " ", " "))
