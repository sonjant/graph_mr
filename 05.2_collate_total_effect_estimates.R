# Purpose: Collate Graph MR total effect estimates.

# Clean environment
rm(list=ls())

# Paths
if (Sys.getenv("GRAPH_MR_ROOT") == "") stop("Set GRAPH_MR_ROOT to the pipeline root directory (see README).")
setwd(Sys.getenv("GRAPH_MR_ROOT"))
source("config/config.R")

# How many cores do you need?
ncores <- 5

# Load packages 
library("tidyverse")
library("vroom") 
library("parallel")
library("data.table")

# Specify a base_dir based on the command line inputs
args <- commandArgs(trailingOnly = TRUE)
method_type <- args[1]
input_type <- args[2]
base_dir <- paste0(method_type, "_", input_type)

# Set the input directory for total effect result files
total_fx_dir <- paste0(results_dir, base_dir, "/total_fx")

# Get unadjusted files (*_total_fx.csv)
totalfx_unadjusted_files <- list.files(total_fx_dir, pattern = "_total_fx\\.csv$", full.names = T, recursive = T)
# Get adjusted files (*_total_fx_agesex.csv)
totalfx_adjusted_files <- list.files(total_fx_dir, pattern = "_total_fx_agesex\\.csv$", full.names = T, recursive = T)

# Read and combine unadjusted files
if (length(totalfx_unadjusted_files) > 0) {
  totalfx_unadjusted <- mclapply(totalfx_unadjusted_files, fread, mc.cores = ncores) 
  totalfx_unadjusted <- rbindlist(totalfx_unadjusted)
} else {
  warning("No unadjusted total_fx files found")
  totalfx_unadjusted <- data.table()
}

# Read and combine adjusted files
if (length(totalfx_adjusted_files) > 0) {
  totalfx_adjusted <- mclapply(totalfx_adjusted_files, fread, mc.cores = ncores) 
  totalfx_adjusted <- rbindlist(totalfx_adjusted)
} else {
  warning("No adjusted total_fx_agesex files found")
  totalfx_adjusted <- data.table()
}

# Set the output file path
output <- paste0(results_dir, base_dir, "/matrices/")
dir.create(output, recursive = T, showWarnings = F)



# Start script ----
writeLines(c(" ", " ", "Started script 05.2 - Graph MR total effect estimate collation.", " ", " "))

# Create R-squared, adjusted R-squared, causal effect estimate, P-value, and F-stat matrices
collate_totalfx_estimates <- function(totalfx, suffix = ""){
  # Define the list of stat_var arguments inside the function
  stat_vars <- c("R_squared", "Causal_estimate", "P_value", "F_stat")
  # Use mclapply inside the function
  mclapply(stat_vars, mc.cores = ncores, FUN = function(stat_var) {
    # Apply dcast.data.table with current stat_var
    dt_cast <- dcast.data.table(totalfx, Exposure ~ Response, value.var = stat_var)
    # Save the resulting data.table to a .csv file using vroom_write
    output_filename <- if (suffix == "") {
      paste0(output, "totalfx_", stat_var, ".csv")
    } else {
      paste0(output, "totalfx_", suffix, "_", stat_var, ".csv")
    }
    vroom_write(dt_cast, output_filename, delim = ",")
    }
  )
}

# Apply the function to unadjusted and adjusted data.tables
if (nrow(totalfx_unadjusted) > 0) {
  writeLines("Collating unadjusted results...")
  collated_results_unadjusted <- collate_totalfx_estimates(totalfx_unadjusted, suffix = "unadjusted")
}

if (nrow(totalfx_adjusted) > 0) {
  writeLines("Collating adjusted results...")
  collated_results_adjusted <- collate_totalfx_estimates(totalfx_adjusted, suffix = "agesex")
}

writeLines(c(" ", " ", "Finished script 05.2 - finished collating Graph MR total effect estimates.", " ", " "))
