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

# Specify a base_dir based on the command line inputs
args <- commandArgs(trailingOnly = TRUE)
method_type <- args[1]
input_type <- args[2]
base_dir <- paste0(method_type, "_", input_type)

# Set the number of cores you require
ncores <- 6

# Set the input_files path to the "output" path from script 03.3
input_list <- list.dirs(path = paste0(ephemeral_dir, "results/", base_dir, "/gwas_results"))
input_list <- grep("complete_cases$", input_list, value = T)

# Set an output file path for the selected instrumental variables
output <- paste0(results_dir, base_dir, "/instruments/")



# Start script ---- ---- ---- ---- ---- ----
writeLines(c(" ", " ", "Started script 03.4 - PLINK data loading and cleaning", " ", " "))

# Bind and save the GWAS result files for each exposure
save_gwas <- function(input_list){
  # For each exposure, bind and save the GWAS result files
  mclapply(input_list, mc.cores = ncores, FUN = function(input_i){
    # Extract the name of the input exposure
    input_name <- basename(dirname(input_i))
    # Create output directory
    dir.create(paste0(output, input_name), showWarnings = F, recursive = T)
    # List, read, and bind the 22 GWAS result files for the current exposure
    gwas_res_filelist <- list.files(path = input_i, pattern = "*glm.linear", full.names = T, recursive = T)
    gwas_res <- rbindlist(l = lapply(gwas_res_filelist, fread, 
                                     select = c("ID", "#CHROM", "POS", "BETA", "SE", "P", "REF", "ALT"),
                                     col.names = c("ID", "CHROM", "POS", "BETA", "SE", "P", "REF", "ALT")), 
                          use.names = T, fill = T) 
    # Filter out rows where the column "ID" contains the value "ID"
    gwas_res <- gwas_res[gwas_res$ID != "ID", ]
    # If GWAS result is empty, create a dataframe, if not, proceed with binding and saving the results
    if (nrow(gwas_res) == 0) {
      # Create an empty dataframe with specified columns and save it as a file so that the number of files remains equal to the # of exposures)
      gwas_res <- data.frame(ID = character(), 
                             CHROM = character(),
                             POS = integer(),
                             BETA = numeric(),
                             SE = numeric(),
                             P = numeric(),
                             REF = character(),
                             ALT = character(),
                             stringsAsFactors = FALSE)
      vroom_write(gwas_res, col_names = T, file = paste0(output, "/", input_name,  "/complete_cases_instruments.tsv"))
      # Also create a clumped_res object with 0 rows
      clumped_res <- data.frame(SNP = character())
      } else {
      # List, read, bind, and save the clumped results for the current exposure
      clumped_res_filelist <- list.files(path = input_i, pattern = "*.clumped.clumps", full.names = T, recursive = T)
      clumped_res <- rbindlist(l = lapply(clumped_res_filelist, fread, select = c("ID")), use.name = T, fill = T)
      clumped_res <- clumped_res %>%
        inner_join(y = gwas_res, by = "ID")
      vroom_write(clumped_res, col_names = T, file = paste0(output, "/", input_name, "/complete_cases_instruments.tsv"))
      }
    # Save the number of SNPs selected by the GWAS and after clumping in a dataframe
    gwas_snps_dim <- data.frame(Exposure = input_name,
                                "n_gwas_snps" = nrow(gwas_res),
                                "n_clumped_snps" = nrow(clumped_res),
                                check.names = F)
    vroom_write(gwas_snps_dim, col_names = T, file = paste0(output, "/", input_name, "/complete_cases_instruments_dims.csv"))
    # Print the current exposure
    writeLines(c(" ", " ", "Joined PLINK results for:", input_i, " ", " "))
    }
  )
}

saved_gwas_res <- save_gwas(input_list)

writeLines(c(" ", " ", "Finished script 03.4 - finished PLINK data loading and cleaning", " ", " "))
