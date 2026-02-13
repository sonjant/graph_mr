# Purpose: Collate the characteristics of out instrumental variables and genetic scores into summary tables

# Clean environment
rm(list=ls())

# Paths
if (Sys.getenv("GRAPH_MR_ROOT") == "") stop("Set GRAPH_MR_ROOT to the pipeline root directory (see README).")
setwd(Sys.getenv("GRAPH_MR_ROOT"))
source("config/config.R")

# Load packages
library("tidyverse")
library("vroom")
library("data.table")
library("readxl")
library("parallel")

# Specify a base_dir based on the command line inputs
args <- commandArgs(trailingOnly = TRUE)
method_type <- args[1]
input_type <- args[2]
base_dir <- paste0(method_type, "_", input_type)

# Specify number of cores
ncores <- 4

# Set the first input to the IV dimension files
input_dims <- paste0(results_dir, base_dir, "/instruments/")

# Set the inputs to the Graph MR total effect results (unadjusted and adjusted)
input_fx_unadjusted <- paste0(results_dir, base_dir, "/total_fx/")
input_fx_adjusted <- paste0(results_dir, base_dir, "/total_fx/")

# Set the output for the summary of IV characteristics (same as the second output in 04.1)
output <- paste0(results_dir, base_dir, "/summary_stats/")
dir.create(output, recursive = T, showWarnings = F)



# Start script ----
writeLines(c(" ", " ", "Started script 06.2 - Collate summary table", " ", " "))

# Process files in batches 
process_files_in_batches <- function(input_path, pattern, batch_size) {
  # List files in batches
  list_files_in_batches <- function(path, pattern, batch_size) {
    all_files <- dir(path, pattern = pattern, recursive = TRUE, full.names = TRUE)
    file_batches <- split(all_files, ceiling(seq_along(all_files) / batch_size))
    return(file_batches)
  }
  dim_stats_batches <- list_files_in_batches(path = input_path, pattern = pattern, batch_size = batch_size)
  combined_df_list <- list()
  # Process each batch
  for (batch in dim_stats_batches) {
    combined_batch_df <- vroom(batch, col_types = cols(), show_col_types = F) 
    combined_df_list <- c(combined_df_list, list(combined_batch_df))
  }
  # Combine all dataframes into one
  combined_df <- bind_rows(combined_df_list)
  return(combined_df)
}

# Summarise SNP dimension stats for each exposure (X) (output from 02.5)
summarised_data <- process_files_in_batches(input_path = input_dims, pattern = "fold[0-9]+_instrument_dimensions\\.csv$", batch_size = 500)

summarised_data <- summarised_data %>%
  dplyr::group_by(Exposure) %>%
  dplyr::summarise(
    mean_n_gwas_snps = mean(n_gwas_snps, na.rm = TRUE),
    sd_n_gwas_snps = sd(n_gwas_snps, na.rm = TRUE),
    median_n_gwas_snps = median(n_gwas_snps, na.rm = TRUE),
    mean_n_clumped_snps = mean(n_clumped_snps, na.rm = TRUE))

# Read the n_empty_folds file and group the stats by Exposure (output from 02.6)
emptyfolds <- vroom(paste0(output, "n_empty_folds.csv"), delim = ",", show_col_types = F) %>% 
  dplyr::mutate("n_folds_containing_instruments" = 10 - N.empty_rows) %>% 
  dplyr::select(-N.empty_rows)

summarised_data <- dplyr::inner_join(summarised_data, emptyfolds, by = "Exposure") 

# Process complete cases instrument dimensions if method_type is 'cftlmr'
if (method_type == 'cftlmr') {
  # Check if complete cases files exist
  cc_files <- dir(input_dims, pattern = "complete_cases_instruments_dims\\.csv$", recursive = TRUE, full.names = TRUE)
  if (length(cc_files) > 0) {
    # Process complete cases files
    complete_cases_summarised_data <- process_files_in_batches(input_path = input_dims, pattern = "complete_cases_instruments_dims\\.csv$", batch_size = 500)
    
    complete_cases_summarised_data <- complete_cases_summarised_data %>%
      dplyr::group_by(Exposure) %>%
      dplyr::summarise(
        cc_n_gwas_snps = mean(n_gwas_snps, na.rm = TRUE),
        cc_n_clumped_snps = mean(n_clumped_snps, na.rm = TRUE))
    
    # Join complete cases stats to summarised_data
    summarised_data <- summarised_data %>%
      dplyr::left_join(complete_cases_summarised_data, by = "Exposure")
  }
}

# Add F-stat, R-squared, and P-value to the GS summary table (unadjusted and adjusted)
# Process unadjusted files
stats.files_unadjusted <- process_files_in_batches(input_path = input_fx_unadjusted, pattern = "_total_fx\\.csv$", batch_size = 500)

stats_unadjusted <- stats.files_unadjusted %>%
  dplyr::mutate(Exposure = str_replace_all(Exposure, "_scores", "")) %>%
  dplyr::filter(Exposure == Response) %>%
  dplyr::select(Exposure, R_squared, P_value, F_stat, Causal_estimate) %>%
  dplyr::rename(
    R_squared_unadjusted = R_squared,
    P_value_unadjusted = P_value,
    F_stat_unadjusted = F_stat,
    Causal_estimate_unadjusted = Causal_estimate
  )

# Process adjusted files
stats.files_adjusted <- process_files_in_batches(input_path = input_fx_adjusted, pattern = "_total_fx_agesex\\.csv$", batch_size = 500)

stats_adjusted <- stats.files_adjusted %>%
  dplyr::mutate(Exposure = str_replace_all(Exposure, "_scores", "")) %>%
  dplyr::filter(Exposure == Response) %>%
  dplyr::select(Exposure, R_squared, P_value, F_stat, Causal_estimate) %>%
  dplyr::rename(
    R_squared_adjusted = R_squared,
    P_value_adjusted = P_value,
    F_stat_adjusted = F_stat,
    Causal_estimate_adjusted = Causal_estimate
  )

# Join both unadjusted and adjusted stats to summarised_data
summarised_data <- summarised_data %>%
  dplyr::left_join(stats_unadjusted, by = "Exposure") %>%
  dplyr::left_join(stats_adjusted, by = "Exposure") 

# Read the chemical metadata
metabolomics_metadata <- vroom(paste0(home, "data/airwave/", base_dir, "/preprocessed/metabolomics_chemical_metadata.tsv"), show_col_types = F) %>% 
  dplyr::select(HMDB_ID, BIOCHEMICAL) %>% 
  dplyr::rename(Exposure = HMDB_ID)

proteomics_metadata <- vroom(paste0(home, "data/airwave/", base_dir, "/preprocessed/proteomics_chemical_metadata.tsv"), show_col_types = F) %>% 
  dplyr::select(UniProt, protein_name) %>%
  dplyr::rename(Exposure = UniProt)

combined_metadata <- dplyr::bind_rows(metabolomics_metadata %>% mutate(Source = "Metabolite"),
                                      proteomics_metadata %>% mutate(Source = "Protein"))

# Create mapping from Exposure to BIOCHEMICAL/protein_name
# For metabolites: use BIOCHEMICAL, for proteins: use protein_name
# If neither BIOCHEMICAL nor protein_name is available, keep the Exposure as is
exposure_to_name_mapping <- combined_metadata %>%
  dplyr::mutate(
    Chemical = dplyr::case_when(
      !is.na(BIOCHEMICAL) ~ BIOCHEMICAL,
      !is.na(protein_name) ~ protein_name,
      TRUE ~ Exposure
    )
  ) %>%
  dplyr::select(Exposure, Chemical) %>%
  dplyr::distinct()

# Join the mapping to summarised_data and replace Exposure with Chemical
summarised_data_named <- summarised_data %>%
  dplyr::left_join(exposure_to_name_mapping, by = "Exposure") %>%
  dplyr::mutate(Exposure = dplyr::coalesce(Chemical, Exposure)) %>%
  dplyr::select(-Chemical) 

# Save the summary attributes of each Genetic score
write.table(summarised_data_named, file = paste0(output, "genetics_scores_attributes_summary.csv"), sep = ",", row.names = F)



#Collate IV dimensions at different stages of the pipeline 
list_files_in_batches <- function(path, pattern, batch_size) {
  all_files <- dir(path, pattern = pattern, recursive = TRUE, full.names = TRUE)
  file_batches <- split(all_files, ceiling(seq_along(all_files) / batch_size))
  return(file_batches)
}

# List all fold-specific SNP dimension files
dim_stats <- list_files_in_batches(path = paste0(home, "results/", base_dir, "/instruments"), pattern = "fold[0-9]+_instrument_dimensions\\.csv$", batch_size = 500)
dim_stats <- do.call(c, dim_stats)
dim_stats <- vroom(dim_stats, show_col_types = F)
dim_stats <- dim_stats %>% 
  dplyr::left_join(exposure_to_name_mapping, by = "Exposure") %>%
  dplyr::mutate(Exposure = dplyr::coalesce(Chemical, Exposure)) %>%
  dplyr::select(-Chemical) 

write.table(dim_stats, file = paste0(output, "instrumental_variables_dimensions.csv"), sep = ",", row.names = F)

# List all complete-case SNP dimension files (only if method_type is 'cftlmr') 
if (method_type == 'cftlmr') {
  cc_dim_stats <- list_files_in_batches(path = paste0(home, "results/", base_dir, "/instruments"), pattern = "complete_cases_instruments_dims.csv", batch_size = 500)
  cc_dim_stats <- do.call(c, cc_dim_stats)
  if (length(cc_dim_stats) > 0) {
    cc_dim_stats <- vroom(cc_dim_stats, show_col_types = F)
    cc_dim_stats <- cc_dim_stats %>% 
      dplyr::left_join(exposure_to_name_mapping, by = "Exposure") %>%
      dplyr::mutate(Exposure = dplyr::coalesce(Chemical, Exposure)) %>%
      dplyr::select(-Chemical) 
    
    write.table(cc_dim_stats, file = paste0(output, "cc_instrumental_variables_dimensions.csv"), sep = ",", row.names = F)
  }
}



writeLines(c(" ", " ", "Finished script 06.2 - We have collated a summary table", " ", " "))
