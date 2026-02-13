# Format Graph MR results into a filtered matrix to plot a network

# Clean environment
rm(list=ls())

# Paths
if (Sys.getenv("GRAPH_MR_ROOT") == "") stop("Set GRAPH_MR_ROOT to the pipeline root directory (see README).")
setwd(Sys.getenv("GRAPH_MR_ROOT"))
source("config/config.R")

# Load packages 
library("tidyverse")
library("vroom") 

# Specify a base_dir based on the command line inputs
args <- commandArgs(trailingOnly = TRUE)
method_type <- args[1]
input_type <- args[2]
base_dir <- paste0(method_type, "_", input_type)

# Set the inputs to the F-statistics and P-values from the X->Y MR regressions (unadjusted and adjusted)
# Unadjusted models
input_fstat_unadjusted <- vroom(paste0(results_dir, base_dir, "/matrices_chem_names/totalfx_unadjusted_F_stat.csv"), .name_repair = "minimal", show_col_types = F)
input_pval_unadjusted <- vroom(paste0(results_dir, base_dir, "/matrices_chem_names/totalfx_unadjusted_P_value.csv"), .name_repair = "minimal", show_col_types = F)

# Adjusted models
input_fstat_adjusted <- vroom(paste0(results_dir, base_dir, "/matrices_chem_names/totalfx_agesex_F_stat.csv"), .name_repair = "minimal", show_col_types = F)
input_pval_adjusted <- vroom(paste0(results_dir, base_dir, "/matrices_chem_names/totalfx_agesex_P_value.csv"), .name_repair = "minimal", show_col_types = F)

# Set and create the output file path
output_path <- paste0(results_dir, base_dir, "/matrices_chem_names/")
dir.create(output_path, recursive = T, showWarnings = F)


# Remove second and later instances of duplicate exposures in columns
remove_second_and_later_duplicates <- function(df) {
  exposures <- df[[1]]
  colnames_exposures <- colnames(df)[-1]
  row_duplicates <- duplicated(exposures)
  col_duplicates <- duplicated(colnames_exposures)
  df_unique_rows <- df[!row_duplicates, ]
  df_unique <- df_unique_rows[, c(TRUE, !col_duplicates)]
  return(df_unique)
}

# Process unadjusted data
input_fstat_unadjusted <- remove_second_and_later_duplicates(input_fstat_unadjusted)
input_pval_unadjusted <- remove_second_and_later_duplicates(input_pval_unadjusted)

# Process adjusted data
input_fstat_adjusted <- remove_second_and_later_duplicates(input_fstat_adjusted)
input_pval_adjusted <- remove_second_and_later_duplicates(input_pval_adjusted)

# Start script ----
writeLines(c(" ", " ", "Started script 07.1 - Format data for graphs.", " ", " "))

# Create TL-CFMR adjacency matrix where exposure-GSes have F-statistic > 10 and Bonferroni-significant P-values for total-effects)
create_adjmatrix <- function(input_fstat, input_pval, suffix = "") {
  # Filter for instruments with F-statistic > 10
  fstat.matrix.data <- input_fstat[, -1]
  fstat.matrix.data <- as.matrix(fstat.matrix.data)
  fstat.matrix.data <- diag(fstat.matrix.data)
  fstats <- data.frame(input_fstat[, 1], fstat.matrix.data)
  colnames(fstats) <- c("Exposure", "Fstatistic")
  fstats <- fstats[order(-fstats$Fstatistic), ]
  fstats <- fstats[(fstats$Fstatistic > 10), ]
  # Create a p-value adjacency matrix with GSes where F-statistic > 10
  pvalues <- input_pval %>% 
    dplyr::inner_join(y = fstats, by = "Exposure") %>% 
    dplyr::filter(!is.na(Exposure)) %>%    
    dplyr::distinct(Exposure, .keep_all = TRUE) %>%  
    tibble::column_to_rownames("Exposure") %>%         
    dplyr::select(-Fstatistic) %>% 
    as.data.frame()
  # Remove self-loops
  for (name in rownames(pvalues)) {  
    if (name %in% colnames(pvalues)) {
      pvalues[name, name] <- NA
    }
  }
  # Function to filter P-values at nominal significance
  nominal_significance <- function(pvalues, threshold = 0.05) {
    adjusted <- pvalues
    adjusted[adjusted >= threshold] <- NA
    return(as.data.frame(adjusted))
  }
  # Function to filter P-values at any given threshold
  adjust_pvalues <- function(pvalues, method, n_tests) {
    adjusted <- apply(pvalues, 2, function(x) p.adjust(x, method = method, n = n_tests))
    adjusted <- as.data.frame(adjusted)
    adjusted[adjusted >= 0.05] <- NA
    return(adjusted)
  }
  # Define the number of tests for each correction
  n_tests <- nrow(pvalues) * ncol(pvalues) - nrow(pvalues)
  # Perform p-value adjustments
  pvalues_nominal <- nominal_significance(pvalues, threshold = 0.05)
  pvalues_bonferroni <- adjust_pvalues(pvalues, method = "bonferroni", n_tests = n_tests)
  # Save row names as Exposure
  pvalues_nominal <- pvalues_nominal %>% rownames_to_column("Exposure")
  pvalues_bonferroni <- pvalues_bonferroni %>% rownames_to_column("Exposure")
  # Save the adjacency matrices with appropriate suffix
  if (suffix == "") {
    write.table(pvalues_nominal, file = paste0(output_path, "fstat10-nominal-totalfx-matrix.csv"), row.names = F, sep = ",") 
    write.table(pvalues_bonferroni, file = paste0(output_path, "fstat10-bonferroni-totalfx-matrix.csv"), row.names = F, sep = ",") 
  } else {
    write.table(pvalues_nominal, file = paste0(output_path, "fstat10-nominal-totalfx-", suffix, "-matrix.csv"), row.names = F, sep = ",") 
    write.table(pvalues_bonferroni, file = paste0(output_path, "fstat10-bonferroni-totalfx-", suffix, "-matrix.csv"), row.names = F, sep = ",") 
  }
}

# Process unadjusted models
writeLines("Processing unadjusted models...")
adjmatrix_unadjusted <- create_adjmatrix(input_fstat_unadjusted, input_pval_unadjusted, suffix = "unadjusted")

# Process adjusted models
writeLines("Processing adjusted models...")
adjmatrix_adjusted <- create_adjmatrix(input_fstat_adjusted, input_pval_adjusted, suffix = "agesex")

writeLines(c(" ", " ", "Finished script 07.1 - Formatted data for graphs.", " ", " "))
