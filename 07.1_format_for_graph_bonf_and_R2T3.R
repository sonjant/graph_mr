# Format Graph MR results into a filtered matrix to plot a network (Bonferroni and R² top-3)

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
input_r2_unadjusted <- vroom(paste0(results_dir, base_dir, "/matrices_chem_names/totalfx_unadjusted_R_squared.csv"), .name_repair = "minimal", show_col_types = F)
input_summary_unadjusted <- vroom(paste0(results_dir, base_dir, "/summary_stats/genetics_scores_attributes_summary.csv"), .name_repair = "minimal", show_col_types = F)

# Adjusted models
input_fstat_adjusted <- vroom(paste0(results_dir, base_dir, "/matrices_chem_names/totalfx_agesex_F_stat.csv"), .name_repair = "minimal", show_col_types = F)
input_pval_adjusted <- vroom(paste0(results_dir, base_dir, "/matrices_chem_names/totalfx_agesex_P_value.csv"), .name_repair = "minimal", show_col_types = F)
input_r2_adjusted <- vroom(paste0(results_dir, base_dir, "/matrices_chem_names/totalfx_agesex_R_squared.csv"), .name_repair = "minimal", show_col_types = F)
input_summary_adjusted <- vroom(paste0(results_dir, base_dir, "/summary_stats/genetics_scores_attributes_summary.csv"), .name_repair = "minimal", show_col_types = F)

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
input_r2_unadjusted <- remove_second_and_later_duplicates(input_r2_unadjusted)

# Process adjusted data
input_fstat_adjusted <- remove_second_and_later_duplicates(input_fstat_adjusted)
input_pval_adjusted <- remove_second_and_later_duplicates(input_pval_adjusted)
input_r2_adjusted <- remove_second_and_later_duplicates(input_r2_adjusted)

# Start script ----
writeLines(c(" ", " ", "Started script 07.1 - Format data for graphs with R² rank top 3 filtering.", " ", " "))

# Compute R² ranks for cross-trait predictive ranking
compute_r2_ranks <- function(r_squared_df, summary_df, fstat_col) {
  # Filter to exposures with F-stat > 10
  idx <- which(summary_df[[fstat_col]] > 10 & !is.na(summary_df[[fstat_col]]))
  exposures_f10 <- summary_df$Exposure[idx]
  
  r_squared_filtered <- r_squared_df[r_squared_df$Exposure %in% exposures_f10, , drop = FALSE]
  
  if (nrow(r_squared_filtered) == 0) {
    return(numeric(0))
  }
  
  # Get exposure names and convert to matrix
  exposure_names <- r_squared_filtered$Exposure
  r_mat <- as.matrix(r_squared_filtered[, setdiff(names(r_squared_filtered), "Exposure"), drop = FALSE])
  rownames(r_mat) <- exposure_names
  outcome_names <- colnames(r_mat)
  
  # For each exposure, compute rank of self-regression R²
  ranks <- numeric(length(exposure_names))
  names(ranks) <- exposure_names
  
  for (i in seq_along(exposure_names)) {
    exp_name <- exposure_names[i]
    
    # Get self-regression R²
    self_r2 <- NA
    if (exp_name %in% outcome_names) {
      col_idx <- which(outcome_names == exp_name)
      if (length(col_idx) > 0) {
        self_r2 <- r_mat[i, col_idx[1]]
      }
    }
    
    # Get all R² values for this exposure
    all_r2 <- as.numeric(r_mat[i, ])
    all_r2 <- all_r2[!is.na(all_r2)]
    
    if (is.na(self_r2) || length(all_r2) == 0) {
      ranks[i] <- NA
      next
    }
    
    # Rank self-regression R² among all R² values (higher R² = better rank)
    sorted_r2 <- sort(all_r2, decreasing = TRUE)
    rank_positions <- which(sorted_r2 == self_r2)
    
    if (length(rank_positions) > 0) {
      ranks[i] <- min(rank_positions)
    } else {
      ranks[i] <- NA
    }
  }
  
  # Remove NA ranks
  ranks <- ranks[!is.na(ranks)]
  return(ranks)
}

# Create TL-CFMR adjacency matrix where exposure-GSes have F-statistic > 10, Bonferroni-significant P-values, and R² rank in top 3
create_adjmatrix <- function(input_fstat, input_pval, input_r2, input_summary, fstat_col, suffix = "") {
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
  
  # Compute R² ranks and filter to top 3 for Bonferroni + R²T3 matrix
  input_r2_clean <- remove_second_and_later_duplicates(input_r2)
  r2_ranks <- compute_r2_ranks(input_r2_clean, input_summary, fstat_col)
  top3_exposures <- names(r2_ranks)[r2_ranks <= 3]
  
  # Save row names as Exposure (convert to column format for saving)
  pvalues_nominal <- pvalues_nominal %>% tibble::rownames_to_column("Exposure")
  pvalues_bonferroni <- pvalues_bonferroni %>% tibble::rownames_to_column("Exposure")
  
  # Filter Bonferroni matrix to top 3 R² rank exposures (rows only)
  if (length(top3_exposures) > 0) {
    # Filter rows (Exposure column) to top 3 R² rank exposures
    pvalues_bonferroni_r2t3 <- pvalues_bonferroni %>%
      dplyr::filter(Exposure %in% top3_exposures)
  } else {
    # If no exposures meet criteria, create empty data frame with same column structure
    pvalues_bonferroni_r2t3 <- pvalues_bonferroni[0, , drop = FALSE]
  }
  
  # Save the adjacency matrices with appropriate suffix
  if (suffix == "") {
    write.table(pvalues_nominal, file = paste0(output_path, "fstat10-nominal-totalfx-matrix.csv"), row.names = F, sep = ",") 
    write.table(pvalues_bonferroni, file = paste0(output_path, "fstat10-bonferroni-totalfx-matrix.csv"), row.names = F, sep = ",") 
    write.table(pvalues_bonferroni_r2t3, file = paste0(output_path, "fstat10-bonferroni-r2T3-totalfx-matrix.csv"), row.names = F, sep = ",") 
  } else {
    write.table(pvalues_nominal, file = paste0(output_path, "fstat10-nominal-totalfx-", suffix, "-matrix.csv"), row.names = F, sep = ",") 
    write.table(pvalues_bonferroni, file = paste0(output_path, "fstat10-bonferroni-totalfx-", suffix, "-matrix.csv"), row.names = F, sep = ",") 
    write.table(pvalues_bonferroni_r2t3, file = paste0(output_path, "fstat10-bonferroni-r2T3-totalfx-", suffix, "-matrix.csv"), row.names = F, sep = ",") 
  }
}

# Process unadjusted models
writeLines("Processing unadjusted models...")
adjmatrix_unadjusted <- create_adjmatrix(input_fstat_unadjusted, input_pval_unadjusted, input_r2_unadjusted, input_summary_unadjusted, "F_stat_unadjusted", suffix = "unadjusted")

# Process adjusted models
writeLines("Processing adjusted models...")
adjmatrix_adjusted <- create_adjmatrix(input_fstat_adjusted, input_pval_adjusted, input_r2_adjusted, input_summary_adjusted, "F_stat_adjusted", suffix = "agesex")

writeLines(c(" ", " ", "Finished script 07.1 - Formatted data for graphs.", " ", " "))
