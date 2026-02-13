# Purpose: Generic preprocessing for Graph MR.
# Takes one exposure matrix + covariates + genotype PCs and produces the standard
# layout (preprocessed/exposures.tsv, gwas_input/...) for any dataset.
# No study-specific logic (no airwave, metabolomics, etc.).

# Clean environment
rm(list = ls())

# Paths
if (Sys.getenv("GRAPH_MR_ROOT") == "") stop("Set GRAPH_MR_ROOT to the pipeline root directory (see README).")
setwd(Sys.getenv("GRAPH_MR_ROOT"))
source("config/config.R")
if (!nzchar(genotype_dir)) stop("Set GRAPH_MR_GENOTYPE_DIR (see README).")

# Load packages
library("tidyverse")
library("vroom")

# Arguments: method_type, dataset (dataset = label for this run, e.g. proteomics, mydata)
args <- commandArgs(trailingOnly = TRUE)
method_type <- args[1]
dataset     <- args[2]
if (is.na(dataset) || length(dataset) < 1L) dataset <- Sys.getenv("GRAPH_MR_DATASET")
if (!nzchar(dataset)) stop("Provide dataset name as second argument or set GRAPH_MR_DATASET.")

base_dir <- paste0(method_type, "_", dataset)
method_types_str <- Sys.getenv("METHOD_TYPES", unset = "")
if (nzchar(method_types_str)) {
  method_types_list <- trimws(strsplit(method_types_str, ",")[[1]])
} else {
  method_types_list <- method_type
}

# User input paths (required)
exposure_matrix_path <- Sys.getenv("GRAPH_MR_EXPOSURE_MATRIX")
covariate_path       <- Sys.getenv("GRAPH_MR_COVARIATE_FILE")
if (!nzchar(exposure_matrix_path)) stop("Set GRAPH_MR_EXPOSURE_MATRIX to path to exposure matrix (sample_id + exposure columns).")
if (!nzchar(covariate_path))       stop("Set GRAPH_MR_COVARIATE_FILE to path to covariate file (sample_id, age, sex).")

# Optional: column names in user files (defaults)
sample_id_col  <- Sys.getenv("GRAPH_MR_SAMPLE_ID_COL",  unset = "barcode")
age_col        <- Sys.getenv("GRAPH_MR_AGE_COL",        unset = "AGE")
sex_col        <- Sys.getenv("GRAPH_MR_SEX_COL",        unset = "SEX")
exposure_meta_path <- Sys.getenv("GRAPH_MR_EXPOSURE_METADATA", unset = "")

# Output directories under data_dir/dataset/
output_preprocessed <- paste0(data_dir, dataset, "/preprocessed/")
output_gwas_input   <- paste0(data_dir, dataset, "/gwas_input/")
for (mt in method_types_list) {
  dir.create(paste0(data_dir, mt, "_", dataset, "/preprocessed/"), recursive = TRUE, showWarnings = FALSE)
  dir.create(paste0(data_dir, mt, "_", dataset, "/gwas_input/"), recursive = TRUE, showWarnings = FALSE)
}

# Read inputs
exposures <- vroom(exposure_matrix_path, show_col_types = FALSE)
covariates <- vroom(covariate_path, show_col_types = FALSE)
gwas_pcs <- vroom(paste0(genotype_dir, "merged.pruned.eigenvec"), show_col_types = FALSE) %>%
  dplyr::rename(barcode = "IID", FID = "#FID")

# Normalize sample ID column name in user data (pipeline uses "barcode" internally)
if (!sample_id_col %in% colnames(exposures)) {
  if ("barcode" %in% colnames(exposures)) sample_id_col <- "barcode"
  else if ("IID" %in% colnames(exposures)) sample_id_col <- "IID"
  else stop("Exposure matrix must have a sample ID column (set GRAPH_MR_SAMPLE_ID_COL if different from 'barcode').")
}
if (sample_id_col != "barcode") exposures <- dplyr::rename(exposures, barcode = all_of(sample_id_col))

if (!age_col %in% colnames(covariates)) age_col <- grep("age|AGE", colnames(covariates), ignore.case = TRUE, value = TRUE)[1]
if (!sex_col %in% colnames(covariates)) sex_col <- grep("sex|gender|SEX", colnames(covariates), ignore.case = TRUE, value = TRUE)[1]
if (is.na(age_col) || is.na(sex_col)) stop("Covariate file must contain age and sex columns.")
covariates <- covariates %>%
  dplyr::select(any_of(c("barcode", "IID", sample_id_col)), all_of(c(age_col, sex_col))) %>%
  dplyr::rename(AGE = all_of(age_col), SEX = all_of(sex_col))
if (!"barcode" %in% colnames(covariates) && "IID" %in% colnames(covariates)) covariates <- dplyr::rename(covariates, barcode = IID)
if (!"barcode" %in% colnames(covariates) && length(sample_id_col) > 0 && sample_id_col %in% colnames(covariates))
  covariates <- dplyr::rename(covariates, barcode = all_of(sample_id_col))

# Merge: exposures + covariates + PCs (restrict to samples in genotype)
merged <- exposures %>%
  dplyr::inner_join(covariates, by = "barcode") %>%
  dplyr::inner_join(gwas_pcs %>% dplyr::select(barcode, FID, PC1, PC2, PC3, PC4, PC5), by = "barcode") %>%
  dplyr::filter(!is.na(barcode), !duplicated(barcode))

# GWAS covariate file (PLINK-style: FID, IID, PC1..PC5, AGE, SEX)
gwas_covariates <- merged %>%
  dplyr::mutate(SEX = case_when(SEX %in% c("M", "male", 1) ~ 1L, SEX %in% c("F", "female", 2) ~ 2L, TRUE ~ NA_integer_)) %>%
  dplyr::select(barcode, FID, PC1, PC2, PC3, PC4, PC5, AGE, SEX) %>%
  dplyr::filter(complete.cases(.))

# Covariates for MR step (05.1): barcode, AGE, SEX
covariates_age_sex <- merged %>% dplyr::select(barcode, AGE, SEX)
if ("SEX" %in% names(covariates_age_sex) && is.character(covariates_age_sex$SEX))
  covariates_age_sex <- dplyr::mutate(covariates_age_sex, SEX = case_when(SEX %in% c("M", "male") ~ 1L, SEX %in% c("F", "female") ~ 2L, TRUE ~ NA_integer_))

# Preprocessed exposures: barcode + standardized exposure columns (exclude non-exposure)
exposure_cols <- setdiff(colnames(merged), c("barcode", "AGE", "SEX", "FID", "PC1", "PC2", "PC3", "PC4", "PC5"))
exposure_only <- merged %>% dplyr::select(barcode, all_of(exposure_cols))
# Standardize continuous columns (leave binary 0/1 as-is)
exposure_standardized <- exposure_only %>%
  dplyr::mutate(across(
    -barcode,
    ~ {
      u <- unique(na.omit(.))
      if (length(u) <= 2L && all(u %in% c(0, 1))) . else as.vector(scale(.))
    }
  ))

# Write outputs to all method_type folders
for (mt in method_types_list) {
  out_pre <- paste0(data_dir, mt, "_", dataset, "/preprocessed/")
  out_gwas <- paste0(data_dir, mt, "_", dataset, "/gwas_input/")
  dir.create(out_pre, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_gwas, recursive = TRUE, showWarnings = FALSE)
  vroom_write(exposure_standardized, paste0(out_pre, "exposures.tsv"), na = "")
  vroom_write(gwas_covariates, paste0(out_gwas, "gwas_covariates_pc1-5_age_sex.txt"))
  vroom_write(covariates_age_sex, paste0(out_gwas, "covariates_age_sex.txt"))
}

# Exposure metadata: exposure_id, display_name (optional: pathway, class for 07.2 colours)
if (nzchar(exposure_meta_path) && file.exists(exposure_meta_path)) {
  meta <- vroom(exposure_meta_path, show_col_types = FALSE)
  if (!"exposure_id" %in% names(meta)) meta <- dplyr::rename(meta, exposure_id = 1)
  if (!"display_name" %in% names(meta) && ncol(meta) >= 2) meta <- dplyr::rename(meta, display_name = 2)
  if (!"display_name" %in% names(meta)) meta$display_name <- meta$exposure_id
} else {
  meta <- data.frame(exposure_id = exposure_cols, display_name = exposure_cols, stringsAsFactors = FALSE)
}
for (mt in method_types_list)
  vroom_write(meta, paste0(data_dir, mt, "_", dataset, "/preprocessed/exposure_metadata.tsv"))

message("Generic preprocessing done for dataset: ", dataset)
