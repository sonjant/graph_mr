# Purpose: Convert HMDBIDs and Uniprot IDs to chemical names

# Clean environment
rm(list=ls())

# Paths
if (Sys.getenv("GRAPH_MR_ROOT") == "") stop("Set GRAPH_MR_ROOT to the pipeline root directory (see README).")
setwd(Sys.getenv("GRAPH_MR_ROOT"))
source("config/config.R")

# Specify number of cores
ncores <- 5

# Load packages 
library("tidyverse")
library("vroom") 
library("parallel")
library("data.table")
library("readxl")

# Arguments: method_type, dataset
args <- commandArgs(trailingOnly = TRUE)
method_type <- args[1]
dataset     <- if (length(args) >= 2L) args[2] else Sys.getenv("GRAPH_MR_DATASET")
if (!nzchar(dataset)) stop("Provide dataset name or set GRAPH_MR_DATASET.")
base_dir <- paste0(method_type, "_", dataset)

# Exposure metadata (generic: exposure_id -> display_name). If missing, exposure IDs used as names.
meta_path <- paste0(data_dir, base_dir, "/preprocessed/exposure_metadata.tsv")
if (file.exists(meta_path)) {
  meta <- vroom(meta_path, show_col_types = FALSE)
  if ("exposure_id" %in% names(meta) && "display_name" %in% names(meta)) {
    combined_metadata <- meta %>% dplyr::rename(Exposure = exposure_id, display_name = display_name)
  } else if ("Exposure" %in% names(meta)) {
    display_col <- setdiff(names(meta), "Exposure")[1]
    if (is.na(display_col)) display_col <- "Exposure"
    combined_metadata <- meta %>% dplyr::rename(display_name = all_of(display_col))
  } else {
    combined_metadata <- data.frame(Exposure = character(), display_name = character())
  }
} else {
  combined_metadata <- data.frame(Exposure = character(), display_name = character())
}

# Graph MR result matrices and output
matrices <- list.files(paste0(results_dir, base_dir, "/matrices"), full.names = TRUE)
output <- paste0(results_dir, base_dir, "/matrices_chem_names/")
dir.create(output, recursive = TRUE, showWarnings = FALSE)



# Start script ----
writeLines(c(" ", " ", "Started script 06.1 - Convert HMDBIDs to chemical names.", " ", " "))

# Translate exposure IDs to display names using exposure_metadata (generic)
convert_names <- function(matrices) {
  if (nrow(combined_metadata) == 0L) {
    for (m in matrices) {
      mat <- vroom(m, show_col_types = FALSE)
      if (any(grepl("_scores", mat$Exposure))) mat$Exposure <- sub("_scores", "", mat$Exposure)
      vroom_write(mat, paste0(output, basename(m)), delim = ",", append = FALSE)
    }
    return(invisible(NULL))
  }
  mclapply(matrices, mc.cores = ncores, FUN = function(matrix) {
    matrix_name <- basename(matrix)
    matrix <- vroom(matrix, show_col_types = FALSE)
    if (any(grepl("_scores", matrix$Exposure))) matrix$Exposure <- sub("_scores", "", matrix$Exposure)
    matrix <- matrix %>%
      dplyr::left_join(combined_metadata, by = "Exposure") %>%
      dplyr::mutate(Exposure = dplyr::coalesce(display_name, Exposure)) %>%
      dplyr::select(-dplyr::any_of("display_name"))
    translated_colnames <- sapply(colnames(matrix), function(col) {
      idx <- match(col, combined_metadata$Exposure)
      if (!is.na(idx) && "display_name" %in% names(combined_metadata))
        combined_metadata$display_name[idx] else col
    })
    colnames(matrix) <- translated_colnames
    vroom_write(matrix, paste0(output, matrix_name), delim = ",", append = FALSE)
    cat("Created", matrix_name, "\n")
  })
}

converted_name_results <- convert_names(matrices)

writeLines(c(" ", " ", "Finished script 06.1 - We have converted HMDBIDs to chemical names.", " ", " "))
