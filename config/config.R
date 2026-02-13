# Graph MR pipeline – central path configuration
# Source this from R scripts after setting GRAPH_MR_ROOT (see README and config/config.example.sh).
# Do not hardcode user-specific paths in scripts; use the variables defined here.

root <- Sys.getenv("GRAPH_MR_ROOT")
if (root == "" || is.na(root)) {
  stop("GRAPH_MR_ROOT environment variable must be set to the pipeline root directory (e.g. the graph_mr_github folder). See README.")
}
root <- normalizePath(root, winslash = "/", mustWork = TRUE)
if (!endsWith(root, "/")) root <- paste0(root, "/")

# Project root (working directory for relative paths like data/<dataset>/)
home <- root
setwd(home)

# Genotype/PLINK data directory (e.g. containing merged.pruned.eigenvec and per-chromosome bed/bim/fam).
# Must be set by the user for steps that need it; no default.
genotype_dir <- Sys.getenv("GRAPH_MR_GENOTYPE_DIR", unset = NA)
if (is.na(genotype_dir)) genotype_dir <- ""
if (nzchar(genotype_dir)) {
  genotype_dir <- normalizePath(genotype_dir, winslash = "/", mustWork = FALSE)
  if (!endsWith(genotype_dir, "/")) genotype_dir <- paste0(genotype_dir, "/")
}

# Ephemeral (intermediate) results directory. Default: <repo>/ephemeral
ephemeral_dir <- Sys.getenv("GRAPH_MR_EPHEMERAL_DIR")
if (ephemeral_dir == "" || is.na(ephemeral_dir)) {
  ephemeral_dir <- paste0(home, "ephemeral/")
} else {
  ephemeral_dir <- normalizePath(ephemeral_dir, winslash = "/", mustWork = FALSE)
  if (!endsWith(ephemeral_dir, "/")) ephemeral_dir <- paste0(ephemeral_dir, "/")
}

# Final results directory. Default: <repo>/results
results_dir <- Sys.getenv("GRAPH_MR_RESULTS_DIR")
if (results_dir == "" || is.na(results_dir)) {
  results_dir <- paste0(home, "results/")
} else {
  results_dir <- normalizePath(results_dir, winslash = "/", mustWork = FALSE)
  if (!endsWith(results_dir, "/")) results_dir <- paste0(results_dir, "/")
}

# Data directory: under this, one folder per dataset (e.g. data/proteomics/, data/mydata/)
# Pipeline reads from data_dir/<dataset>/preprocessed/exposures.tsv and data_dir/<dataset>/gwas_input/
data_dir <- Sys.getenv("GRAPH_MR_DATA_DIR")
if (data_dir == "" || is.na(data_dir)) {
  data_dir <- paste0(home, "data/")
} else {
  data_dir <- normalizePath(data_dir, winslash = "/", mustWork = FALSE)
  if (!endsWith(data_dir, "/")) data_dir <- paste0(data_dir, "/")
}

# Optional: annotation dir for omic-specific files (e.g. proteomics Olink coords under annotation/OlinkProteomics/).
annotation_dir <- Sys.getenv("GRAPH_MR_ANNOTATION_DIR")
if (annotation_dir == "" || is.na(annotation_dir)) annotation_dir <- paste0(home, "data/annotation/")
if (!endsWith(annotation_dir, "/")) annotation_dir <- paste0(annotation_dir, "/")

# Dataset name: label for this run (e.g. proteomics, mydata). No fixed list; user chooses.
# Set via GRAPH_MR_DATASET or pass as second argument to scripts (method_type, dataset).
# Paths use data_dir and dataset, e.g. paste0(data_dir, dataset, "/preprocessed/exposures.tsv")
