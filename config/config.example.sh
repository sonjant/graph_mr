# Example environment variables for the Graph MR pipeline.
# Copy to config.sh, fill in your paths, and source before running:
#   cp config/config.example.sh config/config.sh
#   # edit config/config.sh with your paths
#   source config/config.sh
# config.sh is gitignored so your paths are not committed.

# Required: path to the pipeline repo (graph_mr_github directory)
export GRAPH_MR_ROOT="/path/to/your/graph_mr_github"

# Required for steps that use genotype/PLINK data: directory containing
# merged.pruned.eigenvec and per-chromosome plink files (e.g. merged.chr1.1kg.bed/bim/fam)
export GRAPH_MR_GENOTYPE_DIR="/path/to/genotype/merged/bed"

# Required for 02.8 (extract dosages): directory containing merged.chr*.1kg.vcf.gz files
export GRAPH_MR_VCF_DIR="/path/to/genotype/merged/vcf"

# Optional: directory for large intermediate outputs (GWAS results, etc.).
# Default: ${GRAPH_MR_ROOT}/ephemeral
# export GRAPH_MR_EPHEMERAL_DIR="/path/to/ephemeral"

# Optional: directory for final pipeline results.
# Default: ${GRAPH_MR_ROOT}/results
# export GRAPH_MR_RESULTS_DIR="/path/to/results"

# Dataset name: label for your post-genotype omic (e.g. proteomics, metabolomics, mydata).
# Used as folder name under data/ and in results. No fixed list; you choose the name.
# export GRAPH_MR_DATASET="proteomics"

# Optional: data directory (default: ${GRAPH_MR_ROOT}/data)
# export GRAPH_MR_DATA_DIR="/path/to/data"

# --- For generic preprocessing only (01_preprocess_exposures_generic) ---
# Path to your exposure matrix (sample ID column + one column per exposure) and covariate file (sample ID, age, sex).
# export GRAPH_MR_EXPOSURE_MATRIX="/path/to/exposures.tsv"
# export GRAPH_MR_COVARIATE_FILE="/path/to/covariates_age_sex.txt"
# Optional: exposure metadata (exposure_id, display_name [, pathway/class])
# export GRAPH_MR_EXPOSURE_METADATA="/path/to/exposure_metadata.tsv"
# Optional column names if different from barcode, AGE, SEX:
# export GRAPH_MR_SAMPLE_ID_COL="barcode"
# export GRAPH_MR_AGE_COL="AGE"
# export GRAPH_MR_SEX_COL="SEX"

# Optional: conda environment name and path to conda (if not on PATH)
# export CONDA_ENV="r443"
# export CONDA_HOME="$HOME/miniforge3"
