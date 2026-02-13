#!/usr/bin/env bash
# Run generic preprocessing: exposure matrix + covariates -> standard layout under data/<dataset>/
# Set GRAPH_MR_EXPOSURE_MATRIX, GRAPH_MR_COVARIATE_FILE, and dataset (see docs/INPUT_SPEC.md).
# Usage: METHOD_TYPE=cftlmr INPUT_TYPE=mydata [GRAPH_MR_EXPOSURE_MATRIX=... GRAPH_MR_COVARIATE_FILE=...] ./01_preprocess_data/01_preprocess_exposures_generic.sh

set -e
[ -n "$GRAPH_MR_ROOT" ] || { echo "Error: GRAPH_MR_ROOT must be set (see README)."; exit 1; }
cd "$GRAPH_MR_ROOT"

DATASET="${INPUT_TYPE:-$GRAPH_MR_DATASET}"
[ -n "$DATASET" ] || { echo "Error: Set INPUT_TYPE or GRAPH_MR_DATASET to your dataset name (e.g. proteomics, mydata)."; exit 1; }
METHOD_TYPE="${METHOD_TYPE:-cftlmr}"
export METHOD_TYPES="${METHOD_TYPES:-cftlmr,cfmr}"

Rscript 01_preprocess_data/01_preprocess_exposures_generic.R "$METHOD_TYPE" "$DATASET"
echo "Generic preprocessing done. Run pipeline from 00_clean_PLINK_variants (then 02.1_split_into_folds) with INPUT_TYPE=$DATASET."
