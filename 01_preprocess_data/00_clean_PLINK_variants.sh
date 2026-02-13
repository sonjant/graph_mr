#!/bin/bash
#PBS -N j_00_clean_variants
#PBS -l walltime=01:30:00
#PBS -l select=1:ncpus=4:mem=8gb
#PBS -j oe
#PBS -koed 

# Conda: set CONDA_HOME and CONDA_ENV if needed (e.g. in config.sh)
eval "$("${CONDA_HOME:-$HOME/miniforge3}/bin/conda" shell.bash hook)"
conda activate "${CONDA_ENV:-r443}"

# Check required environment variables
if [ -z "$GRAPH_MR_ROOT" ]; then
    echo "Error: GRAPH_MR_ROOT must be set to the pipeline root directory (see README)."
    exit 1
fi
if [ -z "$GRAPH_MR_GENOTYPE_DIR" ]; then
    echo "Error: GRAPH_MR_GENOTYPE_DIR must be set to the genotype/PLINK directory (see README)."
    exit 1
fi
if [ -z "$METHOD_TYPE" ] || [ -z "$INPUT_TYPE" ]; then
    echo "Error: METHOD_TYPE and INPUT_TYPE must be set"
    exit 1
fi

cd "$GRAPH_MR_ROOT"

# Assign input and output paths
input="${GRAPH_MR_GENOTYPE_DIR}"
output="${GRAPH_MR_EPHEMERAL_DIR:-$GRAPH_MR_ROOT/ephemeral/}results/${METHOD_TYPE}_${INPUT_TYPE}/imputed-beagle-1kg_v3_cleaned/"

echo "Starting PLINK variant cleaning..."
echo "METHOD_TYPE=$METHOD_TYPE"
echo "INPUT_TYPE=$INPUT_TYPE"
echo "Input directory: $input"
echo "Output directory: $output"

# Create output directories
mkdir -p "$output"
echo "Created output directory: $output"

# Use plink2 to fix variant names and filter variants by missingness
for chr in {1..22}; do
  echo "Processing chromosome $chr..."
  plink2 \
    --bfile "${input%/}/merged.chr${chr}.1kg" \
    --set-missing-var-ids '@:\#:\$r:\$a' \
    --new-id-max-allele-len 100 truncate \
    --snps-only \
    --geno 0.1 \
    --make-bed \
    --out "${output}cleaned_chr${chr}"
  if [ $? -ne 0 ]; then
    echo "Error: plink2 failed for chromosome $chr"
    exit 1
  fi
  echo "Completed chromosome $chr"
done

echo "All chromosomes processed successfully!"
