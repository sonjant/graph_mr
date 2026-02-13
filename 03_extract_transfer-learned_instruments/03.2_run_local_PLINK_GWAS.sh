# Purpose: Run GWASes for n number of exposures (X) and k number of folds

#!/bin/bash
#PBS -N j_03_2_run_PLINK_GWAS
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=128:mem=450gb
#PBS -j oe
#PBS -koed 

eval "$("${CONDA_HOME:-$HOME/miniforge3}/bin/conda" shell.bash hook)"
conda activate "${CONDA_ENV:-r443}"

[ -z "$GRAPH_MR_ROOT" ] && { echo "Error: GRAPH_MR_ROOT must be set (see README)."; exit 1; }
DATA_DIR="${GRAPH_MR_DATA_DIR:-$GRAPH_MR_ROOT/data}"
EPHEMERAL_BASE="${GRAPH_MR_EPHEMERAL_DIR:-$GRAPH_MR_ROOT/ephemeral/}results"

cd "${EPHEMERAL_BASE}/${METHOD_TYPE}_${INPUT_TYPE}/"
cp -r imputed-beagle-1kg_v3_cleaned/ "${TMPDIR}"
input="${TMPDIR}/imputed-beagle-1kg_v3_cleaned/"

cd "${DATA_DIR}/${METHOD_TYPE}_${INPUT_TYPE}/"
cp -r gwas_input "${TMPDIR}/gwas_input/"
input_participants="${TMPDIR}/gwas_input/"

cd "${EPHEMERAL_BASE}/${METHOD_TYPE}_${INPUT_TYPE}"
cp -r gwas_results "${TMPDIR}/gwas_results/"
output="${TMPDIR}/gwas_results/"

cd "${TMPDIR}"

# Start script ---- ---- ---- ---- ---- ----

# Calculate N if not provided (fallback to counting from exposure_list.txt)
if [ -z "$N" ]; then
    N=$(($(wc -l < ${input_participants}exposure_list.txt) + 1))
fi

# IN PARALLEL - run a GWAS for all 22 chromosomes
for i in $(seq 1 $N)
do
    exposure=$(sed -n ${i}p < ${input_participants}exposure_list.txt) # Create a variable for the exposure name by indexing each name in the .txt file.
    echo "starting exposure $exposure"
    
    mkdir -p "${output}${exposure}/complete_cases"
    
    parallel plink2   --bfile ${input}cleaned_chr{1} \
                      --ci 0.95 \
                      --covar ${input_participants}gwas_covariates_pc1-5_age_sex.txt \
                      --covar-variance-standardize \
                      --glm hide-covar cols=+a1freq,+beta omit-ref \
                      --maf 0.005 \
                      --pfilter 5e-8 \
                      --pheno ${input_participants}${exposure}/complete_cases_gwas_participants.txt \
                      --out ${output}${exposure}/complete_cases/chr{1} ::: {1..22}
                      
    echo "Finished exposure $exposure. Will now copy results back to ephemeral"
    cp -r "${output}${exposure}" "${EPHEMERAL_BASE}/${METHOD_TYPE}_${INPUT_TYPE}/gwas_results/"
done

echo -e "\nFinished running GWASs for exposures 1 to $N \n"