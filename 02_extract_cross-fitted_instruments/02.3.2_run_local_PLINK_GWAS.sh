# Purpose: Run GWASes for n number of exposures (X) and k number of folds

#!/bin/bash
#PBS -N j_0232_run_PLINK_GWAS
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

# Calculate the range for this script (second third: N/3+1 to 2*N/3)
N_THIRD=$((N / 3))
START=$((N_THIRD + 1))
END=$((2 * N_THIRD))

# IN PARALLEL - run a GWAS for all 22 chromosomes on folds 1 to n
# Split gwas_results into 3 chunks
# Do processing then copy the contents of gwas_results{i} back to RDS

for i in $(seq $START $END)
do
    exposure=$(sed -n ${i}p < ${input_participants}exposure_list.txt | tr -d '\r\n' | xargs) # Create a variable for the exposure name by indexing each name in the .txt file, trim whitespace and newlines
    
    # Skip if exposure is empty
    if [ -z "$exposure" ]; then
        echo "Warning: Exposure at line $i is empty, skipping..."
        continue
    fi
    
    echo "starting exposure $exposure"

    for n in {1..10}
    do
        # Create output directory if it doesn't exist
        mkdir -p "${output}${exposure}/fold${n}"
        
        parallel plink2 \
            --bfile "${input}cleaned_chr{1}" \
            --ci 0.95 \
            --covar "${input_participants}gwas_covariates_pc1-5_age_sex.txt" \
            --covar-variance-standardize \
            --glm hide-covar cols=+a1freq,+beta omit-ref \
            --maf 0.005 \
            --pfilter 5e-8 \
            --pheno "${input_participants}${exposure}/fold${n}_gwas_participants.txt" \
            --out "${output}${exposure}/fold${n}/chr{1}" ::: {1..22}
    done

    echo "finished exposure $exposure. Will now copy results back to ephemeral"
    cp -r "${output}${exposure}" "${EPHEMERAL_BASE}/${METHOD_TYPE}_${INPUT_TYPE}/gwas_results/"
done

echo -e "\nFinished running GWASs for exposures $START to $END \n"


