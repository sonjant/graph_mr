#PBS -N j_02_1_split_into_folds
#PBS -l walltime=09:00:00
#PBS -l select=1:ncpus=10:mem=10gb
#PBS -j oe
#PBS -koed 

eval "$("${CONDA_HOME:-$HOME/miniforge3}/bin/conda" shell.bash hook)"
conda activate "${CONDA_ENV:-r443}"

[ -z "$GRAPH_MR_ROOT" ] && { echo "Error: GRAPH_MR_ROOT must be set (see README)."; exit 1; }
cd "$GRAPH_MR_ROOT"
Rscript 02_extract_cross-fitted_instruments/02.1_split_into_folds.R "$METHOD_TYPE" "$INPUT_TYPE"