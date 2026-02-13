#PBS -N j_02_5_load_clean_PLINK_results
#PBS -l walltime=03:00:00
#PBS -l select=1:ncpus=3:mem=4gb
#PBS -j oe
#PBS -koed 

eval "$("${CONDA_HOME:-$HOME/miniforge3}/bin/conda" shell.bash hook)"
conda activate "${CONDA_ENV:-r443}"

[ -z "$GRAPH_MR_ROOT" ] && { echo "Error: GRAPH_MR_ROOT must be set (see README)."; exit 1; }
cd "$GRAPH_MR_ROOT"
Rscript 02_extract_cross-fitted_instruments/02.5_load_clean_PLINK_results.R "$METHOD_TYPE" "$INPUT_TYPE"
