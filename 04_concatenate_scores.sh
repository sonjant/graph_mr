#PBS -N j_04_concatenate_scores
#PBS -l walltime=09:00:00
#PBS -l select=1:ncpus=8:mem=10gb
#PBS -j oe
#PBS -koed 

eval "$("${CONDA_HOME:-$HOME/miniforge3}/bin/conda" shell.bash hook)"
conda activate "${CONDA_ENV:-r443}"

[ -z "$GRAPH_MR_ROOT" ] && { echo "Error: GRAPH_MR_ROOT must be set (see README)."; exit 1; }
cd "$GRAPH_MR_ROOT"
Rscript 04_concatenate_scores.R "$METHOD_TYPE" "$INPUT_TYPE" 