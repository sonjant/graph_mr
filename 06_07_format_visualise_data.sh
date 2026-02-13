#PBS -N j_06_07_format_visualise_data
#PBS -l walltime=09:00:00
#PBS -l select=1:ncpus=8:mem=20gb
#PBS -j oe
#PBS -koed 

eval "$("${CONDA_HOME:-$HOME/miniforge3}/bin/conda" shell.bash hook)"
conda activate "${CONDA_ENV:-r443}"

[ -z "$GRAPH_MR_ROOT" ] && { echo "Error: GRAPH_MR_ROOT must be set (see README)."; exit 1; }
cd "$GRAPH_MR_ROOT"

Rscript 06.1_convert_hmdbid_to_chem_name.R "$METHOD_TYPE" "$INPUT_TYPE"
Rscript 06.2_collate_instrument_summary_table.R "$METHOD_TYPE" "$INPUT_TYPE"
Rscript 07.1_format_for_graph.R "$METHOD_TYPE" "$INPUT_TYPE"
Rscript 07.2_format_adjmat_to_cytoscape.R "$METHOD_TYPE" "$INPUT_TYPE" 
