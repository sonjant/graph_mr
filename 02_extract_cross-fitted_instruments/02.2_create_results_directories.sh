# Purpose: Create result directories for the GWAS results.

#!/bin/bash
#PBS -N j_02_2_create_res_directories
#PBS -l walltime=09:00:00
#PBS -l select=1:ncpus=2:mem=2gb
#PBS -j oe
#PBS -koed 

[ -z "$GRAPH_MR_ROOT" ] && { echo "Error: GRAPH_MR_ROOT must be set (see README)."; exit 1; }
cd "$GRAPH_MR_ROOT"

# Set the number of folds, k, you wish to use for generating cross-fitted genetic risk scores
k=10

# Set the file path of your split exposure files that were the output in script 02.1
DATA_DIR="${GRAPH_MR_DATA_DIR:-$GRAPH_MR_ROOT/data}"
input_participants="${DATA_DIR}/${METHOD_TYPE}_${INPUT_TYPE}/gwas_input/"

echo "Using input directory: $input_participants"

# Check if METHOD_TYPES environment variable is set (for creating directories for multiple methods)
if [ -n "$METHOD_TYPES" ]; then
  # Parse comma-separated method types
  IFS=',' read -ra METHOD_ARRAY <<< "$METHOD_TYPES"
  METHOD_ARRAY=("${METHOD_ARRAY[@]// /}")  # Remove whitespace
else
  # Use single METHOD_TYPE if not specified
  METHOD_ARRAY=("$METHOD_TYPE")
fi

# Create output directories for all method types
for MT in "${METHOD_ARRAY[@]}"; do
  output="${GRAPH_MR_EPHEMERAL_DIR:-$GRAPH_MR_ROOT/ephemeral/}results/${MT}_${INPUT_TYPE}/gwas_results/"
  echo "Creating directories for method: $MT"



# Start script ---- ---- ---- ---- ---- ----
# Search the directory of exposures and add their names to a .txt file
for dir in ${input_participants}/*/; do 
  basename "$dir"
done > ${input_participants}exposure_list.txt

  # For each exposure, create result directories in the ephemeral directory for n number of folds and a folder for complete cases.
  while read name; do 
    mkdir -p "$output$name"
    mkdir -p "$output$name/complete_cases"
    for i in $(seq 1 $k); do 
      mkdir -p "$output$name/fold$i" 
    done
  done < ${input_participants}exposure_list.txt
done

# Remove the trailing newline at the metabolite text file
truncate -s -1 ${input_participants}exposure_list.txt

echo "PLINK results folders created for all method types"