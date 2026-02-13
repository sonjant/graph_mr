#!/bin/sh
#PBS -N j_03_7_extract_dosages
#PBS -l walltime=09:00:00
#PBS -l select=1:ncpus=22:mem=10gb
#PBS -j oe
#PBS -koed 

# Load conda environment
eval "$("${CONDA_HOME:-$HOME/miniforge3}/bin/conda" shell.bash hook)"
conda activate "${CONDA_ENV:-r443}"

[ -z "$GRAPH_MR_ROOT" ] && { echo "Error: GRAPH_MR_ROOT must be set (see README)."; exit 1; }
[ -z "$GRAPH_MR_VCF_DIR" ] && { echo "Error: GRAPH_MR_VCF_DIR must be set (see README)."; exit 1; }

RESULTS_BASE="${GRAPH_MR_RESULTS_DIR:-$GRAPH_MR_ROOT/results}"
EPHEMERAL_BASE="${GRAPH_MR_EPHEMERAL_DIR:-$GRAPH_MR_ROOT/ephemeral/}results"

snp_list="${RESULTS_BASE}/${METHOD_TYPE}_${INPUT_TYPE}/instrument_dosages/transfer_learned_instruments/all_unique_instruments_list.txt"
vcf_dir="${GRAPH_MR_VCF_DIR}"
output_dir="${EPHEMERAL_BASE}/${METHOD_TYPE}_${INPUT_TYPE}/instrument_dosages/transfer_learned_instruments"
mkdir -p "$output_dir"

# Function to extract dosages from one chromosome 
process_chr() {
  chr=$1
  echo "Processing chromosome ${chr}..."
  input_vcf="${vcf_dir}/merged.chr${chr}.1kg.vcf.gz"
  selected_variants="${output_dir}/selected_variants.chr${chr}.vcf.gz"
  dosage_file="${output_dir}/dosages.chr${chr}.txt"
  bcftools view -R "$snp_list" "$input_vcf" -Oz -o "$selected_variants"
  bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t[%DS\t]\n' "$selected_variants" > "$dosage_file"
  echo "Finished processing chromosome ${chr}."
}

export -f process_chr
export snp_list
export vcf_dir
export output_dir
export log_file

# Run the processing in parallel and log job details
parallel process_chr ::: {1..22}

# Combine all the dosage files into one file
cat ${output_dir}/dosages.chr{1..22}.txt > ${output_dir}/all_instrument_dosages.txt

echo "All chromosomes processed and combined into dosage-file.txt"

