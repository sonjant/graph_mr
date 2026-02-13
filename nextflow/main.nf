/*
 * Graph MR Nextflow pipeline (cftlmr + cfmr). Plug-and-play: genotype + one dataset.
 * Usage (from pipeline root): nextflow run nextflow/main.nf -c nextflow/nextflow.config --root $PWD --genotype_dir /path/to/plink --vcf_dir /path/to/vcf --dataset mydata
 */

nextflow.enable.dsl = 2

def data_dir = params.data_dir ?: "${params.root}/data"
def ephemeral_dir = params.ephemeral_dir ?: "${params.root}/ephemeral"
def results_dir = params.results_dir ?: "${params.root}/results"

def env_export = """
  export GRAPH_MR_ROOT="${params.root}"
  export GRAPH_MR_GENOTYPE_DIR="${params.genotype_dir}"
  export GRAPH_MR_VCF_DIR="${params.vcf_dir}"
  export GRAPH_MR_DATA_DIR="${data_dir}"
  export GRAPH_MR_EPHEMERAL_DIR="${ephemeral_dir}"
  export GRAPH_MR_RESULTS_DIR="${results_dir}"
  export METHOD_TYPE="${params.primary_method}"
  export INPUT_TYPE="${params.dataset}"
  export METHOD_TYPES="${params.method_types}"
"""
def env_graph_mr = """
  export GRAPH_MR_ROOT="${params.root}"
  export GRAPH_MR_GENOTYPE_DIR="${params.genotype_dir}"
  export GRAPH_MR_VCF_DIR="${params.vcf_dir}"
  export GRAPH_MR_DATA_DIR="${data_dir}"
  export GRAPH_MR_EPHEMERAL_DIR="${ephemeral_dir}"
  export GRAPH_MR_RESULTS_DIR="${results_dir}"
  export INPUT_TYPE="${params.dataset}"
  export METHOD_TYPES="${params.method_types}"
"""

// ---------------------------------------------------------------------------
// Common pipeline: 00 -> 02.1 -> 02.2 (no study-specific preprocess; use 01_preprocess_exposures_generic.sh first if needed)
// ---------------------------------------------------------------------------

process clean_plink_variants {
  tag "common"
  output: path('done.trigger'), optional: true
  script:
  """
  ${env_export}
  cd ${params.root} && bash 01_preprocess_data/00_clean_PLINK_variants.sh
  touch done.trigger
  """
}

process split_into_folds {
  tag "common"
  input: path('trigger')
  output: path('done.trigger'), optional: true
  script:
  """
  ${env_export}
  cd ${params.root} && bash 02_extract_cross-fitted_instruments/02.1_split_into_folds.sh
  touch done.trigger
  """
}

process create_results_directories {
  tag "common"
  input: path('trigger')
  output: path('done.trigger'), optional: true
  script:
  """
  ${env_export}
  cd ${params.root} && bash 02_extract_cross-fitted_instruments/02.2_create_results_directories.sh
  touch done.trigger
  """
}

process get_exposure_count {
  tag "common"
  input: path('trigger')
  output: stdout
  script:
  """
  EXPOSURE_LIST="${data_dir}/${params.primary_method}_${params.dataset}/gwas_input/exposure_list.txt"
  if [[ ! -f "\$EXPOSURE_LIST" ]]; then
    echo "Waiting for exposure_list.txt..." >&2
    sleep 30
  fi
  N=\$((\$(wc -l < "\$EXPOSURE_LIST") + 1))
  echo \$N
  """
}

process run_local_PLINK_GWAS_1 {
  tag "common"
  input: val(n)
  output: path('done.trigger'), optional: true
  script:
  """
  ${env_export}
  export N="$n"
  cd ${params.root} && bash 02_extract_cross-fitted_instruments/02.3.1_run_local_PLINK_GWAS.sh
  touch done.trigger
  """
}

process run_local_PLINK_GWAS_2 {
  tag "common"
  input: val(n)
  output: path('done.trigger'), optional: true
  script:
  """
  ${env_export}
  export N="$n"
  cd ${params.root} && bash 02_extract_cross-fitted_instruments/02.3.2_run_local_PLINK_GWAS.sh
  touch done.trigger
  """
}

process run_local_PLINK_GWAS_3 {
  tag "common"
  input: val(n)
  output: path('done.trigger'), optional: true
  script:
  """
  export METHOD_TYPE="${params.primary_method}"
  export INPUT_TYPE="${params.dataset}"
  export METHOD_TYPES="${params.method_types}"
  export N="$n"
  cd ${params.root} && bash 02_extract_cross-fitted_instruments/02.3.3_run_local_PLINK_GWAS.sh
  touch done.trigger
  """
}

process run_PLINK_clumping {
  tag "common"
  input: path('t1'), path('t2'), path('t3')
  output: path('done.trigger'), optional: true
  script:
  """
  ${env_export}
  cd ${params.root} && bash 02_extract_cross-fitted_instruments/02.4_run_PLINK_clumping.sh
  touch done.trigger
  """
}

process load_clean_PLINK_results {
  tag "common"
  input: path('trigger'), val(task_id)
  output: path('done.trigger'), optional: true
  script:
  """
  export METHOD_TYPE="${params.primary_method}"
  export INPUT_TYPE="${params.dataset}"
  export METHOD_TYPES="${params.method_types}"
  export PBS_ARRAY_INDEX="$task_id"
  cd ${params.root} && bash 02_extract_cross-fitted_instruments/02.5_load_clean_PLINK_results.sh
  touch done.trigger
  """
}

process create_instruments_list {
  tag "common"
  input: path('trigger')
  output: path('done.trigger'), optional: true
  script:
  """
  ${env_export}
  cd ${params.root} && bash 02_extract_cross-fitted_instruments/02.6_create_instruments_list.sh
  touch done.trigger
  """
}

process create_unique_instruments_list {
  tag "common"
  input: path('trigger')
  output: path('done.trigger'), optional: true
  script:
  """
  ${env_export}
  cd ${params.root} && bash 02_extract_cross-fitted_instruments/02.7_create_unique_instruments_list.sh
  touch done.trigger
  """
}

process extract_dosages {
  tag "common"
  input: path('trigger')
  output: path('done.trigger'), optional: true
  script:
  """
  ${env_graph_mr}
  export METHOD_TYPE="${params.primary_method}"
  cd ${params.root} && bash 02_extract_cross-fitted_instruments/02.8_extract_dosages.sh
  touch done.trigger
  """
}

process create_cross_fitted_scores {
  tag "common"
  input: path('trigger'), val(batch_id)
  output: path('done.trigger'), optional: true
  script:
  """
  ${env_export}
  export PBS_ARRAY_INDEX="$batch_id"
  cd ${params.root} && bash 02_extract_cross-fitted_instruments/02.9_create_cross-fitted_scores.sh
  touch done.trigger
  """
}

// ---------------------------------------------------------------------------
// Transfer learning (cftlmr only): JOB13–JOB20
// ---------------------------------------------------------------------------

process create_missing_complete_lists {
  tag "cftlmr"
  output: path('done.trigger'), optional: true
  script:
  """
  ${env_export}
  cd ${params.root} && bash 03_extract_transfer-learned_instruments/03.1_create_missing_and_complete_cases_lists.sh
  touch done.trigger
  """
}

process run_local_PLINK_GWAS_tl {
  tag "cftlmr"
  input: path('trigger'), val(n)
  output: path('done.trigger'), optional: true
  script:
  """
  ${env_graph_mr}
  export METHOD_TYPE="cftlmr"
  export N="$n"
  cd ${params.root} && bash 03_extract_transfer-learned_instruments/03.2_run_local_PLINK_GWAS.sh
  touch done.trigger
  """
}

process run_PLINK_clumping_tl {
  tag "cftlmr"
  input: path('trigger')
  output: path('done.trigger'), optional: true
  script:
  """
  export METHOD_TYPE="cftlmr"
  export INPUT_TYPE="${params.dataset}"
  cd ${params.root} && bash 03_extract_transfer-learned_instruments/03.3_run_PLINK_clumping.sh
  touch done.trigger
  """
}

process load_clean_PLINK_results_tl {
  tag "cftlmr"
  input: path('trigger')
  output: path('done.trigger'), optional: true
  script:
  """
  ${env_graph_mr}
  export METHOD_TYPE="cftlmr"
  cd ${params.root} && bash 03_extract_transfer-learned_instruments/03.4_load_clean_PLINK_results.sh
  touch done.trigger
  """
}

process create_instruments_list_tl {
  tag "cftlmr"
  input: path('trigger')
  output: path('done.trigger'), optional: true
  script:
  """
  export METHOD_TYPE="cftlmr"
  export INPUT_TYPE="${params.dataset}"
  cd ${params.root} && bash 03_extract_transfer-learned_instruments/03.5_create_instruments_list.sh
  touch done.trigger
  """
}

process create_unique_instruments_list_tl {
  tag "cftlmr"
  input: path('trigger')
  output: path('done.trigger'), optional: true
  script:
  """
  ${env_graph_mr}
  export METHOD_TYPE="cftlmr"
  cd ${params.root} && bash 03_extract_transfer-learned_instruments/03.6_create_unique_instruments_list.sh
  touch done.trigger
  """
}

process extract_dosages_tl {
  tag "cftlmr"
  input: path('trigger')
  output: path('done.trigger'), optional: true
  script:
  """
  export METHOD_TYPE="cftlmr"
  export INPUT_TYPE="${params.dataset}"
  cd ${params.root} && bash 03_extract_transfer-learned_instruments/03.7_extract_dosages.sh
  touch done.trigger
  """
}

process create_transfer_learned_scores {
  tag "cftlmr"
  input: path('trigger'), val(batch_id)
  output: path('done.trigger'), optional: true
  script:
  """
  ${env_graph_mr}
  export METHOD_TYPE="cftlmr"
  export PBS_ARRAY_INDEX="$batch_id"
  cd ${params.root} && bash 03_extract_transfer-learned_instruments/03.8_create_transfer-learned_scores.sh
  touch done.trigger
  """
}

// ---------------------------------------------------------------------------
// Method-specific: concatenate scores, total effect models, format/visualise
// ---------------------------------------------------------------------------

process join_cftlmr_triggers {
  tag "cftlmr"
  input: val(cf_done), val(tl_done)
  output: path('done.trigger'), optional: true
  script: "touch done.trigger"
}

process emit_trigger {
  tag "trigger"
  input: val(any)
  output: path('done.trigger'), optional: true
  script: "touch done.trigger"
}

process concatenate_scores {
  tag "${method_type}"
  input: path('trigger'), val(method_type)
  output: path('done.trigger'), optional: true
  script:
  """
  export METHOD_TYPE="$method_type"
  export INPUT_TYPE="${params.dataset}"
  cd ${params.root} && bash 04_concatenate_scores.sh
  touch done.trigger
  """
}

process concatenate_scores_cftlmr {
  tag "cftlmr"
  input: path('trigger')
  output: path('done.trigger'), optional: true
  script:
  """
  ${env_graph_mr}
  export METHOD_TYPE="cftlmr"
  cd ${params.root} && bash 04_concatenate_scores.sh
  touch done.trigger
  """
}

process concatenate_scores_cfmr {
  tag "cfmr"
  input: path('trigger')
  output: path('done.trigger'), optional: true
  script:
  """
  export METHOD_TYPE="cfmr"
  export INPUT_TYPE="${params.dataset}"
  cd ${params.root} && bash 04_concatenate_scores.sh
  touch done.trigger
  """
}

process run_total_effect_models {
  tag "${method_type}"
  input: path('trigger'), val(method_type), val(task_id)
  output: path('done.trigger'), optional: true
  script:
  """
  ${env_graph_mr}
  export METHOD_TYPE="$method_type"
  export PBS_ARRAY_INDEX="$task_id"
  cd ${params.root} && bash 05.1_run_total_effect_models.sh
  touch done.trigger
  """
}

process collate_total_effect_estimates {
  tag "${method_type}"
  input: path('trigger'), val(method_type)
  output: path('done.trigger'), optional: true
  script:
  """
  ${env_graph_mr}
  export METHOD_TYPE="$method_type"
  cd ${params.root} && bash 05.2_collate_total_effect_estimates.sh
  touch done.trigger
  """
}

process collate_total_effect_estimates_cftlmr {
  tag "cftlmr"
  input: path('trigger')
  output: path('done.trigger'), optional: true
  script:
  """
  export METHOD_TYPE="cftlmr"
  export INPUT_TYPE="${params.dataset}"
  cd ${params.root} && bash 05.2_collate_total_effect_estimates.sh
  touch done.trigger
  """
}

process collate_total_effect_estimates_cfmr {
  tag "cfmr"
  input: path('trigger')
  output: path('done.trigger'), optional: true
  script:
  """
  ${env_graph_mr}
  export METHOD_TYPE="cfmr"
  cd ${params.root} && bash 05.2_collate_total_effect_estimates.sh
  touch done.trigger
  """
}

process format_visualise_data {
  tag "cftlmr"
  input: path('trigger')
  output: path('done.trigger'), optional: true
  script:
  """
  export METHOD_TYPE="cftlmr"
  export INPUT_TYPE="${params.dataset}"
  cd ${params.root} && bash 06_07_format_visualise_data.sh
  touch done.trigger
  """
}

process format_visualise_data_cfmr {
  tag "cfmr"
  input: path('trigger')
  output: path('done.trigger'), optional: true
  script:
  """
  ${env_graph_mr}
  export METHOD_TYPE="cfmr"
  cd ${params.root} && bash 06_07_format_visualise_data.sh
  touch done.trigger
  """
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

workflow {
  // --- Common: 00 -> 02.1 -> 02.2 (run 01_preprocess_exposures_generic.sh first if needed)
  clean_plink_variants()
  split_into_folds(clean_plink_variants.out)
  create_results_directories(split_into_folds.out)

  // --- Get N after 02.2
  get_exposure_count(create_results_directories.out)

  def n_val = get_exposure_count.out.map { it.trim().toInteger() }
  n_val.view { "N (exposures) = $it" }

  // --- Cross-fitting: 02.3.1–02.3.3 in parallel, then 02.4
  run_local_PLINK_GWAS_1(n_val)
  run_local_PLINK_GWAS_2(n_val)
  run_local_PLINK_GWAS_3(n_val)
  run_PLINK_clumping(run_local_PLINK_GWAS_1.out, run_local_PLINK_GWAS_2.out, run_local_PLINK_GWAS_3.out)

  // --- 02.5 array (1..N): each task gets PBS_ARRAY_INDEX
  def idx_ch = n_val.flatMap { n -> channel(1..n) }
  load_clean_PLINK_results(run_PLINK_clumping.out.combine(idx_ch))

  // --- 02.6 -> 02.7 -> 02.8
  def load_done = load_clean_PLINK_results.out.collect()
  create_instruments_list(emit_trigger(load_done).out)
  create_unique_instruments_list(create_instruments_list.out)
  extract_dosages(create_unique_instruments_list.out)

  // --- 02.9 array: A = ceil((N-1)/10)
  def batch_a_ch = n_val.flatMap { n ->
    def a = (int)((n - 1) / 10) + (((n - 1) % 10 > 0) ? 1 : 0)
    channel(1..a)
  }
  create_cross_fitted_scores(extract_dosages.out.combine(batch_a_ch))

  // --- Transfer learning (cftlmr): 03.1 -> ... -> 03.8
  create_missing_complete_lists()
  run_local_PLINK_GWAS_tl(create_missing_complete_lists.out, n_val)
  run_PLINK_clumping_tl(run_local_PLINK_GWAS_tl.out)
  load_clean_PLINK_results_tl(run_PLINK_clumping_tl.out)
  create_instruments_list_tl(load_clean_PLINK_results_tl.out)
  create_unique_instruments_list_tl(create_instruments_list_tl.out)
  extract_dosages_tl(create_unique_instruments_list_tl.out)
  def batch_b_ch = n_val.flatMap { n ->
    def b = (int)((n - 1) / 50) + (((n - 1) % 50 > 0) ? 1 : 0)
    channel(1..b)
  }
  create_transfer_learned_scores(extract_dosages_tl.out.combine(batch_b_ch))

  // --- CFTLMR: 04 (after both 02.9 and 03.8) -> 05.1 (array) -> 05.2 -> 06_07
  def cf_done = create_cross_fitted_scores.out.collect()
  def tl_done = create_transfer_learned_scores.out.collect()
  join_cftlmr_triggers(cf_done.combine(tl_done))
  concatenate_scores_cftlmr(join_cftlmr_triggers.out)
  def task_idx = n_val.flatMap { n -> channel(1..n) }
  def cftlmr_05_ch = concatenate_scores_cftlmr.out.combine(channel.of('cftlmr').combine(task_idx)).map { (t, mt_id) -> tuple(t, mt_id[0], mt_id[1]) }
  run_total_effect_models(cftlmr_05_ch)
  // CFMR 04 and 05.1
  concatenate_scores_cfmr(emit_trigger(cf_done).out)
  def cfmr_05_ch = concatenate_scores_cfmr.out.combine(channel.of('cfmr').combine(task_idx)).map { (t, mt_id) -> tuple(t, mt_id[0], mt_id[1]) }
  run_total_effect_models(cfmr_05_ch)
  // Both methods' 05.1 complete; then 05.2 (collate) and 06_07 (format/visualise) for each
  def all_05_done = run_total_effect_models.out.collect()
  emit_trigger(all_05_done).out.into { collate_cftlmr_trigger; collate_cfmr_trigger }
  collate_total_effect_estimates_cftlmr(collate_cftlmr_trigger)
  collate_total_effect_estimates_cfmr(collate_cfmr_trigger)
  format_visualise_data(collate_total_effect_estimates_cftlmr.out)
  format_visualise_data_cfmr(collate_total_effect_estimates_cfmr.out)
}