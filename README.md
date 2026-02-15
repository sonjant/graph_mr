# Graph MR pipeline

![Graph MR logo](docs/Colourful-Gradient.svg)  
*Logo by [Wei-Ting Chia](https://chiaweiting.cargo.site/About).*

Mendelian randomization pipeline for multi-omics (e.g. proteomics, metabolomics, transcriptomics) with cross-fitting and optional transfer learning. **Plug-and-play**: you need **genotype data** and **one post-genotype omic dataset**; no study-specific names (e.g. no "airwave" or "tripleintersection") are hardcoded. Suitable for HPC clusters using PBS; also runnable via Nextflow.

## Plug-and-play: what you need

1. **Genotype data**: PLINK format (+ VCF for dosage steps), and PCs (e.g. `merged.pruned.eigenvec`).
2. **One exposure dataset**: a matrix with sample IDs (matching genotype) and exposure columns (e.g. proteins, metabolites, phenotypes).
3. **Covariates**: age and sex per sample (same sample ID as exposures and genotype).

See **`docs/INPUT_SPEC.md`** for the exact input layout and the **generic preprocessing** script that converts your files into the pipeline format.

## Configuration (required before running)

Paths and environment are controlled by **environment variables** so the pipeline is portable and shareable. No user-specific paths are hardcoded.

### 1. Set environment variables

Copy the example config and edit with your paths:

```bash
cp config/config.example.sh config/config.sh
# Edit config/config.sh with your paths
source config/config.sh
```

**Required:**

| Variable | Description |
|----------|-------------|
| `GRAPH_MR_ROOT` | Full path to this pipeline directory (the repo root). |
| `GRAPH_MR_GENOTYPE_DIR` | Directory containing genotype/PLINK data: `merged.pruned.eigenvec` and per-chromosome files e.g. `merged.chr1.1kg.bed`/`.bim`/`.fam`. |
| `GRAPH_MR_VCF_DIR` | Directory containing VCF files for dosage extraction (e.g. `merged.chr1.1kg.vcf.gz`). Required for steps 02.8 and 03.7. |

**Optional:**

| Variable | Description | Default |
|----------|-------------|--------|
| `GRAPH_MR_EPHEMERAL_DIR` | Base directory for large intermediate outputs (GWAS, dosages). | `$GRAPH_MR_ROOT/ephemeral` |
| `GRAPH_MR_RESULTS_DIR` | Directory for final pipeline results. | `$GRAPH_MR_ROOT/results` |
| `GRAPH_MR_DATASET` | Your dataset name (e.g. `proteomics`, `mydata`). Used as folder under `data/` and in results. | Set per run or use `INPUT_TYPE` in scripts. |
| `GRAPH_MR_DATA_DIR` | Directory containing one folder per dataset. | `$GRAPH_MR_ROOT/data` |
| `CONDA_HOME` | Path to conda installation. | `$HOME/miniforge3` |
| `CONDA_ENV` | Conda environment name (must have R and required packages). | `r443` |

Add `config/config.sh` to `.gitignore` so your local paths are not committed.

### 2. Directory layout

- **Config**: `config/config.R` (sourced by R scripts) and `config/config.example.sh` (copy to `config/config.sh`, edit, then source before running).
- **Data**: Under `$GRAPH_MR_ROOT/data/` (or `GRAPH_MR_DATA_DIR`) create **one folder per dataset** (any name you choose, e.g. `data/proteomics/`, `data/mydata/`). Each folder must follow the layout in `docs/INPUT_SPEC.md` (e.g. `preprocessed/exposures.tsv`, `gwas_input/covariates_age_sex.txt`). Use the **generic preprocessing** script to generate this from your exposure matrix + covariates.
- **Results**: By default, `results/` and `ephemeral/` are under `GRAPH_MR_ROOT`. Override with `GRAPH_MR_RESULTS_DIR` and `GRAPH_MR_EPHEMERAL_DIR` if you want them elsewhere.

### 3. Running the pipeline

**Option A: Shell (PBS or local)**  
Set `INPUT_TYPE` to your **dataset name** (any string, e.g. `proteomics`, `mydata`). This is the only "input type": a label for the folder under `data/`.

```bash
source config.sh
export METHOD_TYPE="cftlmr"
export INPUT_TYPE="proteomics"   # or your dataset name
export METHOD_TYPES="cftlmr,cfmr"

# If you haven't prepared data yet: run generic preprocessing once (set GRAPH_MR_EXPOSURE_MATRIX, GRAPH_MR_COVARIATE_FILE)
./01_preprocess_data/01_preprocess_exposures_generic.sh

# Then run the pipeline (e.g. from 00_clean_PLINK_variants, then 02.1_split_into_folds, etc.)
./01_preprocess_data/00_clean_PLINK_variants.sh
./02_extract_cross-fitted_instruments/02.1_split_into_folds.sh
# ... remaining steps
```

**Option B: Nextflow**  
Run the pipeline with one command; Nextflow sets env vars from params and runs the first stages (preprocess optional → clean PLINK → split folds → create dirs). See *Nextflow* below.

R scripts source `config/config.R` from the pipeline root and expect **working directory =** `GRAPH_MR_ROOT`. Shell scripts `cd` to `GRAPH_MR_ROOT` before calling R.

### Nextflow (plug-and-play)

A single Nextflow pipeline lives in **`nextflow/`**. From the pipeline root run:

```bash
nextflow run nextflow/main.nf -c nextflow/nextflow.config --root "$(pwd)" --genotype_dir /path/to/plink --vcf_dir /path/to/vcf --dataset mydata
```

- **Required params**: `--root`, `--genotype_dir`, `--vcf_dir`, `--dataset` (your dataset name).
- **Optional**: `--data_dir`, `--ephemeral_dir`, `--results_dir`; `--method_types 'cftlmr,cfmr'`.
- **Data**: Prepare the standard layout first (run `01_preprocess_exposures_generic.sh` or create it manually; see `docs/INPUT_SPEC.md`). Params `--run_preprocess`, `--exposure_matrix`, `--covariate_file` are reserved for a future optional preprocess step.

Use a params file: `nextflow run nextflow/main.nf -c nextflow/nextflow.config -params-file nextflow/params.example.json` (edit the paths in the JSON). The workflow runs the full pipeline: 00 → 02.1 → 02.2 → GWAS → clumping → instruments → dosages → scores → 04 → 05.1/05.2 → 06_07.

## Requirements

- R with packages: tidyverse, vroom, readxl, data.table, parallel, gtools, and (for some steps) ensembldb, EnsDb.Hsapiens.v86, AnnotationFilter, GenomicRanges, igraph
- PLINK 2, bcftools (for genotype/dosage steps)
- Conda (or similar) to manage the R environment
- PBS (or run scripts directly without qsub)

## Licence and citation

Unless otherwise indicated, its contents are licensed under an Attribution-ShareAlike 4.0 International (CC BY-SA 4.0)