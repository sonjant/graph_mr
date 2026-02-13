# Graph MR: generic input specification (plug-and-play)

Any user with **genotype data** and **one post-genotype omic dataset** (e.g. proteomics, metabolomics, transcriptomics, phenomics) can run Graph MR. Dataset and study names are not hardcoded.

## What you need

1. **Genotype data** (PLINK format + optional VCF for dosages), and principal components (e.g. `merged.pruned.eigenvec`).
2. **One exposure dataset**: a matrix with **sample IDs** (matching genotype) and **exposure columns** (one column per exposure, e.g. protein, metabolite, or phenotype).
3. **Covariates**: age and sex for each sample (and optionally other covariates). Same sample ID as in the exposure matrix.

## Standard directory layout (after preprocessing)

The pipeline expects a **single dataset** (any name you choose, e.g. `proteomics`, `mydata`, `airwave_metabolomics`) under a data directory. You can either:

- **Option A**: Run the **generic preprocessing** script (see below) on your raw exposure matrix + covariates; it will create this layout, or  
- **Option B**: Create the layout yourself.

Layout:

```
<data_dir>/<dataset>/
├── preprocessed/
│   ├── exposures.tsv          # Required: sample_id column + one column per exposure (standardized recommended)
│   └── exposure_metadata.tsv  # Optional: exposure_id, display_name [, pathway/class]
└── gwas_input/
    ├── gwas_covariates_pc1-5_age_sex.txt   # Required for GWAS (created by preprocessing)
    └── covariates_age_sex.txt             # Required for MR step: columns barcode, AGE, SEX
```

- **`exposures.tsv`**: First column must be the **sample ID** used in genotype (e.g. `barcode` or `IID`). The pipeline uses the name `barcode` internally; preprocessing can rename your ID column to `barcode`. Remaining columns = exposure variables (one per protein/metabolite/phenotype/etc.).
- **`exposure_metadata.tsv`** (optional): Maps exposure IDs to human-readable names for figures. Columns: `exposure_id` (matches column names in exposures.tsv), `display_name`; optional: `pathway`, `class`, etc.
- **`gwas_covariates_pc1-5_age_sex.txt`**: Used by PLINK GWAS. Must include sample ID (FID, IID or barcode), PCs 1–5, AGE, SEX. Format expected by your GWAS step (see preprocessing script).
- **`covariates_age_sex.txt`**: Tab- or comma-separated: `barcode`, `AGE`, `SEX` (SEX: 1/2 or M/F). Used in the total-effect MR step (05.1).

## Generic preprocessing script

Use `01_preprocess_data/01_preprocess_exposures_generic.R` to go from your files to the layout above.

**Inputs (set via environment variables or arguments):**

- **Exposure matrix**: One table with sample ID column and exposure columns (e.g. protein levels, metabolite levels). No missing values in sample ID.
- **Covariate file**: Sample ID, age, sex. Same sample ID as in the exposure matrix and as in genotype.
- **Genotype PCs**: e.g. from PLINK `--pca`; file with sample ID and PC1, PC2, … (e.g. `merged.pruned.eigenvec`).

**What the script does:**

- Merges exposure matrix, covariates, and PCs on sample ID.
- Filters to samples present in genotype (PCs file).
- Writes standardized (or optional) exposure matrix → `preprocessed/exposures.tsv` (with sample ID column named `barcode`).
- Builds GWAS covariate file (PCs + age + sex) → `gwas_input/gwas_covariates_pc1-5_age_sex.txt`.
- Writes `gwas_input/covariates_age_sex.txt` (barcode, AGE, SEX) for the MR step.
- Optionally writes a minimal `exposure_metadata.tsv` (exposure_id = column name, display_name = column name) if you don’t provide one.

**No study-specific logic**: no references to “airwave”, “metabolomics”, “proteomics”, “tripleintersection”, or “phenomolecular”. The **dataset** is just a label (e.g. `--dataset proteomics`) used in paths.

## Running the pipeline

1. Set **GRAPH_MR_ROOT**, **GRAPH_MR_GENOTYPE_DIR**, **GRAPH_MR_VCF_DIR** (for dosage steps) as in README.
2. Set **GRAPH_MR_DATASET** to your dataset name (e.g. `proteomics`, `mydata`). This is the only “input type”: a label for the folder under `data/`.
3. Either:
   - Run generic preprocessing so that `data/<GRAPH_MR_DATASET>/preprocessed/exposures.tsv` and `data/<GRAPH_MR_DATASET>/gwas_input/` exist, or
   - Create that structure yourself.
4. Run the pipeline (or Nextflow) with `dataset=<your_dataset>` (and method type as usual). All steps will read/write using `data/<dataset>/` and `results/<method>_<dataset>/`; no hardcoded dataset names.

## Summary

- **One post-genotype dataset** = one folder under `data/<dataset>/` with the layout above.
- **Dataset name** = arbitrary label (e.g. proteomics, transcriptomics, my_study).
- **Genotype + one such dataset** is enough to run Graph MR in a plug-and-play way; no need for “airwave”, “tripleintersection”, or “phenomolecular” in the code or config.
