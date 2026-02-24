# FusAphid_results

## Overview

This repository contains a complete MOFA2-based multi-omics analysis workflow for aphid and honeydew samples from an NIV exposure study. It includes:

- the input workbook used to drive analysis
- the full analysis script (`Input/MOFA-CORE.R`)
- annotation caches and derived annotation tables
- a trained MOFA model (`.hdf5`)
- global and per-factor downstream outputs (plots and tables)

The repo is currently in a "results snapshot + rerunnable script" state rather than a minimal source-only package.

## Current repository state

### Top-level contents

- `Input/`: workbook + scripts
- `Annotations/`: generated/cached KEGG, lipid, and FELLA resources
- `Results/`: generated model and analysis outputs
- `renv/`: renv bootstrap files (no committed `renv.lock`)
- `MOFA-Aphid.R`: legacy wrapper script in repo root

### Active analysis views in the current workbook

Source: `Input/Aphids and honeydew from exposure to niv MOFA.xlsx`, sheet `view_metadata`.

| View | Type | Kegg_Species | View_Group | Active | Features (rows) | Samples (columns excluding `Name`) |
| --- | --- | --- | --- | --- | ---: | ---: |
| LCMS-Aphid | Metabolic | api | Aphid | Y | 617 | 24 |
| LCMS-BuchHam | Metabolic | buc,hde | Aphid | Y | 99 | 24 |
| LCMS-Ham | Metabolic | hde | Aphid | Y | 2 | 24 |
| LCMS-Buch | Metabolic | buc | Aphid | Y | 35 | 24 |
| LCMS-Honeydew | Metabolic | api | Honeydew | Y | 512 | 24 |
| OrbiSIMS-MetaProt | OrbiSIMS | api | Honeydew | Y | 13 | 24 |
| OrbiSIMS-MetaDeProt | OrbiSIMS | api | Honeydew | Y | 20 | 24 |

There are currently no active `RNAseq` views in this workbook snapshot.

### Sample grouping snapshot

Source: workbook sheet `sample_metadata`.

- Total samples: 24
- `category` groups: `Mock` (6), `DONwt` (6), `NIVwt` (6), `NIVko` (6)
- `group`: `group1` for all 24 samples

### Existing generated outputs currently present

- `Results/Full_MOFA_model.hdf5` exists
- `Results/Per_Factor_Analysis` contains `Factor1` to `Factor4`
- `Annotations/` already contains:
  - `KEGG_gene_annotations.csv` (legacy combined naming)
  - `KEGG_metabolite_annotations.csv`
  - `KEGG_metabolite_annotations_api.csv`
  - `KEGG_metabolite_annotations_buc.csv`
  - `KEGG_metabolite_annotations_hde.csv`
  - `Lipid.annotations.csv`
- KEGG cache species currently present under `Annotations/kegg_annotation_cache/`:
  - `api`, `buc`, `hde`, and `taes` pathway cache files
- FELLA database folders currently present under `Annotations/fella_db/`:
  - `api`, `buc`, `hde`, `taes` (with graph folders for `api`, `buc`, `hde`)
- Result file counts by area:
  - `Corrilation_and_Variance_plots`: 6 files
  - `Covariate_factor_comparisons`: 9 files
  - `PCA`: 18 files
  - `Significance_testing`: 13 files
  - `Per_Factor_Analysis`: 185 files

## Execution entry points

### Recommended entry point

Use `Input/MOFA-CORE.R` directly.

Run from repository root:

```bash
Rscript -e 'input_file <- "Input/Aphids and honeydew from exposure to niv MOFA.xlsx"; OLD_NORMALISATION <- TRUE; source("Input/MOFA-CORE.R")'
```

### Legacy wrappers (not recommended without edits)

Both of the following contain hard-coded legacy paths outside this repository:

- `MOFA-Aphid.R` (repo root)
- `Input/MOFA-Aphid.R`

They currently include:

- hard-coded `dir=".../Rumiana_Ray_NIV_paper/Aphid"`
- `source("../Scripts/MOFA-CORE.R")`

So they will not run as-is in this repo unless path values are updated.

## Script workflow in detail

The `Input/MOFA-CORE.R` script defines functions first, then executes the full workflow near the bottom.

### Stage 1: Environment setup

Function:

- `setup_project_environment(use_renv = TRUE, renv_path = "renv")`

Execution call in this repo:

- `setup_project_environment(use_renv = FALSE)`

Behavior:

- installs missing CRAN/Bioconductor packages
- loads all required libraries into the session
- can optionally initialize/restore renv if enabled

### Stage 2: Metadata and input loading

Functions:

- `import_view_metadata(input_file)`
- `import_sample_metadata(input_file, view_metadata)`
- `load_and_clean_input_data(input_file, view_metadata)`

Behavior:

- reads `view_metadata` and `sample_metadata`
- reads all sheets where `view_metadata$Active == "Y"`
- cleans `Name` values:
  - replaces Greek letters with text equivalents (`alpha`, `beta`, etc.)
  - strips `%` and `,`
  - transliterates to ASCII (`iconv(..., "ASCII//TRANSLIT")`)
- moves `Name` to row names
- aligns sample column order across views to match the first active view

### Stage 3: Annotation and cache construction/loading

Functions:

- `generate_lipid_annotations(...)`
- `generate_kegg_annotation_lists(...)`
- `annotate_transcriptomic_data_with_kegg(...)`
- `annotate_metabolomic_data_with_kegg(...)`
- `load_or_build_fella_database(...)`

Behavior:

- lipid annotations:
  - uses `rgoslin::parseLipidNames`
  - calls RefMet REST API endpoint for extra lipid fields
  - writes/loads `Annotations/Lipid.annotations.csv`
- KEGG mapping:
  - generates or loads cached species pathway lists from `Annotations/kegg_annotation_cache/*_pathways.rds`
  - builds `TERM2GENE` and `TERM2NAME` for genes and compounds per species
- transcriptomic annotation:
  - merges feature IDs to KEGG gene IDs
  - writes species files `Annotations/KEGG_gene_annotations_<species>.csv`
  - in this repository snapshot, there are no active transcriptomic (`RNAseq`) views, so this step may be skipped
- metabolite annotation:
  - cleans compound names (uses substring after `|`)
  - queries `KEGGREST::keggFind("compound", ...)`
  - writes `Annotations/KEGG_metabolite_annotations_<species>.csv`
- FELLA:
  - loads or builds per-species KEGG graph and diffusion matrix
  - stored under `Annotations/fella_db/<species>/graph/`

### Stage 4: Normalization (two modes)

Selection is controlled by global `OLD_NORMALISATION`.

- `TRUE` -> `normalise_data_mofa_v1(...)`
- `FALSE` -> `normalise_data_mofa_v2(...)`

#### v1 (`normalise_data_mofa_v1`)

Defaults:

- `na_limit = 1e-7`

Behavior summary:

- for views starting `RNASeq` or `LCMS`, near-zero values become `NA`
- other views:
  - negative values become `NA`
  - positive values offset by `0.1 * min_nonzero`
- `log2` transform
- row-wise z-score
- removes duplicate row names and all-zero/all-NA rows
- appends `_<view_name>` suffix to feature IDs

#### v2 (`normalise_data_mofa_v2`)

Defaults:

- `group_col = "category"`
- `top_n = 1000`
- `noise_sd = 1e-8`
- `noise_seed = 1`
- `retry_with_noise = TRUE`

Behavior summary:

- filters invalid/zero rows, applies `log2(x + 1)`
- supervised ranking with limma moderated F if possible
- falls back to row variance if limma is not viable
- optional tiny-noise jitter for constant rows
- keeps top `top_n` features per view
- row-wise z-score and feature suffixing with `_<view_name>`

### Stage 5: MOFA model initialization/training/loading

Function:

- `initialise_and_run_mofa(mofa_list, output_file = "Results/Full_MOFA_model.hdf5")`

Key training options set in code:

- max factors capped at 6 (`model_opts$num_factors <= 6`)
- `train_opts$drop_factor_threshold = 0.01`
- `train_opts$convergence_mode = "slow"`
- `train_opts$maxiter = 10000`
- training uses `run_mofa(..., use_basilisk = TRUE)`

Model file behavior:

- if `Results/Full_MOFA_model.hdf5` exists, script loads it
- otherwise it trains and writes that file

### Stage 6: Metadata assignment and color maps

Functions:

- `generate_sample_metadata(...)`
- `assign_sample_metadata_to_mofa(...)`
- `define_group_colours(...)`
- `define_view_colours(...)`

Notes:

- if sheet-derived `sample_metadata` is available, it is converted to factors
- metadata columns considered numeric-like text are excluded from color mapping

### Stage 7: MOFA post-processing and statistics

Functions:

- `extract_kegg_weights(...)`
- `run_limma_statistics(...)`
- `standardise_feature_names(...)`

Important behavior:

- `extract_kegg_weights` creates KEGG-only per-view/per-species weight matrices in:
  - `MOFAobject.trained@expectations$KEGG_weights`
- `run_limma_statistics` runs per-view, per-metadata-column limma tests and keeps rows with:
  - FDR (`adj.P.Val`) <= 0.05

### Stage 8: Global visual outputs

Functions and output roots:

- `generate_mofa_variance_plots(...)`
  - `Results/Corrilation_and_Variance_plots/Variance_explained.png`
- `generate_mofa_factor_matrix_plot(...)`
  - `Results/PCA/`
- `generate_mofa_factor_comparison_plot(...)`
  - `Results/Covariate_factor_comparisons/`
- `generate_mofa_factor_correlation_plot(...)`
  - `Results/Corrilation_and_Variance_plots/`

### Stage 9: Significance/regression diagnostics

Functions:

- `save_weight_tables(...)`
- `plot_scatter_pval_vs_weight(...)`
- `plot_violin_pval_vs_weight(...)`

Outputs:

- `Results/Significance_testing/tables/Weights_AllViews_AllFactors__<meta_col>.csv`
- `Results/Significance_testing/regression_testing/...`
- `Results/Significance_testing/violin_plots/...`

### Stage 10: Per-factor outputs

Functions:

- `plot_top_weights_heatmap_faceted_by_view(...)`
- `plot_top_weights_view(factor_name)`
- `plot_top_weights_per_factor(factor_name, top_n = 10)`
- `plot_scatter_per_view(factor_name)`
- `plot_heatmap_per_view(factor_name, top_features = 25)`
- `create_weight_outputs(factor_name, view_metadata)`

Output root:

- `Results/Per_Factor_Analysis/Factor*/...`

### Stage 11: KEGG enrichment and network/pathway outputs

Functions:

- `perform_kegg_enrichment(...)`
- `plot_enrichment_network(...)`
- `generate_simple_pathway_graph(...)`

Key thresholds/parameters from code:

- skip cluster if KEGG ID overlap < 10
- GSEA: `pvalueCutoff = 0.05`, `pAdjustMethod = "none"`, `minGSSize = 10`, `maxGSSize = 500`, `nPermSimple = 10000`
- pathway plots generated for enrichment rows with `pvalue < 0.05`
- enrichment network keeps top 100 edges by absolute enrichment score and uses `pvalue < 0.05`

### Stage 12: FELLA enrichment

Function:

- `run_fella_enrichment(...)`

Behavior summary:

- runs per metabolomic view and per species listed in `Kegg_Species`
- writes enrichment table:
  - `Results/Per_Factor_Analysis/<Factor>/<View>/FELLA_Enrichment_<species>.csv`
- writes diffusion graph image:
  - `Results/Per_Factor_Analysis/<Factor>/<View>/FELLA_enrichment_graph_top50_<species>.png`

## Input workbook contract

The pipeline assumes these workbook conventions.

### Required sheets

1. `view_metadata`
2. all active data views named in `view_metadata$View` where `Active == "Y"`
3. `sample_metadata` (required for the current code path)

### `view_metadata` required columns

- `View`
- `Type`
- `Kegg_Species`
- `Active`

Columns also used/retained in current file:

- `Pathview_Species`
- `View_Group`

### Data view sheet shape

- first column must be `Name`
- remaining columns are sample IDs
- sample columns should be in the same order across active views (script enforces this by reordering)

### `sample_metadata` expected orientation

The script expects:

- first column named `Name` listing metadata variable names
- sample IDs as remaining columns
- then transposes this to sample rows x covariate columns

If this orientation is changed, `import_sample_metadata()` and downstream metadata usage must be updated.

## Output map (directory-level reference)

### `Results/Corrilation_and_Variance_plots/`

- `Variance_explained.png`
- `Factor_vs_Sample_corrilation.png`
- `Factor_vs_Covariate_corrilation.png`
- `Factor_vs_Covariate_corrilation_pval.png`
- `Factor_vs_Covariate_corrilation_pval.csv`
- `Factor_vs_Factor_corrilation.png`

### `Results/PCA/`

- combined matrix plot:
  - `<meta_col>_Factor_Pair_PCA_Matrix.png`
- individual upper-triangle sample PCA panels:
  - `individual_data_loading_plots/<meta_col>_<i>_vs_<j>.png`
- individual lower-triangle feature-weight panels:
  - `individual_sample_PCA_plots/<meta_col>_<i>_vs_<j>_weights.png`
- diagonal density plots:
  - `individual_histograms/<meta_col><i>_density.png`

### `Results/Covariate_factor_comparisons/`

- combined per-covariate panels:
  - `<meta_col>_Factor_Comparison.png`
- per-factor violin plots:
  - `individual_violin_plots/<meta_col>_Factor<k>_violin.png`
- per-factor mean-vs-weight scatter:
  - `individual_data_loading_plots/<meta_col>_Factor<k>_mean_vs_weight.png`

### `Results/Significance_testing/`

- filtered table exports:
  - `tables/Weights_AllViews_AllFactors__<meta_col>.csv`
- regression testing:
  - `regression_testing/R2_scatter__<meta_col>.png`
  - `regression_testing/individual_R2_barplots/R2__<meta_col>_<pcol>_Factor<k>.png`
  - `regression_testing/individual_R2_scatter_plots/Scatter__<meta_col>_<pcol>_Factor<k>.png`
- violin summaries:
  - `violin_plots/Violin__<meta_col>.png`

### `Results/Per_Factor_Analysis/Factor*/`

Global per-factor files:

- `All_Weights.csv`
- `Top-weights-bar.png`
- `Top_100_edges_KEGG_cnet.png`
- `KEGG_Enrichment_Weighted.csv`
- `KEGG_Enrichment.png`

Per-view subfolders (examples: `LCMS-Aphid`, `OrbiSIMS-MetaProt`) typically include:

- `Top_Weights.png`
- `Heatmap.png`
- `scatterplots/<meta_col>_Scatter.png`
- `FELLA_Enrichment_<species>.csv` (metabolic views)
- `FELLA_enrichment_graph_top50_<species>.png` (metabolic views)
- `KEGG_pathway_plots/*.png` (when significant pathways exist)

## Configuration knobs you can change safely

Common run-time controls:

- `input_file`: workbook path
- `OLD_NORMALISATION`: `TRUE` (v1) or `FALSE` (v2)

Function defaults worth tuning:

- `normalise_data_mofa_v2(top_n = 1000, group_col = "category")`
- `initialise_and_run_mofa(output_file = "Results/Full_MOFA_model.hdf5")`
- `plot_top_weights_per_factor(top_n = 10)`
- `plot_heatmap_per_view(top_features = 25)`

To force retraining only:

- remove or rename `Results/Full_MOFA_model.hdf5` before running

To rebuild annotations/caches:

- remove relevant files under `Annotations/` before running (script will regenerate where possible)

## External services and network usage

Depending on what is already cached, the script may call:

- KEGG REST (`KEGGREST`) for pathway and compound mappings
- RefMet REST endpoint at `metabolomicsworkbench.org` for lipid annotations

If caches already exist in `Annotations/`, those are reused and network calls may be reduced.

## Dependency details

`setup_project_environment()` installs/loads packages at runtime.

### CRAN packages listed in script

`data.table`, `ggplot2`, `psych`, `dplyr`, `ggpubr`, `GGally`, `readxl`, `readr`, `tidyr`, `stringi`, `openxlsx`, `reshape2`, `writexl`, `textclean`, `stringr`, `viridis`, `grid`, `gtable`, `gridExtra`, `gridGraphics`, `stats`, `clusterProfiler`, `pathview`, `magick`, `tibble`, `tidygraph`, `ggraph`, `igraph`, `scales`, `ggtext`, `purrr`, `XML`, `png`, `httr`, `httr2`, `jsonlite`, `VIM`, `preprocessCore`, `matrixStats`, `ggforce`

### Bioconductor packages listed in script

`tximport`, `DESeq2`, `MOFA2`, `KEGGREST`, `MSnbase`, `clusterProfiler`, `ggkegg`, `rgoslin`, `fgsea`, `sva`, `limma`, `FELLA`, `biomaRt`, `ComplexHeatmap`

## Reproducibility notes

- `renv/` exists but no `renv.lock` is committed in this repository, so exact package versions are not pinned here.
- several outputs are conditionally skipped if files already exist (for example `Variance_explained.png`) or if data is insufficient (for example low-feature views for some plots).
- some analyses depend on mutable external resources (KEGG/RefMet).

## Known caveats in the current codebase

- wrapper scripts use hard-coded legacy paths and should be treated as templates.
- output folder/file names intentionally include historical spelling (`Corrilation`) because that is what the script writes.
- workbook has inactive view-name variants with hyphen/underscore differences; this does not currently break execution because those rows are inactive.
- the script is stateful and relies on shared objects (`metadata`, `group.colors`, etc.) created in execution order, so calling functions ad hoc requires setting those objects first.

## Practical rerun modes

### Fast rerun using existing model

```bash
Rscript -e 'input_file <- "Input/Aphids and honeydew from exposure to niv MOFA.xlsx"; OLD_NORMALISATION <- TRUE; source("Input/MOFA-CORE.R")'
```

### Rerun with v2 normalization

```bash
Rscript -e 'input_file <- "Input/Aphids and honeydew from exposure to niv MOFA.xlsx"; OLD_NORMALISATION <- FALSE; source("Input/MOFA-CORE.R")'
```

### Retrain MOFA model

Rename or remove `Results/Full_MOFA_model.hdf5`, then run one of the commands above.

## Development status

- This repository currently focuses on scripted analysis execution and generated outputs.
- There are no automated tests, CI workflows, or packaged command-line interfaces in the current tree.
