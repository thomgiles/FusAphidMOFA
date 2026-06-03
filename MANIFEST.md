# Repository Manifest

This manifest records the files that make this analysis deposit interpretable as a standalone computational artifact.

## Source Code

- `R/MOFA-CORE.R`: source-order aggregator for the function library.
- `R/001_*.R` to `R/038_*.R`: one file per function from the publication function library.
- `scripts/MOFA-Aphid.R`: executable publication runner.
- `tests/run_unit_tests.R`: lightweight unit test runner for function-level checks.
- `tests/check_repository.R`: lightweight validation script for the deposited files.

## Input Data

- `Input/Aphids and honeydew from exposure to niv MOFA.xlsx`: workbook containing `view_metadata`, `sample_metadata`, and the active analysis views.

Active workbook views in this deposit:

| View | Type | KEGG species | Group | Features | Samples |
| --- | --- | --- | --- | ---: | ---: |
| LCMS-Aphid | Metabolic | api | Aphid | 617 | 24 |
| LCMS-BuchHam | Metabolic | buc,hde | Aphid | 99 | 24 |
| LCMS-Ham | Metabolic | hde | Aphid | 2 | 24 |
| LCMS-Buch | Metabolic | buc | Aphid | 35 | 24 |
| LCMS-Honeydew | Metabolic | api | Honeydew | 512 | 24 |
| OrbiSIMS-MetaProt | OrbiSIMS | api | Honeydew | 13 | 24 |
| OrbiSIMS-MetaDeProt | OrbiSIMS | api | Honeydew | 20 | 24 |

## Annotation and Cache Files

- `Annotations/KEGG_gene_annotations.csv`
- `Annotations/KEGG_metabolite_annotations.csv`
- `Annotations/KEGG_metabolite_annotations_api.csv`
- `Annotations/KEGG_metabolite_annotations_buc.csv`
- `Annotations/KEGG_metabolite_annotations_hde.csv`
- `Annotations/Lipid.annotations.csv`
- `Annotations/kegg_annotation_cache/`: cached KEGG pathway resources.
- `Annotations/fella_db/`: cached FELLA graph and diffusion resources.

## Primary Results

- `Results/Full_MOFA_model.hdf5`: trained MOFA model.
- `Results/Corrilation_and_Variance_plots/`: variance explained and factor correlation outputs.
- `Results/PCA/`: factor pair plots and loadings/sample projections.
- `Results/Covariate_factor_comparisons/`: factor/covariate comparison figures.
- `Results/Significance_testing/`: significance and regression diagnostic outputs.
- `Results/Per_Factor_Analysis/`: per-factor feature, KEGG, and FELLA outputs.

This deposit currently contains 219 files under `Results/` and 49 files under `Annotations/`.
