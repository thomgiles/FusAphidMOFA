# FusAphidMOFA

This repository is the publication archive for the Fusarium-aphid multi-omics MOFA analysis. It contains the analysis code, the input workbook, curated annotation caches, and the deposited result tree used for the manuscript.

The manuscript workbook is the demonstration input for rerunning the published analysis.

The main workflow lives in `MOFA-CORE.R`. `MOFA-Aphid.R` is the top-level runner used to execute the analysis from the repository root.

## Repository Layout

| Path | Purpose |
| --- | --- |
| `MOFA-CORE.R` | Main analysis file containing the functions and the publication workflow. |
| `MOFA-Aphid.R` | Runner script that sources the core workflow. |
| `Input/Aphids and honeydew from exposure to niv MOFA.xlsx` | Analysis workbook. |
| `Annotations/` | KEGG, lipid, and FELLA annotation caches. |
| `Results/` | Deposited analysis outputs and trained model. |
| `renv/` | Project-local R environment support. |
| `.Rprofile` | Activates the local `renv` setup when the project is opened in R. |
| `DESCRIPTION` | Machine-readable package dependency list. |
| `MANIFEST.md` | File-level deposit inventory. |
| `LICENSE` | GNU GPL-3 license notice. |

Current snapshot sizes:

- `Annotations/`: 49 files
- `Results/`: 219 files

## System Requirements

Tested locally with:

- macOS / Darwin 25.3.0, arm64
- R 4.4.2 (2024-10-31)

The workflow uses CRAN and Bioconductor packages declared in `DESCRIPTION`. Key dependencies include `MOFA2`, `limma`, `DESeq2`, `clusterProfiler`, `KEGGREST`, `FELLA`, `ggplot2`, `readxl`, `dplyr`, `pathview`, and related plotting and data-processing packages.

No special hardware is required. GPU acceleration is not required. A normal desktop or laptop should be sufficient, although the full run can take a while and may use substantial memory during MOFA fitting.

The repository should include a final `renv.lock` once the current environment build is committed. Use that lockfile for exact dependency pinning before final archival.

## Installation

From a clean machine:

1. Install R 4.4.2 or later.
2. Clone this repository.
3. Open the project from the repository root so that `.Rprofile` can activate the local `renv` setup.
4. Run the analysis runner once to let the project environment resolve and load the required packages.

If you prefer to prepare the environment manually, install the packages listed in `DESCRIPTION` and then run the workflow from the repository root.

Typical setup time on a normal desktop is on the order of tens of minutes if packages must be installed or compiled from source.

## Running the Analysis

Run the publication workflow from the repository root:

```bash
Rscript MOFA-Aphid.R
```

The runner clears `./Results` before execution and then rebuilds the analysis outputs in place.

Important:

- run this from the repository root
- back up `Results/` first if you want to keep the deposited snapshot unchanged
- the workflow assumes the workbook is present at `Input/Aphids and honeydew from exposure to niv MOFA.xlsx`

## Input Workbook

The workbook is the analysis input used for the publication. It must contain:

- `view_metadata`
- `sample_metadata`
- one data sheet per active view listed in `view_metadata`

Expected structure:

- `view_metadata` columns include `View`, `Type`, `Kegg_Species`, and `Active`
- the `sample_metadata` sheet is used to assemble sample-level metadata
- each active data sheet should have `Name` as the first column
- the remaining columns should be sample IDs shared across the active views

## Outputs

The published outputs live under `Results/`. Key directories include:

- `Results/Full_MOFA_model.hdf5`
- `Results/Corrilation_and_Variance_plots/`
- `Results/PCA/`
- `Results/Covariate_factor_comparisons/`
- `Results/Significance_testing/`
- `Results/Per_Factor_Analysis/`

The exact output map is captured in `MANIFEST.md`.

## Reproducibility Notes

- `MOFA-CORE.R` contains the analysis implementation used for the publication archive.
- `MOFA-Aphid.R` is the top-level entrypoint for rerunning the workflow.
- `Annotations/` contains cached KEGG, lipid, and FELLA resources so the analysis can run without rebuilding everything from scratch.
- If those caches are removed, the workflow may need network access and will take longer.
- The committed `renv.lock` should be used for exact package version pinning once the environment build is complete.
- The analysis is intended to be run from the repository root because the code uses relative paths.

## License

This repository is licensed under GNU GPL-3. See `LICENSE` and `DESCRIPTION` for the current license declaration.
