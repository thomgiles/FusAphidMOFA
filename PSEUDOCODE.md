# Pseudocode for the FusAphidMOFA Workflow

This document explains, in plain workflow terms, what the analysis does from start to finish.

It is meant for the manuscript, supplement, and deposit review. It is not executable code.

## What the workflow is doing

The analysis takes one workbook containing several data types from the aphid and honeydew study. It then:

- reads the metadata that describes the workbook
- loads each active data sheet
- cleans the feature names and sample layout
- adds biological annotation from KEGG, lipid nomenclature, and FELLA networks
- turns the data into a format that MOFA can use
- fits the MOFA model or loads a previously trained one
- adds statistical summaries and pathway labels
- writes a large set of figures, tables, and model files into `Results/`

The core methods are based on MOFA/MOFA+, limma, DESeq2, clusterProfiler, KEGG, pathview, FELLA, ggkegg, and rgoslin [1-11].

## Overall execution

```text
START

Set the input workbook path to:
  `Input/Aphids and honeydew from exposure to niv MOFA.xlsx`

Choose the normalisation branch:
  legacy mode or updated mode

If a fresh rerun is needed:
  delete the current `Results/` directory

Load the workflow from `MOFA-CORE.R`

Prepare the R environment:
  - set the CRAN mirror if needed
  - restore the project library when `renv.lock` is available
  - install any missing packages
  - load the required packages into the session

Continue to workbook import, model fitting, and output generation
END
```

## 1. Read the workbook structure

```text
Read the `view_metadata` sheet

Use `view_metadata` to find:
  - which sheets are active
  - what type each sheet is
  - which KEGG species code to use

Read the `sample_metadata` sheet name from `view_metadata`
Read the `sample_metadata` sheet itself

Convert the first column into sample names
Transpose the table so each row is one sample
```

What this means:

- `view_metadata` is the workbook map.
- `sample_metadata` gives the sample labels and group information used later in plots and comparisons.

## 2. Load the active data sheets

```text
Find every workbook sheet marked as active

For each active sheet:
  - read the sheet into a table
  - keep the first column as the feature name column
  - make sure the sample columns line up across all views
  - clean the values where needed
  - store the result in a named list
```

What this means:

- Each omics layer is treated as one view.
- The workflow makes sure the same samples are aligned across all views before modelling starts.

## 3. Build the annotation resources

```text
For lipid data:
  - reuse stored lipid annotations when they exist
  - otherwise create them from the lipid names
  - standardise lipid naming with rgoslin / Goslin [8]

For KEGG:
  - build the gene-to-pathway table
  - build the pathway-to-name table

For transcriptomic data:
  - map transcript identifiers to KEGG-compatible gene identifiers

For metabolomic data:
  - clean metabolite names
  - match them to KEGG compounds

For FELLA:
  - load cached pathway graphs when they exist
  - otherwise build the graph and diffusion resources
  - save the result for later reuse
```

What this means:

- The data are being translated into pathway language so later plots can say not just "this feature changed", but "this feature points to this pathway".
- KEGG-based lookups and pathway maps use KEGG-related resources and tools [3, 5, 6, 9].

## 4. Normalise the data

```text
IF legacy normalisation is selected:
  use the historical workflow from the publication run
ELSE:
  use the updated normalisation routine
  include sample metadata where needed

Return a list of matrices ready for MOFA
```

What this means:

- This step turns different data types into a form MOFA can compare fairly.
- The goal is to keep the same sample order and make the data suitable for joint modelling.

## 5. Fit or load the MOFA model

```text
Create the MOFA input object from the normalised matrices

IF a trained HDF5 model already exists:
  load it
ELSE:
  fit a new model
  save it to `Results/Full_MOFA_model.hdf5`
```

What this means:

- MOFA learns hidden factors that explain the main patterns shared across the data types and the patterns unique to each one [1, 2].
- The fitted model is the centre of the analysis.

## 6. Prepare metadata and visual styling

```text
Create sample-level metadata from the sample names

Attach that metadata to the MOFA object

Create colour schemes for:
  - experimental groups
  - data views

Create a cleaned version of the MOFA object for plotting
```

What this means:

- These steps do not change the model itself.
- They make the output figures easier to read and interpret.

## 7. Add statistics and pathway labels

```text
Map KEGG labels onto the factor loadings

Run limma-style statistical testing on the factor weights [4]

Store the statistical results for later plots and tables
```

What this means:

- limma is used to help describe which features are associated with the factors [4].
- The factor weights are now linked to interpretable biological labels.

## 8. Produce the summary figures

```text
Plot how much variation each factor explains

Plot the factor matrix for samples versus factors

Plot comparisons of factor values across sample groups

Plot factor-by-factor and factor-by-covariate relationships

Write the main weight tables to disk

Plot summary scatter plots relating statistical significance to factor weight

Plot summary violin plots relating statistical significance to factor weight

Plot a faceted heatmap of the strongest weights by view
```

What this means:

- These are the broad summary figures for the manuscript.
- They show what the hidden factors look like across the whole study.

## 9. Work through each factor one by one

```text
Find the list of latent factors in the trained MOFA object

FOR each factor:
  - show the strongest features from each view
  - show the strongest features in a ranked bar chart
  - show scatter plots for the factor-specific features
  - show a heatmap of the strongest features
  - write factor-specific weight tables
  - run KEGG enrichment for that factor
  - draw an enrichment network
  - generate pathway diagrams for the most informative signals
  - run FELLA enrichment and save the network outputs
```

What this means:

- This is the "zoom in" stage.
- Each factor is examined on its own so the model can be interpreted biologically rather than just mathematically.

The enrichment and pathway visualisation stages use KEGG, pathview, clusterProfiler, FELLA, and ggkegg-style network graphics [3, 5-7, 9].

## 10. Write the output tree

```text
Save the main model file
Save the global summary plots
Save the PCA and comparison plots
Save the significance-testing outputs
Save the per-factor tables, heatmaps, scatter plots, pathway plots, and FELLA outputs
```

What this means:

- The result tree is designed to be archivable and easy to inspect.
- The file inventory is documented in `MANIFEST.md`.

## 11. Behaviour of the runner script

`MOFA-Aphid.R` is the top-level runner.

```text
Set the workbook path
Choose the normalisation mode
Delete the old `Results/` directory
Source `MOFA-CORE.R`
Zip the repository snapshot after the run
```

What this means:

- The runner is a thin wrapper.
- The real analysis logic lives in `MOFA-CORE.R`.

## 12. Reproducibility notes

```text
Run everything from the repository root

Keep the cached `Annotations/` directory when possible

Use the committed `renv.lock` to reproduce the package environment

Do not change the workbook layout unless the import logic is updated at the same time

Keep `MANIFEST.md` aligned with the actual `Results/` tree
```

What this means:

- The analysis depends on relative paths.
- The cached annotations make reruns faster and more stable.
- The lockfile is what pins the package versions.

## High-level flow

```text
Workbook
  -> metadata read
  -> active sheets loaded
  -> feature names cleaned
  -> annotation resources built
  -> data normalised
  -> MOFA fitted
  -> statistics and pathway labels added
  -> summary figures generated
  -> factor-specific plots and enrichment generated
  -> output tree written to disk
```

## References

[1] Argelaguet, R., Velten, B., Arnol, D., et al. Multi-Omics Factor Analysis - a framework for unsupervised integration of multi-omics data sets. *Molecular Systems Biology* 14(6):e8124 (2018). DOI: 10.15252/msb.20178124

[2] Argelaguet, R., Cuomo, A. S. E., Stegle, O., et al. MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data. *Genome Biology* 21:111 (2020). DOI: 10.1186/s13059-020-02015-1

[3] Kanehisa, M. and Goto, S. KEGG: Kyoto Encyclopedia of Genes and Genomes. *Nucleic Acids Research* 28(1):27-30 (2000). DOI: 10.1093/nar/28.1.27

[4] Ritchie, M. E., Phipson, B., Wu, D., et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research* 43(7):e47 (2015). DOI: 10.1093/nar/gkv007

[5] Luo, W. and Brouwer, C. Pathview: an R/Bioconductor package for pathway-based data integration and visualization. *Bioinformatics* 29(14):1830-1831 (2013). DOI: 10.1093/bioinformatics/btt285

[6] Yu, G., Wang, L.-G., Han, Y., and He, Q.-Y. clusterProfiler: an R package for comparing biological themes among gene clusters. *OMICS* 16(5):284-287 (2012). DOI: 10.1089/omi.2011.0118

[7] Castellano-Escuder, P., Andrés-León, E., Rojas, A., and Gajardo-Viñas, S. FELLA: an R package to enrich metabolomics data. *BMC Bioinformatics* 20:22 (2019). DOI: 10.1186/s12859-018-2566-9

[8] Kopczynski, D., Hoffmann, N., Peng, B., and Ahrends, R. Goslin: A Grammar of Succinct Lipid Nomenclature. *Analytical Chemistry* 92(16):10957-10960 (2020). DOI: 10.1021/acs.analchem.0c01690

[9] Maeda, N., et al. ggkegg: analysis and visualization of KEGG data utilizing the grammar of graphics. *Bioinformatics* 39(10):btad622 (2023). DOI: 10.1093/bioinformatics/btad622

[10] Love, M. I., Huber, W., and Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15:550 (2014). DOI: 10.1186/s13059-014-0550-8

[11] Gu, Z., Eils, R., and Schlesner, M. Complex heatmaps reveal patterns and correlations in multidimensional genomic data. *Bioinformatics* 32(18):2847-2849 (2016). DOI: 10.1093/bioinformatics/btw313

[12] Package documentation and citations for `KEGGREST`, `readxl`, and `openxlsx` were used as implementation references where no single canonical journal article was the main citation source.
