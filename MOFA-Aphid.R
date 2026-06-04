## MOFA script for Aphids ##

input_file="Input/Aphids and honeydew from exposure to niv MOFA.xlsx"

OLD_NORMALISATION=TRUE
unlink(recursive = T, "./Results")
source("./MOFA-CORE.R");

###########################################################
###########################################################
# CODE EXECUTION
###########################################################
###########################################################

# ==== Initialise the environment ====

# Sets up project environment, optionally with renv support
setup_project_environment(use_renv = TRUE)


# ==== Load and pre-process sample metadata ====

# Load raw Excel sheets, clean up 'Name' columns, drop unused columns,
# harmonize column ordering across views

view_metadata <- import_view_metadata(input_file)
sample_metadata <- import_sample_metadata(input_file, view_metadata)
data_list <- load_and_clean_input_data(input_file, view_metadata)

# Load existing or generate new lipid annotations from lipidomic sheet names
LIPID_ANNOTATIONS <- generate_lipid_annotations(view_metadata, data_list)

# Downloads or loads TERM2GENE and TERM2NAME for gene and metabolite KEGG pathways
KEGG_list <- generate_kegg_annotation_lists(view_metadata)
TERM2GENE_list <- KEGG_list$TERM2GENE
TERM2NAME_list <- KEGG_list$TERM2NAME

# Maps transcriptomic features (e.g. Ensembl IDs) to KEGG gene IDs via Ensembl → Entrez conversion
gene_annotation <- annotate_transcriptomic_data_with_kegg(data_list, KEGG_list$TERM2GENE, view_metadata)

# Matches compound names (post-pipe cleaned) to KEGG compound IDs using keggFind()
metabolite_annotation <- annotate_metabolomic_data_with_kegg(data_list,  KEGG_list$TERM2GENE, view_metadata)

# Builds or loads FELLA graph and diffusion matrix for pathway enrichment
fella_list <- load_or_build_fella_database(view_metadata)

# ==== Data normalisation and MOFA model training ====

# Normalise (impute, and re-structure data for MOFA)"

if(OLD_NORMALISATION){
  print("using v1 normalization method")
  mofa_list <- normalise_data_mofa_v1(data_list)
} else {
  mofa_list <- normalise_data_mofa_v2(data_list,sample_metadata)
}

# Create and train a new MOFA model, or load existing hdf5
MOFAobject.trained <- initialise_and_run_mofa(mofa_list)

# ==== Sample metadata and colour setup for visualisation ====

# Create metadata from sample names (e.g. categories), binary indicators, groupings
metadata <- generate_sample_metadata(colnames(mofa_list[[1]]))

# Attach sample metadata to MOFA object
MOFAobject.trained <- assign_sample_metadata_to_mofa(MOFAobject.trained, metadata)

# Define colour palettes per group and per view
group.colors <- define_group_colours(metadata)

view.colors <- define_view_colours(names(MOFAobject.trained@data))

# ==== Post-process MOFA model  ====

# Add KEGG annotations
MOFAobject.trained <- extract_kegg_weights(MOFAobject.trained,
                                           view_metadata,
                                           metabolite_annotation,
                                           gene_annotation)
# Add Stats annotations
MOFAobject.trained <-  run_limma_statistics(MOFAobject.trained, group.colors)

# Create a MOFA object with standardised feature names (for plotting)
MOFAobject.trained.cleanNames <- standardise_feature_names(MOFAobject.trained)

# ==== Visualisation of MOFA global outputs ====

# Plot proportion of variance explained
generate_mofa_variance_plots(MOFAobject.trained, view.colors)

# Plot factor matrices (samples x factors)
generate_mofa_factor_matrix_plot(MOFAobject.trained, view.colors, group.colors)

# Plot comparisons of factor values across sample groups
generate_mofa_factor_comparison_plot(MOFAobject.trained, group.colors, view.colors, metadata)

# Plot factor–factor and factor–covariate correlation matrices
generate_mofa_factor_correlation_plot(MOFAobject.trained, view.colors, sample_metadata)

# Save weight tables
save_weight_tables(MOFAobject.trained, group.colors)

# Generate scatter plots comaring DE to factor weights
plot_scatter_pval_vs_weight(view.colors)

# Generate Violin plots comaring DE to factor weights
plot_violin_pval_vs_weight(view.colors)

plot_top_weights_heatmap_faceted_by_view(MOFAobject.trained)

# ==== Per-factor MOFA output generation ====

# Extract list of latent factor names from MOFA
factors <- colnames(get_factors(MOFAobject.trained)[["group1"]])

# Plot top-ranked features per view
for (factor_name in factors) {
  message(paste("plotting top weights for", factor_name))
  plot_top_weights_view(factor_name)
}

# Plot top-ranked barchart
for (factor_name in factors) {
  message(paste("plotting top weights for", factor_name))
  plot_top_weights_per_factor(factor_name, top_n = 10)
}

# Plot scatter plots of top feature weights
for (factor_name in factors) {
  message(paste("plotting scatter for", factor_name))
  plot_scatter_per_view(factor_name)
}

# Plot heatmap of top-ranked features
for (factor_name in factors) {
  message(paste("plotting heatmap for", factor_name))
  plot_heatmap_per_view(factor_name)
}

# Save all weights in tidy format
for (factor_name in factors) {
  message(paste("creating weight outputs for", factor_name))
  create_weight_outputs(factor_name, view_metadata)
}

# Perform KEGG enrichment
for (factor_name in factors) {
  message(paste("performing KEGG enrichment for", factor_name))
  perform_kegg_enrichment(factor_name, MOFAobject.trained, view_metadata, TERM2GENE_list, TERM2NAME_list)
}

# Plot enrichment network -> needs work
for (factor_name in factors) {
  message(paste("plotting enrichment network for", factor_name))
  plot_enrichment_network(factor_name, gene_annotation, metabolite_annotation, view.colors)
}


for (factor_name in factors) {
  message(paste("generating singular pathway visualisations for", factor_name))
  generate_simple_pathway_graph(factor_name, view_metadata, MOFAobject.trained)
}


# Perform FELLA enrichment
for (factor_name in factors) {
  message(paste("running FELLA enrichment for", factor_name))
  run_fella_enrichment(factor_name, view_metadata, fella_list, TERM2GENE_list, MOFAobject.trained)
}

###########################################################
###########################################################
# Tidy up files and create zip for sharing
###########################################################
###########################################################

outfile <- paste0("../",basename(dir),"-",format(Sys.Date(), "%d-%B-%Y"),".zip")

utils::zip(zipfile = outfile, files = ".", flags = "-r9Xq")

# Remove the copied file
unlink("./MOFA-CORE.R")