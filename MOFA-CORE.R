###########################################################
###########################################################
# FUNCTIONS
###########################################################
###########################################################

setup_project_environment <- function(use_renv = TRUE,
                                      renv_path = "renv") {
  print("Initialising project environment...")
  
  if (isTRUE(getOption("repos")[["CRAN"]]) || identical(getOption("repos")[["CRAN"]], "@CRAN@")) {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
  }
  
  if (use_renv) {
    if (!requireNamespace("renv", quietly = TRUE)) {
      install.packages("renv")
      print("Installed 'renv' for virtual environment support.")
    } else {
      print("'renv' already installed.")
    }
    
    renv_lock <- if (file.exists(file.path(renv_path, "renv.lock"))) file.path(renv_path, "renv.lock") else "renv.lock"
    if (file.exists(renv_lock)) {
      print("Restoring existing renv environment...")
      renv::restore(prompt = FALSE)
    } else {
      print("No renv.lock found; continuing without restore.")
    }
  }
  
  # Ensure BiocManager is available
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
    print("Installed 'BiocManager'.")
  }
  
  # CRAN packages
  cran_packages <- c(
    "data.table",
    "ggplot2",
    "mvtnorm",
    "psych",
    "dplyr",
    "ggpubr",
    "GGally",
    "readxl",
    "readr",
    "tidyr",
    "stringi",
    "openxlsx",
    "reshape2",
    "writexl",
    "textclean",
    "stringr",
    "viridis",
    "grid",
    "gtable",
    "gridExtra",
    "gridGraphics",
    "stats",
    "magick",
    "tibble",
    "tidygraph",
    "ggraph",
    "igraph",
    "scales",
    "ggtext",
    "purrr",
    "XML",
    "png",
    "httr",
    "httr2",
    "jsonlite",
    "VIM",
    "preprocessCore",
    "matrixStats",
    "ggforce"
  )
  
  # Bioconductor packages
  bioc_packages <- c(
    "tximport",
    "DESeq2",
    "MOFA2",
    "KEGGREST",
    "MSnbase",
    "clusterProfiler",
    "pathview",
    "preprocessCore",
    "ggkegg",
    "rgoslin",
    "fgsea",
    "sva",
    "limma",
    "FELLA",
    "biomaRt",
    "ComplexHeatmap"
  )
  
  # Install missing CRAN packages
  install_cran <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      print(paste("Installing CRAN package:", pkg))
      install.packages(pkg, dependencies = TRUE)
    } else {
      print(paste("CRAN package already installed:", pkg))
    }
  }
  
  # Install missing Bioconductor packages
  install_bioc <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      print(paste("Installing Bioconductor package:", pkg))
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      print(paste("Bioconductor package already installed:", pkg))
    }
  }
  
  print("Checking and installing CRAN packages...")
  sapply(unique(cran_packages), install_cran)
  
  print("Checking and installing Bioconductor packages...")
  sapply(bioc_packages, install_bioc)
  
  # Load all packages
  print("Loading all required packages into session...")
  invisible(lapply(c(cran_packages, bioc_packages), library, character.only = TRUE))
  
  print("All packages successfully installed and loaded.")
}

######### GET VIEW METADATA #########
import_view_metadata <- function(input_file){
  print("Reading metadata sheet...")
  
  view_metadata <- as.data.frame(readxl::read_excel(input_file, sheet = "view_metadata"))
  return(view_metadata)
}

######### GET SAMPLE METADATA #########
import_sample_metadata <- function(input_file, view_metadata){
  print("Reading metadata sheet...")
  
  metadata_sheet <- view_metadata$View[view_metadata$Type =="sample_metadata"]
  if(!is.null(metadata_sheet)){
    metadata_list <- as.data.frame(readxl::read_excel(input_file, sheet = metadata_sheet))
    
    # Set first column ("Name") as rownames
    rownames(metadata_list) <- metadata_list$Name
    metadata_list <- metadata_list[ , -1, drop = FALSE]  # drop the Name column
    metadata_list_t <- as.data.frame(t(metadata_list))
    
    return(metadata_list_t)
  } 
}


###################### Load and Process Data #######################
load_and_clean_input_data <- function(input_file, view_metadata) {
  print("Loading and processing input data...")
  
  sheet_names <- subset(view_metadata, Active == "Y", select = View)$View
  
  print(sheet_names)
  
  # Read Excel sheets and create named list
  print("Reading input sheets...")
  data_list <- setNames(lapply(sheet_names, function(sheet) {
    as.data.frame(readxl::read_excel(input_file, sheet = sheet))
  }), sheet_names)
  
  print(head(data_list))
  
  # Define Greek letter replacement function
  replace_greek_letters <- function(text) {
    greek_map <- c(
      "α" = "alpha",
      "β" = "beta",
      "γ" = "gamma",
      "δ" = "delta",
      "ε" = "epsilon",
      "ζ" = "zeta",
      "η" = "eta",
      "θ" = "theta",
      "ι" = "iota",
      "κ" = "kappa",
      "λ" = "lambda",
      "μ" = "mu",
      "ν" = "nu",
      "ξ" = "xi",
      "ο" = "omicron",
      "π" = "pi",
      "ρ" = "rho",
      "σ" = "sigma",
      "τ" = "tau",
      "υ" = "upsilon",
      "φ" = "phi",
      "χ" = "chi",
      "ψ" = "psi",
      "ω" = "omega",
      "ς" = "sigma",
      "%" = "",
      "," = ""
    )
    for (char in names(greek_map)) {
      text <- gsub(char, greek_map[[char]], text, fixed = TRUE)
    }
    return(text)
  }
  
  # Clean 'Name' column
  clean_name_column <- function(df) {
    df$Name <- sapply(df$Name, function(x) {
      x <- replace_greek_letters(x)
      x <- iconv(x,
                 from = "UTF-8",
                 to = "ASCII//TRANSLIT",
                 sub = "")
      return(x)
    })
    return(df)
  }
  data_list <- lapply(data_list, clean_name_column)
  
  # Format dataframes: set row names and drop first column
  adjust_dataframe <- function(df) {
    df <- df %>% `row.names<-`(df$Name) %>% dplyr::select(-1)
    return(df)
  }
  data_list <- lapply(data_list, adjust_dataframe)
  
  # Get the reference column order from the first view
  ref_order <- colnames(data_list[[1]])
  
  # Reorder all views to match the reference
  data_list <- lapply(data_list, function(df)
    df[, ref_order])
  
  print("Data loading and initial cleaning complete.")
  return(data_list)
}






######### PROCESS LIPID annotations ###########

generate_lipid_annotations <- function(view_metadata,
                                       data_list,
                                       lipid_annotation_file = "Annotations/Lipid.annotations.csv") {
  
  
  lipidomic_sheets <- subset(view_metadata, Active == "Y" & Type == "Lipid", select = View)$View
  
  if (length(lipidomic_sheets) == 0) {
    message("No lipidomic sheets provided. Skipping lipid annotation.")
    return(NULL)
  }
  
  print("Processing lipid annotations...")
  
  if (file.exists(lipid_annotation_file)) {
    message("Loading existing lipid annotations from file...")
    
    if (file.info(lipid_annotation_file)$size < 5) {
      print("Existing lipid annotation file is empty. Proceeding without annotations.")
      return(data.frame())
    } else {
      print("Lipid annotations loaded.")
      return(read.csv(lipid_annotation_file))
    }
  }
  
  print("No annotation file found. Building lipid annotations from scratch...")
  
  lipid_annotations <- data.frame(
    name = character(),
    CleanNames = character(),
    ShortAnnotation = character(),
    Description = character(),
    stringsAsFactors = FALSE
  )
  
  for (lipid_sheet in lipidomic_sheets) {
    print(paste("Processing sheet:", lipid_sheet))
    lipid_names <- rownames(data_list[[lipid_sheet]])
    
    for (lipid_name in lipid_names) {
      clean_lipid_name <- gsub("'", "", lipid_name)
      clean_lipid_names <- unlist(strsplit(clean_lipid_name, "/"))
      GoslinDF <- parseLipidNames(clean_lipid_names)
      GoslinDF$RefMet_ID <- NULL
      
      required_cols <- c(
        "Normalized.Name",
        "Message",
        "Grammar",
        "Mass",
        "Sum.Formula",
        "Extended.Species.Name",
        "Lipid.Maps.Main.Class",
        "Lipid.Maps.Category",
        "Functional.Class.Abbr",
        "RefMet_ID"
      )
      
      for (col in required_cols) {
        if (!(col %in% colnames(GoslinDF))) {
          GoslinDF[[col]] <- NA
        }
      }
      
      for (i in 1:nrow(GoslinDF)) {
        original.Name <- GoslinDF[i, "Original.Name"]
        query_url <- paste0(
          "https://www.metabolomicsworkbench.org/rest/refmet/name/",
          URLencode(original.Name, reserved = TRUE),
          "/all"
        )
        response <- httr::GET(query_url)
        
        if (httr::status_code(response) == 200) {
          json_data <- httr::content(response, as = "text", encoding = "UTF-8")
          parsed_data <- jsonlite::fromJSON(json_data, flatten = TRUE)
          
          if (length(parsed_data) > 1) {
            GoslinDF[i, "Normalized.Name"] <- parsed_data$name
            GoslinDF[i, "Message"] <- "NA"
            GoslinDF[i, "Grammar"] <- "refmet"
            GoslinDF[i, "Mass"] <- parsed_data$exactmass
            GoslinDF[i, "Sum.Formula"] <- parsed_data$formula
            GoslinDF[i, "Extended.Species.Name"] <- "ST"
            GoslinDF[i, "Lipid.Maps.Main.Class"] <- parsed_data$main_class
            GoslinDF[i, "Functional.Class.Abbr"] <- parsed_data$sub_class
            GoslinDF[i, "RefMet_ID"] <- parsed_data$refmet_id
          }
        }
      }
      
      GoslinDF$Functional.Class.Abbr <- gsub("\\[|\\]", "", GoslinDF$Functional.Class.Abbr)
      
      GoslinDF <- GoslinDF %>%
        dplyr::summarise(across(everything(), ~ paste(unique(.), collapse = " | ")))
      
      new_annotation <- data.frame(
        name = lipid_name,
        CleanNames = GoslinDF[, "Normalized.Name"],
        ShortAnnotation = GoslinDF[, "Lipid.Maps.Main.Class"],
        Description = paste0(
          "MainClass=[",
          GoslinDF[, "Lipid.Maps.Main.Class"],
          "],ExtendedClass=[",
          GoslinDF[, "Extended.Species.Name"],
          "],FunctionalClass=[",
          GoslinDF[, "Functional.Class.Abbr"],
          "],Mass=[",
          GoslinDF[, "Mass"],
          "],Formula=[",
          GoslinDF[, "Sum.Formula"],
          "]"
        ),
        stringsAsFactors = FALSE
      )
      
      lipid_annotations <- rbind(lipid_annotations, new_annotation)
    }
  }
  
  print("Finalising lipid annotations...")
  lipid_annotations <- unique(lipid_annotations)
  write.csv(lipid_annotations, file = lipid_annotation_file, row.names = FALSE)
  print("Lipid annotations written to file.")
  
  return(lipid_annotations)
}


######################  Get the KEGG pathway data #######################

# ---- KEGG annotation generator (species → {gene, metabolite}) ----
generate_kegg_annotation_lists <- function(view_metadata,
                                           cache_dir = "Annotations/kegg_annotation_cache") {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  
  kegg_species <- unique(unlist(strsplit(na.omit(subset(view_metadata, Active == "Y", select = Kegg_Species)$Kegg_Species), ",")))
  
  TERM2GENE_list <- list()
  TERM2NAME_list <- list()
  
  for (sp in kegg_species) {
    message("Processing KEGG species: ", sp)
    pathways_file <- file.path(cache_dir, paste0(sp, "_pathways.rds"))
    
    # ---- Retrieve or load pathways ----
    if (file.exists(pathways_file)) {
      species_pathways <- readRDS(pathways_file)
      message("Loaded ", length(species_pathways), " pathways for ", sp)
    } else {
      species_pathways <- KEGGREST::keggList("pathway", sp)
      saveRDS(species_pathways, pathways_file)
      message("Downloaded ", length(species_pathways), " pathways for ", sp)
    }
    
    species_ids <- sub("path:", "", names(species_pathways))
    map_ids <- sub(paste0("^", sp), "map", species_ids)
    
    TERM2NAME <- data.frame(
      ID = species_ids,
      Name = as.character(species_pathways),
      stringsAsFactors = FALSE
    )
    
    # Replace the organism name in TERM2NAME$Name
    TERM2NAME$Name <- sub(" - .*", paste0(" - ", sp), TERM2NAME$Name)
    
    translator <- data.frame(
      species = species_ids,
      map = map_ids,
      stringsAsFactors = FALSE
    )
    
    # ---- Metabolite mappings ----
    compound2map <- KEGGREST::keggLink("pathway", "compound")
    compound_ids <- sub("^cpd:", "", names(compound2map))
    compound_pathways <- sub("^path:", "", compound2map)
    
    compound_links <- data.frame(
      map = compound_pathways,
      KEGG_ID = compound_ids,
      stringsAsFactors = FALSE
    )
    
    compound_links_merged <- merge(translator, compound_links, by = "map")
    TERM2GENE_metabolite <- compound_links_merged[, c("species", "KEGG_ID")]
    colnames(TERM2GENE_metabolite) <- c("Pathway", "KEGG_ID")
    
    # ---- Gene mappings ----
    gene2pathway_raw <- KEGGREST::keggLink("pathway", sp)
    gene_ids <- sub(paste0("^", sp, ":"), "", names(gene2pathway_raw))
    pathway_ids <- sub("^path:", "", gene2pathway_raw)
    
    TERM2GENE_gene <- data.frame(
      Pathway = pathway_ids,
      KEGG_ID = gene_ids,
      stringsAsFactors = FALSE
    )
    
    # ---- Store for this species ----
    TERM2GENE_list[[sp]] <- list(
      gene = TERM2GENE_gene,
      metabolite = TERM2GENE_metabolite
    )
    
    TERM2NAME_list[[sp]] <- list(
      gene = TERM2NAME,
      metabolite = TERM2NAME
    )
  }
  
  return(list(TERM2GENE = TERM2GENE_list, TERM2NAME = TERM2NAME_list))
}


######################  Annotate Transcriptomics sheets #######################
# ---- Annotate transcriptomic data with KEGG gene IDs ----
annotate_transcriptomic_data_with_kegg <- function(
    data_list,
    TERM2GENE_list,
    view_metadata,
    outdir = "Annotations"
) {
  # Derive transcriptomic sheets from view_metadata
  transcriptomic_sheets <- subset(view_metadata, Active == "Y" & Type == "RNAseq", select = View)$View
  if (length(transcriptomic_sheets) == 0) {
    warning("No active transcriptomic views found in view_metadata.")
    return(NULL)
  }
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  gene_annotation <- list()
  
  # Collect all unique genes across transcriptomic sheets
  gene_list <- unlist(lapply(transcriptomic_sheets, function(sheet) {
    genes <- rownames(data_list[[sheet]])
    gsub(paste0("_", sheet), "", genes)
  }))
  gene_list <- unique(gene_list)
  message("Total unique gene identifiers to annotate: ", length(gene_list))
  
  # Extract KEGG species associated with transcriptomic views
  transcriptomic_species <- unique(na.omit(subset(
    view_metadata,
    View %in% transcriptomic_sheets & Type == "RNAseq",
    select = Kegg_Species
  )$Kegg_Species))
  
  # Flatten comma-separated species
  transcriptomic_species <- unique(unlist(strsplit(transcriptomic_species, ",")))
  
  if (length(transcriptomic_species) == 0) {
    warning("No KEGG species codes found in view_metadata for transcriptomic views.")
    return(NULL)
  }
  
  # Loop over species
  for (sp in transcriptomic_species) {
    message("Annotating genes for species: ", sp)
    output_file <- file.path(outdir, paste0("KEGG_gene_annotations_", sp, ".csv"))
    
    # If cached, load and return
    if (file.exists(output_file)) {
      message("Loading existing gene annotations for ", sp, " from file.")
      gene_annotation[[sp]] <- utils::read.csv(output_file, stringsAsFactors = FALSE)
      next
    } else {
      # Extract KEGG TERM2GENE table for this species
      if (!("gene" %in% names(TERM2GENE_list[[sp]]))) {
        warning("No gene TERM2GENE entry for species ", sp)
        next
      }
      kegg_df <- TERM2GENE_list[[sp]]$gene
      kegg_df$EntrezID <- sub(paste0("^", sp, ":"), "", kegg_df$KEGG_ID)
      
      # Assume input genes are already Entrez IDs
      ann <- merge(
        data.frame(gene = gene_list, stringsAsFactors = FALSE),
        kegg_df,
        by.x = "gene",
        by.y = "EntrezID"
      )
      ann <- ann[, c("gene", "KEGG_ID")]
      
      message("Matched KEGG IDs for ", nrow(ann), " genes in ", sp)
      utils::write.csv(ann, output_file, row.names = FALSE)
      gene_annotation[[sp]] <- ann
      
    }
  }
  

  
  return(gene_annotation)
}

######################  Annotate Metabolite sheets #######################

# ---- Annotate metabolomic data with KEGG compound IDs ----
annotate_metabolomic_data_with_kegg <- function(
    data_list,
    TERM2GENE_list,
    view_metadata,
    outdir = "Annotations"
) {
  # Derive metabolomic sheets from view_metadata
  metabolomic_views <- subset(view_metadata, Active == "Y" & Type == "Metabolic", select = View)$View
  if (length(metabolomic_views) == 0) {
    warning("No active metabolomic views found in view_metadata.")
    return(NULL)
  }
  
  # Derive species list from metabolomic views
  metabolomic_species <- unique(na.omit(subset(
    view_metadata,
    View %in% metabolomic_views & Type == "Metabolic",
    select = Kegg_Species
  )$Kegg_Species))
  metabolomic_species <- unique(unlist(strsplit(metabolomic_species, ",")))
  
  if (length(metabolomic_species) == 0) {
    warning("No KEGG species codes found in view_metadata for metabolomic views.")
    return(NULL)
  }
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  metabolite_annotation <- list()
  
  # Step 1: Clean metabolite names across sheets
  compound_list <- unlist(lapply(metabolomic_views, function(sheet) {
    compounds <- rownames(data_list[[sheet]])
    compounds <- gsub(paste0("_", sheet), "", compounds)
    compounds <- trimws(sub(".*\\|", "", compounds))  # Keep part after pipe
    return(compounds)
  }))
  compound_list <- unique(compound_list)
  message("Total unique cleaned compound names: ", length(compound_list))
  
  # Helper: query KEGG for a single metabolite
  get_kegg_id <- function(metabolite_name, max_retries = 2) {
    attempt <- 1
    while (attempt <= max_retries) {
      result <- tryCatch({
        KEGGREST::keggFind("compound", metabolite_name)
      }, error = function(e) {
        return(NULL)
      })
      if (!is.null(result) && length(result) > 0) {
        return(names(result)[1])  # First match
      }
      attempt <- attempt + 1
    }
    return(NA)
  }
  
  # Step 2: For each species, generate annotation
  for (sp in metabolomic_species) {
    message("Annotating metabolites for species: ", sp)
    output_file <- file.path(outdir, paste0("KEGG_metabolite_annotations_", sp, ".csv"))
    
    # If cached, load
    if (file.exists(output_file)) {
      message("Loading existing metabolite annotations for ", sp)
      metabolite_annotation[[sp]] <- utils::read.csv(output_file, stringsAsFactors = FALSE)
      next
    }
    
    # Add progress bar
    n <- length(compound_list)
    pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
    
    kegg_ids <- vector("character", n)
    for (i in seq_along(compound_list)) {
      kegg_ids[i] <- get_kegg_id(compound_list[i])
      utils::setTxtProgressBar(pb, i)
    }
    close(pb)
    
    # Combine results
    ann <- data.frame(
      Metabolite = compound_list,
      KEGG_ID = kegg_ids,
      stringsAsFactors = FALSE
    )
    
    # Filter out compounds without a KEGG match
    before_filter <- nrow(ann)
    ann <- ann[!is.na(ann$KEGG_ID), ]
    after_filter <- nrow(ann)
    
    message("Annotated ", after_filter, " of ", before_filter, " metabolites for ", sp)
    
    # Clean KEGG ID formatting
    ann$KEGG_ID <- sub("^cpd:", "", ann$KEGG_ID)
    
    # Save and return per species
    utils::write.csv(ann, output_file, row.names = FALSE)
    metabolite_annotation[[sp]] <- ann
  }
  
  return(metabolite_annotation)
}



##############################
# Get KEGG graph data for FELLA analysis (multi-species support)
##############################
load_or_build_fella_database <- function(view_metadata, fella_dir = "Annotations/fella_db") {
  
  
  # Ensure base directory exists
  if (!dir.exists(fella_dir)) {
    dir.create(fella_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Derive metabolomic sheets from view_metadata
  metabolomic_views <- subset(view_metadata, Active == "Y" & Type == "Metabolic", select = View)$View
  if (length(metabolomic_views) == 0) {
    warning("No active metabolomic views found in view_metadata.")
    return(NULL)
  }
  
  metabolomic_species <- unique(na.omit(subset(
    view_metadata,
    View %in% metabolomic_views & Type == "Metabolic",
    select = Kegg_Species
  )$Kegg_Species))
  metabolomic_species <- unique(unlist(strsplit(metabolomic_species, ",")))
  
  # Container for per-species FELLA objects
  fella_list <- list()
  
  for (sp in metabolomic_species) {
    message("Processing FELLA database for species: ", sp)
    
    # Define per-species paths
    fella_db_dir <- file.path(fella_dir, sp)
    graph_file   <- file.path(fella_db_dir, paste0(sp,"_kegg_graph.rds"))
    
    # Ensure species directory exists
    if (!dir.exists(fella_db_dir)) {
      dir.create(fella_db_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    # Load or create KEGG graph
    if (file.exists(graph_file)) {
      message("  Loading KEGG graph for ", sp, " from disk...")
      kegg_graph <- readRDS(graph_file)
    } else {
      message("  Building KEGG graph for ", sp, " from KEGG REST...")
      kegg_graph <- FELLA::buildGraphFromKEGGREST(organism = sp)
      saveRDS(kegg_graph, graph_file)
    }
    
    # Load or build FELLA database
    if (dir.exists(fella_db_dir) &&
        file.exists(file.path(paste0(fella_db_dir,"/graph"), "diffusion.matrix.RData"))) {
      message("  Loading FELLA database for ", sp, " from disk...")
      fella.data <- FELLA::loadKEGGdata(
        databaseDir = paste0(fella_db_dir,"/graph"),
        internalDir = FALSE,
        loadMatrix  = "diffusion"
      )
    } else {
      message("  Building FELLA database for ", sp, "...")
      FELLA::buildDataFromGraph(
        keggdata.graph = kegg_graph,
        databaseDir    = paste0(fella_db_dir,"/graph"),
        internalDir    = FALSE,
        matrices       = c("diffusion"),
        normality      = c("diffusion"),
        niter          = 100
      )
      
      fella.data <- FELLA::loadKEGGdata(
        databaseDir = paste0(fella_db_dir,"/graph"),
        internalDir = FALSE,
        loadMatrix  = "diffusion"
      )
    }
    
    # Store result under species key
    fella_list[[sp]] <- fella.data
  }
  
  return(fella_list)
}



##############################
# Normalise data
##############################

normalise_data_mofa_v1 <- function(data_list,
                                na_limit = 1e-7) {
  stopifnot(is.list(data_list))
  message("Normalising data...")
  
  result <- vector("list", length(data_list))
  names(result) <- names(data_list)
  
  for (view_name in names(data_list)) {
    mat <- as.matrix(data_list[[view_name]])
    mode(mat) <- "numeric"
    
    treat_zeros_as_NA <- startsWith(view_name, "RNASeq") ||
      startsWith(view_name, "LCMS")
    
    ## modality-specific zero handling
    if (treat_zeros_as_NA) {
      mat[abs(mat) <= na_limit] <- NA
    } else {
      mat[mat < 0] <- NA
      if (any(mat > 0, na.rm = TRUE)) {
        min_nonzero <- min(mat[mat > 0], na.rm = TRUE)
        mat <- mat + 0.1 * min_nonzero
      }
    }
    
    ## log2 transform
    mat <- suppressWarnings(log2(mat))
    mat[!is.finite(mat)] <- NA
    
    ## row-wise z-scoring
    mat <- t(apply(mat, 1, function(x) {
      if (sum(!is.na(x)) > 1) as.numeric(scale(x)) else rep(NA_real_, length(x))
    }))
    rownames(mat) <- paste0(rownames(data_list[[view_name]]), "_", view_name)
    colnames(mat) <- colnames(data_list[[view_name]])
    
    ## remove duplicated rows by name
    dup_idx <- duplicated(rownames(mat))
    if (any(dup_idx)) {
      message(
        "Removing ", sum(dup_idx), " duplicated rows: ",
        paste(unique(rownames(mat)[dup_idx]), collapse = ", ")
      )
      mat <- mat[!dup_idx, , drop = FALSE]
    }
    
    ## remove rows with all 0/NA
    keep <- apply(mat, 1, function(x) any(!is.na(x) & x != 0))
    removed <- sum(!keep)
    if (removed > 0) {
      message("Removed ", removed, " rows with all values 0/NA.")
    }
    mat <- mat[keep, , drop = FALSE]
    
    result[[view_name]] <- mat
    
    message(sprintf(
      "Processed view '%s': %d → %d features after filtering.",
      view_name,
      nrow(data_list[[view_name]]),
      nrow(mat)
    ))
  }
  
  result
}

##############################
# Normalise data for MOFA (v3)
# - log2(x + 1)
# - supervised filtering using limma moderated F when possible
# - if limma fails / no residual df: add tiny noise + retry, else fallback to row variance
# - row-wise z-scoring using matrixStats
##############################

normalise_data_mofa_v2 <- function(data_list,
                                   sample_metadata,
                                   group_col = "category",
                                   top_n = 1000L,
                                   noise_sd = 1e-8,
                                   noise_seed = 1L,
                                   retry_with_noise = TRUE) {
  stopifnot(is.list(data_list))
  stopifnot(is.data.frame(sample_metadata))
  stopifnot(group_col %in% colnames(sample_metadata))
  
  requireNamespace("limma")
  requireNamespace("matrixStats")
  
  message("Normalising data...")
  
  result <- list()
  
  add_tiny_noise_to_constant_rows <- function(mat, noise_sd, seed) {
    # Adds noise only to rows with zero variance (ignoring NAs)
    rsd <- matrixStats::rowSds(mat, na.rm = TRUE)
    const <- is.finite(rsd) & rsd == 0
    if (!any(const)) return(list(mat = mat, n = 0L))
    
    set.seed(seed)
    # scale noise to magnitude (avoid all-zero scale)
    scale_vec <- pmax(1, abs(matrixStats::rowMeans2(mat[const, , drop = FALSE], na.rm = TRUE)))
    # recycle scale per row across columns
    noise <- matrix(rnorm(sum(const) * ncol(mat), mean = 0, sd = noise_sd),
                    nrow = sum(const), ncol = ncol(mat))
    noise <- noise * scale_vec
    mat[const, ] <- mat[const, ] + noise
    
    list(mat = mat, n = as.integer(sum(const)))
  }
  
  for (view_name in names(data_list)) {
    message("Processing view: ", view_name)
    
    mat <- data_list[[view_name]]
    mat <- as.matrix(mat)
    mode(mat) <- "numeric"
    
    # Basic cleaning
    mat[mat < 0] <- NA
    mat <- mat[rowSums(!is.na(mat) & mat != 0) > 0, , drop = FALSE]
    
    # Log2 transform
    mat <- log2(mat + 1)
    mat[!is.finite(mat)] <- NA
    
    if (nrow(mat) < 2L) {
      warning("Skipping view ", view_name, ": fewer than 2 features after filtering")
      next
    }
    
    # Check sample alignment (strongly recommended)
    if (ncol(mat) != nrow(sample_metadata)) {
      warning(
        "View ", view_name,
        ": ncol(mat) = ", ncol(mat),
        " but nrow(sample_metadata) = ", nrow(sample_metadata),
        ". Ensure columns are samples in the same order as sample_metadata."
      )
    }
    
    group <- factor(sample_metadata[[group_col]])
    design <- model.matrix(~ 0 + group)
    
    # Helper: unsupervised ranking fallback
    rank_by_variance <- function(mat) {
      v <- matrixStats::rowVars(mat, na.rm = TRUE)
      v[is.na(v)] <- -Inf
      v
    }
    
    # Try supervised limma ranking when residual df exists
    rank_scores <- NULL
    limma_ok <- TRUE
    
    # residual df = n_samples - rank(design)
    res_df <- ncol(mat) - qr(design)$rank
    if (!is.finite(res_df) || res_df <= 0) {
      limma_ok <- FALSE
    } else {
      fit <- limma::lmFit(mat, design)
      fit2 <- tryCatch(
        limma::eBayes(fit),
        error = function(e) e
      )
      
      if (inherits(fit2, "error")) {
        if (isTRUE(retry_with_noise)) {
          # Add tiny noise to constant rows then retry once
          jittered <- add_tiny_noise_to_constant_rows(mat, noise_sd = noise_sd, seed = noise_seed)
          if (jittered$n > 0L) {
            mat2 <- jittered$mat
            fit_retry <- limma::lmFit(mat2, design)
            fit_retry2 <- tryCatch(
              limma::eBayes(fit_retry),
              error = function(e) e
            )
            if (!inherits(fit_retry2, "error")) {
              mat <- mat2
              fit2 <- fit_retry2
            } else {
              limma_ok <- FALSE
            }
          } else {
            limma_ok <- FALSE
          }
        } else {
          limma_ok <- FALSE
        }
      }
      
      if (limma_ok) {
        rank_scores <- fit2$F
        if (is.null(rank_scores)) limma_ok <- FALSE
      }
    }
    
    if (!limma_ok) {
      # If rows are constant, optionally jitter them so z-scoring doesn't collapse everything
      if (isTRUE(retry_with_noise)) {
        jittered <- add_tiny_noise_to_constant_rows(mat, noise_sd = noise_sd, seed = noise_seed)
        if (jittered$n > 0L) mat <- jittered$mat
      }
      rank_scores <- rank_by_variance(mat)
    }
    
    rank_scores[is.na(rank_scores)] <- -Inf
    keep_idx <- order(rank_scores, decreasing = TRUE)[seq_len(min(top_n, length(rank_scores)))]
    mat <- mat[keep_idx, , drop = FALSE]
    
    # Row-wise z-scoring
    row_means <- matrixStats::rowMeans2(mat, na.rm = TRUE)
    row_sds   <- matrixStats::rowSds(mat, na.rm = TRUE)
    
    ok <- is.finite(row_sds) & row_sds > 0
    mat_scaled <- matrix(NA_real_, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    
    mat_scaled[ok, ] <- sweep(
      sweep(mat[ok, , drop = FALSE], 1, row_means[ok], "-"),
      1, row_sds[ok], "/"
    )
    
    # Final cleanup + unique feature names
    if (!is.null(rownames(mat_scaled))) {
      rownames(mat_scaled) <- paste0(rownames(mat_scaled), "_", view_name)
    }
    
    mat_scaled <- mat_scaled[rowSums(!is.na(mat_scaled)) > 1, , drop = FALSE]
    
    result[[view_name]] <- mat_scaled
    
    message(sprintf(
      "View '%s': %d features retained (%s)",
      view_name, nrow(mat_scaled),
      if (limma_ok) "limma F" else "row variance fallback"
    ))
  }
  
  result
}



######################  Load and run MOFA #######################

initialise_and_run_mofa <- function(mofa_list, output_file = "Results/Full_MOFA_model.hdf5") {
  
  if (!dir.exists("Results/")) dir.create("Results/")
  
  print("Initialising and running MOFA...")
  
  message("Creating MOFA object...")
  MOFAobject <- create_mofa(mofa_list)
  
  message("Setting MOFA options...")
  data_opts  <- get_default_data_options(MOFAobject)
  model_opts <- get_default_model_options(MOFAobject)
  
  if(model_opts$num_factors >= 6){
    model_opts$num_factors <- 6
  }
  
  train_opts <- get_default_training_options(MOFAobject)
  
  train_opts$drop_factor_threshold <- 0.01
  
  train_opts$convergence_mode <- "slow"
  train_opts$maxiter <- 10000
  
  mefisto_opts <- get_default_mefisto_options(MOFAobject)
  
  
  message("Preparing MOFA model...")
  MOFAobject <- prepare_mofa(
    MOFAobject,
    data_options     = data_opts,
    mefisto_options = mefisto_opts,
    model_options    = model_opts,
    training_options = train_opts
  )
  
  if (file.exists(output_file)) {
    message("Loading existing MOFA model from file...")
    MOFAobject.trained <- load_model(output_file)
  } else {
    message("No existing MOFA model found. Training new model...")
    MOFAobject.trained <- run_mofa(MOFAobject, outfile = output_file, use_basilisk = TRUE)
    message("MOFA training complete.")
  }
  
  message("Annotating metadata with sample categories...")
  MOFAobject.trained@samples_metadata$category <- as.factor(MOFAobject.trained@samples_metadata$group)
  
  return(MOFAobject.trained)
}

######################  Post-process #######################

message("Generating kegg DF...")

# ---- Extract KEGG-only weights per species ----
extract_kegg_weights <- function(MOFAobject.trained,
                                 view_metadata,
                                 metabolite_annotation = NULL,
                                 gene_annotation = NULL) {
  message("Extracting and scaling MOFA factor weights...")
  MOFAobject.trained@expectations$adjusted_W <- get_weights(MOFAobject.trained, scale = TRUE)
  print("MOFA initialisation and training complete.")
  
  # --- Identify transcriptomic and metabolomic sheets ---
  transcriptomic_sheets <- subset(view_metadata, Active == "Y" & Type == "RNAseq", select = View)$View
  metabolomic_views   <- subset(view_metadata, Active == "Y" & Type == "Metabolic", select = View)$View
  
  # --- KEGG annotation maps (per species) ---
  annotation_maps <- list()
  if (!is.null(metabolite_annotation)) {
    message("Using metabolite annotations (species-specific).")
    annotation_maps$metabolite <- metabolite_annotation
  }
  if (!is.null(gene_annotation)) {
    message("Using gene annotations (species-specific).")
    annotation_maps$gene <- gene_annotation
  }
  
  # --- Prepare output slot ---
  MOFAobject.trained@expectations$KEGG_weights <- list()
  
  # --- Filter weights with KEGG IDs, restricted by view<->species mapping ---
  message("Extracting KEGG-only weights from MOFA object...")
  weights <- MOFAobject.trained@expectations$adjusted_W
  
  # Loop per view
  for (df_name in names(weights)) {
    # Which species belong to this view?
    sp_list <- unique(unlist(strsplit(
      view_metadata$Kegg_Species[view_metadata$View == df_name & view_metadata$Active == "Y"], ","
    )))
    sp_list <- trimws(na.omit(sp_list))
    
    if (length(sp_list) == 0) {
      message("Skipping view: ", df_name, " (no active species assigned).")
      next
    }
    
    df <- weights[[df_name]]
    clean_names <- toupper(trimws(gsub("_[^_]*$", "", sub("^.*\\|", "", rownames(df)))))
    row_map <- data.frame(original = rownames(df), cleaned = clean_names, stringsAsFactors = FALSE)
    
    # For each allowed species for this view
    for (sp in sp_list) {
      ann <- NULL
      if (df_name %in% metabolomic_views && !is.null(annotation_maps$metabolite[[sp]])) {
        ann <- annotation_maps$metabolite[[sp]]
        by_col <- "Metabolite"
        ann[[by_col]] <- toupper(ann[[by_col]])
      } else if (df_name %in% transcriptomic_sheets && !is.null(annotation_maps$gene[[sp]])) {
        ann <- annotation_maps$gene[[sp]]
        by_col <- "gene"
        ann[[by_col]] <- toupper(ann[[by_col]])
      }
      
      matched_rows <- if (!is.null(ann)) {
        merge(row_map, ann, by.x = "cleaned", by.y = by_col, all.x = FALSE)
      } else NULL
      
      if (!is.null(matched_rows) && nrow(matched_rows) > 0) {
        df_filtered <- df[matched_rows$original, , drop = FALSE]
        rownames(df_filtered) <- gsub(paste0("^", sp, ":"), "", matched_rows$KEGG_ID)
        df_filtered <- unique(df_filtered)
        slot_name <- paste0(df_name, "_", sp)
        MOFAobject.trained@expectations$KEGG_weights[[slot_name]] <- df_filtered
        
        message(sprintf("Species: %s | View: %s | Matched: %d / %d",
                        sp, df_name, nrow(matched_rows), nrow(df)))
      } else {
        message(sprintf("No KEGG mapping found for species %s in view: %s", sp, df_name))
      }
    }
  }
  
  print("Finished extracting and storing KEGG-only weights per species.")
  return(MOFAobject.trained)
}


# ---- Run limma statistics per grouping variable ----
run_limma_statistics <- function(MOFAobject.trained, group.colors) {
  message("Generating statistics per grouping variable in group.colors...")
  view_names <- names(MOFAobject.trained@data)
  MOFAobject.trained@expectations$Statistics <- list()
  
  for (view_name in view_names) {
    message(sprintf("Processing view: %s", view_name))
    raw_data <- as.data.frame(MOFAobject.trained@data[[view_name]]$group1)
    
    # --- Skip empty or all-NA views ---
    if (ncol(raw_data) < 3) {
      message(sprintf("  Skipping view %s (fewer than 3 samples)", view_name))
      next
    }
    if (all(is.na(raw_data))) {
      message(sprintf("  Skipping view %s (all values are NA)", view_name))
      next
    }
    
    feature_ids <- rownames(raw_data)
    samples_meta <- MOFAobject.trained@samples_metadata
    samples_meta <- samples_meta[match(colnames(raw_data), samples_meta$sample), ]
    
    for (meta_col in names(group.colors)) {
      if (!meta_col %in% colnames(samples_meta)) {
        message(sprintf("  Skipping %s (not found in metadata)", meta_col))
        next
      }
      
      # Subset to samples with non-NA covariate and non-NA data
      keep_samples <- !is.na(samples_meta[[meta_col]])
      sub_raw <- raw_data[, keep_samples, drop = FALSE]
      sub_meta <- samples_meta[keep_samples, , drop = FALSE]
      
      # --- Skip if all features NA after filtering ---
      if (all(is.na(sub_raw))) {
        message(sprintf("  Skipping %s (all values NA after filtering)", meta_col))
        next
      }
      
      # Remove features that are entirely NA
      sub_raw <- sub_raw[rowSums(!is.na(sub_raw)) > 0, , drop = FALSE]
      if (nrow(sub_raw) == 0 || ncol(sub_raw) < 3) {
        message(sprintf("  Skipping %s (insufficient data after NA removal)", meta_col))
        next
      }
      
      message(sprintf("  Grouping variable: %s (%d samples)", meta_col, ncol(sub_raw)))
      group_labels <- factor(sub_meta[[meta_col]], levels = names(group.colors[[meta_col]]))
      
      if (nlevels(group_labels) < 2 || ncol(sub_raw) <= nlevels(group_labels)) {
        message(sprintf("  Skipping %s (insufficient residual degrees of freedom)", meta_col))
        next
      }
      
      # --- Design and fitting ---
      safe_levels <- make.names(levels(group_labels))
      design <- model.matrix(~ 0 + group_labels)
      colnames(design) <- safe_levels
      
      fit <- tryCatch(
        limma::lmFit(sub_raw, design),
        error = function(e) {
          message(sprintf("  Skipping %s (lmFit failed: %s)", meta_col, e$message))
          return(NULL)
        }
      )
      if (is.null(fit)) next
      
      key <- paste(view_name, meta_col, sep = "_")
      results_df <- data.frame(row.names = rownames(sub_raw))
      
      tt <- NULL
      if (nlevels(group_labels) == 2) {
        contrast_matrix <- limma::makeContrasts(
          contrasts = paste0(safe_levels[1], " - ", safe_levels[2]),
          levels = design
        )
        fit2 <- tryCatch(
          {
            fit2 <- limma::contrasts.fit(fit, contrast_matrix)
            limma::eBayes(fit2)
          },
          error = function(e) {
            message(sprintf("  Skipping %s (eBayes failed: %s)", meta_col, e$message))
            return(NULL)
          }
        )
        if (is.null(fit2)) next
        tt <- suppressWarnings(limma::topTable(fit2, coef = 1, number = Inf, sort.by = "none"))
      } else {
        fit2 <- tryCatch(
          limma::eBayes(fit),
          error = function(e) {
            message(sprintf("  Skipping %s (eBayes failed: %s)", meta_col, e$message))
            return(NULL)
          }
        )
        if (is.null(fit2)) next
        tt <- suppressWarnings(limma::topTable(fit2, coef = 1:ncol(design), number = Inf, sort.by = "none", p.value = 1))
      }
      
      # --- Extract results if available ---
      results_df$logFC <- NA
      results_df$Pval  <- NA
      results_df$FDR   <- NA
      if (!is.null(tt) && nrow(tt) > 0) {
        if ("logFC" %in% colnames(tt)) results_df[rownames(tt), "logFC"] <- tt$logFC
        if ("P.Value" %in% colnames(tt)) results_df[rownames(tt), "Pval"] <- tt$P.Value
        if ("adj.P.Val" %in% colnames(tt)) results_df[rownames(tt), "FDR"] <- tt$adj.P.Val
      }
      colnames(results_df) <- paste(paste(levels(group_labels), collapse = " vs "), colnames(results_df))
      
      keep <- !is.na(results_df[, 3]) & results_df[, 3] <= 0.05
      results_df <- results_df[keep, , drop = FALSE]
      
      MOFAobject.trained@expectations$Statistics[[key]] <- results_df
    }
  }
  
  message("Finished generating statistics for all views and grouping variables.")
  return(MOFAobject.trained)
}


######################  Add metadata to the experiment #######################
generate_sample_metadata <- function(sample_names) {
  
  if (!is.null(sample_metadata)) {
    for (col in colnames(sample_metadata)) {
      sample_metadata[[col]] <- as.factor(sample_metadata[[col]])
    }
    
    return(sample_metadata)
  } else {
    message("Generating metadata from sample names...")
    metadata <- data.frame(matrix(ncol = 0, nrow = length(sample_names)))
    rownames(metadata) <- sample_names
    metadata$sample <- rownames(metadata)
    metadata$category <- substr(metadata$sample, 1, nchar(metadata$sample) - 3)
    metadata$category <- factor(metadata$category, levels = unique(metadata$category))
    return(metadata)
  }
  
}

assign_sample_metadata_to_mofa <- function(MOFAobject.trained, metadata) {
  
  
  message("Assigning group and binary indicators to MOFA object...")
  MOFAobject.trained@samples_metadata$group <- metadata$category
  for (cat in levels(metadata$category)) {
    metadata[[cat]] <- as.integer(metadata$category == cat)
  }
  samples_metadata(MOFAobject.trained) <- metadata
  return(MOFAobject.trained)
}

define_group_colours <- function(metadata,
                                 base.colors = c(
                                   "#bcf60c", "#3cb44b", "#ffe119", "#f58231", "#99c0f0",
                                   "#46f0f0", "#911eb4", "#3a13eb", "#e6194b", "#f032e6",
                                   "#fabebe", "#008080", "#e6beff", "#9a6324", "#fffac8"
                                 ),
                                 skip_cols = c("sample", "group"),
                                 rotate_by = c("name", "index"),
                                 max_levels = 50L,
                                 numeric_like_prop = 0.9) {
  
  rotate_by <- match.arg(rotate_by)
  nbase <- length(base.colors)
  
  offset_for <- function(col_name, k) {
    if (rotate_by == "index") (k - 1) %% nbase else (sum(utf8ToInt(col_name)) %% nbase)
  }
  
  # TRUE if a (factor/character) column is really continuous numbers stored as text
  is_numeric_like_text <- function(x, prop = 0.9) {
    if (!(is.factor(x) || is.character(x))) return(FALSE)
    
    y <- as.character(x)
    y <- y[!is.na(y)]
    if (length(y) == 0L) return(TRUE)
    
    y <- trimws(y)
    y <- y[nzchar(y)]
    if (length(y) == 0L) return(TRUE)
    
    suppressWarnings(num <- as.numeric(y))
    mean(!is.na(num)) >= prop
  }
  
  group.colors <- list()
  cols <- setdiff(names(metadata), skip_cols)
  
  for (k in seq_along(cols)) {
    col <- cols[k]
    x <- metadata[[col]]
    
    # ignore true numeric columns outright
    if (is.numeric(x)) next
    
    # only consider categorical text columns
    if (!(is.factor(x) || is.character(x))) next
    
    # drop "continuous disguised as text" (your current failure mode)
    if (is_numeric_like_text(x, prop = numeric_like_prop)) next
    
    lvls <- if (is.factor(x)) levels(x) else unique(x[!is.na(x)])
    lvls <- as.character(lvls)
    lvls <- lvls[!is.na(lvls) & nzchar(trimws(lvls))]
    n <- length(lvls)
    
    if (n <= 1L) next
    if (!is.null(max_levels) && n > max_levels) next
    
    off <- offset_for(col, k)
    idx <- ((seq_len(n) + off - 1) %% nbase) + 1
    colset <- base.colors[idx]
    names(colset) <- lvls
    
    group.colors[[col]] <- colset
  }
  
  group.colors
}

define_view_colours <- function(view_names) {
  category.colors <- c(
    "#bcf60c",
    "#3cb44b",
    "#ffe119",
    "#f58231",
    "#99c0f0",
    "#46f0f0",
    "#911eb4",
    "#3a13eb",
    "#e6194b",
    "#f032e6",
    "#fabebe",
    "#008080",
    "#e6beff",
    "#9a6324",
    "#fffac8"
  )
  n <- length(view_names)
  view.colors <- category.colors[1:n]
  names(view.colors) <- view_names
  return(view.colors)
}

standardise_feature_names <- function(MOFAobject.trained, max_length = 30) {
  message("Standardising feature names...")
  
  truncate_name <- function(x) {
    ifelse(
      nchar(x) > max_length,
      paste0(substr(x, 1, max_length - 3), "..."),
      x
    )
  }
  
  # Create cleanNames without view suffix
  MOFAobject.trained@features_metadata$cleanNames <- mapply(
    function(f, v) sub(paste0("_", v, "$"), "", f),
    MOFAobject.trained@features_metadata$feature,
    MOFAobject.trained@features_metadata$view
  )
  
  # Apply truncation
  MOFAobject.trained@features_metadata$truncatedNames <- 
    truncate_name(MOFAobject.trained@features_metadata$cleanNames)
  
  # Copy
  MOFAobject.trained.cleanNames <- MOFAobject.trained
  
  # Update expectations$W
  MOFAobject.trained.cleanNames@expectations$W <- lapply(
    names(MOFAobject.trained.cleanNames@expectations$W),
    function(view) {
      W <- MOFAobject.trained.cleanNames@expectations$W[[view]]
      if (!is.null(W) && (is.matrix(W) || is.data.frame(W) || inherits(W, "Matrix"))) {
        rn <- rownames(W)
        rn <- sub(paste0("_", view, "$"), "", rn)
        rn <- truncate_name(rn)
        rownames(W) <- rn
      }
      W
    }
  )
  names(MOFAobject.trained.cleanNames@expectations$W) <- 
    names(MOFAobject.trained@expectations$W)
  
  # Update data matrices (nested lists: view -> group)
  MOFAobject.trained.cleanNames@data <- lapply(
    names(MOFAobject.trained@data),
    function(view) {
      view_list <- MOFAobject.trained@data[[view]]
      if (!is.null(view_list) && is.list(view_list)) {
        view_list <- lapply(view_list, function(D) {
          if (!is.null(D) && (is.matrix(D) || is.data.frame(D) || inherits(D, "Matrix"))) {
            rn <- rownames(D)
            rn <- sub(paste0("_", view, "$"), "", rn)
            rn <- truncate_name(rn)
            rownames(D) <- rn
          }
          D
        })
      }
      view_list
    }
  )
  names(MOFAobject.trained.cleanNames@data) <- names(MOFAobject.trained@data)
  
  # Replace feature with truncated names
  MOFAobject.trained.cleanNames@features_metadata$feature <- 
    MOFAobject.trained.cleanNames@features_metadata$truncatedNames
  
  return(MOFAobject.trained.cleanNames)
}





########################################
# VARIANCE PLOTS
########################################

# ---- Variance explained plots ----
generate_mofa_variance_plots <- function(MOFAobject.trained,
                                         view.colors,
                                         outdir = "Results/Corrilation_and_Variance_plots") {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  output_file <- file.path(outdir, "Variance_explained.png")
  if (file.exists(output_file)) {
    message("File already exists: ", output_file)
    return(invisible(NULL))
  }
  
  message("Generating variance explained plots from MOFA object...")
  
  # Get variance data
  variance_data <- get_variance_explained(MOFAobject.trained)
  
  # -------------------------
  # Factor-wise variance
  # -------------------------
  group_vs_view_variance_df <- as.data.frame(variance_data$r2_per_factor$group1)
  number_of_factors <- nrow(group_vs_view_variance_df)
  group_vs_view_variance_df$Factor <- rownames(group_vs_view_variance_df)
  
  variance_long <- reshape2::melt(
    group_vs_view_variance_df,
    id.vars = "Factor",
    variable.name = "View",
    value.name = "Variance"
  )
  variance_long$Variance <- variance_long$Variance / number_of_factors
  
  # -------------------------
  # View-level variance
  # -------------------------
  view_variance_df <- as.data.frame(variance_data$r2_total$group1)
  view_variance_df$View <- rownames(view_variance_df)
  
  view_variance_long <- reshape2::melt(
    view_variance_df,
    id.vars = "View",
    variable.name = "Factor",
    value.name = "VarianceExplained"
  )
  view_variance_long$Variance <- view_variance_long$VarianceExplained / number_of_factors
  
  # Ensure consistent order
  ordered_views <- sort(unique(variance_long$View))
  variance_long$View <- factor(variance_long$View, levels = ordered_views)
  view_variance_long$View <- factor(view_variance_long$View, levels = ordered_views)
  
  # -------------------------
  # Plot 1: Heatmap of factor-wise variance
  # -------------------------
  Totalvariance <- ggplot2::ggplot(variance_long, ggplot2::aes(x = View, y = Factor, fill = Variance)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = "white", high = "purple", name = "Variance (R²)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14, colour = "black"),
      axis.text.y = ggplot2::element_text(size = 12, colour = "black"),
      axis.text.x = ggplot2::element_text(size = 12, angle = 90, hjust = 1, face = "bold", colour = "black"),
      axis.title = ggplot2::element_text(size = 15, face = "bold", colour = "black"),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(title = " ", x = "", y = "MOFA Factors")
  
  # -------------------------
  # Plot 2: Barplot of total variance per view
  # -------------------------
  wheatvarExp <- ggplot2::ggplot(view_variance_long, ggplot2::aes(x = View, y = Variance, fill = View)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(values = view.colors) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14, colour = "black"),
      axis.text.y = ggplot2::element_text(size = 12, colour = "black"),
      axis.text.x = ggplot2::element_text(size = 12, angle = 90, hjust = 1, face = "bold", colour = "black"),
      axis.title = ggplot2::element_text(size = 15, face = "bold", colour = "black"),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = "Variance Explained by MOFA Views and Factors",
      x = "", y = "Variance (R²)", fill = "View"
    )
  
  # -------------------------
  # Combine and save
  # -------------------------
  message("Saving combined variance plot to: ", output_file)
  combined_plot <- gridExtra::grid.arrange(wheatvarExp, Totalvariance, ncol = 2, nrow = 1)
  
  ggplot2::ggsave(filename = output_file, plot = combined_plot, width = 10, height = 5)
  
  message("Variance explained plots saved.")
}



########################################
# PCA PLOTS
########################################

generate_mofa_factor_matrix_plot <- function(MOFAobject.trained,
                                             view.colors,
                                             group.colors,
                                             output_dir = "Results/PCA",
                                             ridge = 1e-6,
                                             ellipse_scale = 2,
                                             padding = 1.05,
                                             min_axis = 1e-2) {
  
  output_dir
  output_dir_sample <- paste0(output_dir,"/individual_sample_PCA_plots")
  output_dir_data <- paste0(output_dir,"/individual_data_loading_plots")
  output_dir_histo <- paste0(output_dir,"/individual_histograms")
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!dir.exists(output_dir_sample))  dir.create(output_dir_sample, recursive = TRUE)
  if (!dir.exists(output_dir_data))  dir.create(output_dir_data, recursive = TRUE)
  if (!dir.exists(output_dir_histo))  dir.create(output_dir_histo, recursive = TRUE)
  
  message("Generating PCA matrix of MOFA factors for each grouping variable...")
  
  factor_numbers <- 1:ncol(get_weights(MOFAobject.trained)[[1]])
  weights <- get_weights(MOFAobject.trained)
  
  combined_weights <- do.call(rbind, lapply(names(weights), function(view_name) {
    df <- as.data.frame(weights[[view_name]])
    df$Feature <- rownames(df)
    df$View <- view_name
    df
  }))
  
  melted_weights <- reshape2::melt(
    combined_weights,
    id.vars = c("Feature", "View"),
    variable.name = "Factor",
    value.name = "Weight"
  )
  
  feature_matrix <- reshape2::dcast(
    melted_weights,
    Feature + View ~ Factor,
    value.var = "Weight"
  )
  
  FactorZ <- MOFAobject.trained@expectations$Z$group1
  FactorZ <- as.data.frame(FactorZ)
  FactorZ$sample <- rownames(FactorZ)
  FactorZ <- dplyr::inner_join(metadata, FactorZ, by = c("sample" = "sample"))
  
  for (meta_col in names(group.colors)) {
    message("Processing grouping variable: ", meta_col)
    
    plot_matrix <- matrix(list(),
                          nrow = length(factor_numbers),
                          ncol = length(factor_numbers))
    
    for (i in factor_numbers) {
      for (j in factor_numbers) {
        
        # ---------- UPPER TRIANGLE ----------
        if (i < j) {
          
          plot_df <- FactorZ[, c("sample", meta_col,
                                 paste0("Factor", i),
                                 paste0("Factor", j))]
          colnames(plot_df)[2] <- "Group"
          plot_df <- na.omit(plot_df)
          plot_df$Group <- factor(plot_df$Group,
                                  levels = names(group.colors[[meta_col]]))
          
          ellipse_df <- plot_df |>
            dplyr::group_by(Group) |>
            dplyr::group_modify(~{
              X <- as.matrix(.x[, c(paste0("Factor", i), paste0("Factor", j))])
              if (nrow(X) < 3) return(data.frame(x = numeric(0), y = numeric(0)))
              
              mu <- colMeans(X)
              S  <- stats::cov(X, use = "pairwise.complete.obs") + diag(ridge, 2)
              eig <- eigen(S, symmetric = TRUE)
              
              # enforce minimal variance to avoid line-collapse
              axes <- sqrt(pmax(eig$values, min_axis))
              theta <- seq(0, 2*pi, length.out = 200)
              circle <- cbind(cos(theta), sin(theta))
              
              A <- eig$vectors %*% diag(axes) %*% t(eig$vectors)
              pts <- t(mu + ellipse_scale * padding * t(circle %*% A))
              
              data.frame(x = pts[, 1], y = pts[, 2])
            }) |>
            dplyr::ungroup()
          
          plot_matrix[[i, j]] <- ggplot(
            plot_df,
            aes_string(x = paste0("Factor", i),
                       y = paste0("Factor", j),
                       color = "Group",
                       fill  = "Group")
          ) +
            geom_point(size = 4) +
            geom_polygon(
              data = ellipse_df,
              aes(x = x, y = y, color = Group, fill = Group),
              inherit.aes = FALSE,
              alpha = 0.5,
              show.legend = FALSE
            ) +
            scale_color_manual(values = group.colors[[meta_col]]) +
            scale_fill_manual(values = group.colors[[meta_col]]) +
            labs(title = paste("Sample PCA:", meta_col,
                               "- Factor", i, "vs", j)) +
            scale_x_continuous(expand = expansion(mult = 0.15)) +
            scale_y_continuous(expand = expansion(mult = 0.15)) +
            coord_cartesian(clip = "off") +
            theme_classic() +
            ggplot2::theme(
              plot.title = ggplot2::element_text(face = "bold", size = 14, colour = "black"),
              axis.text.y = ggplot2::element_text(size = 12, colour = "black"),
              axis.text.x = ggplot2::element_text(size = 12, angle = 90, hjust = 1, face = "bold", colour = "black"),
              axis.title = ggplot2::element_text(size = 15, face = "bold", colour = "black"),
              panel.grid = ggplot2::element_blank()
            )
          
          # save individual (UPPER) plot -> output_dir_data
          indiv_file <- file.path(
            output_dir_data,
            paste0(meta_col, "_", i, "_vs_", j, ".png")
          )
          ggplot2::ggsave(indiv_file, plot = plot_matrix[[i, j]],
                          width = 5, height = 5, units = "in", dpi = 150, device = "png")
          
          # ---------- DIAGONAL ----------
        } else if (i == j) {
          
          factor_values <- as.data.frame(
            get_factors(MOFAobject.trained)$group1[, i]
          )
          colnames(factor_values) <- "FactorValue"
          factor_values$Sample <- rownames(factor_values)
          factor_values$SampleType <- stringr::str_replace(
            factor_values$Sample, "_\\d+$", ""
          )
          
          plot_matrix[[i, j]] <- ggplot(
            factor_values,
            aes(x = FactorValue, fill = SampleType)
          ) +
            geom_density(alpha = 0.6) +
            scale_fill_manual(values = group.colors[[meta_col]]) +
            theme_classic() +
            labs(title = paste("Density Plot:", meta_col,
                               "- Factor", i),
                 x = "Factor Value", y = "Density") +
            ggplot2::theme(
              plot.title = ggplot2::element_text(face = "bold", size = 14, colour = "black"),
              axis.text.y = ggplot2::element_text(size = 12, colour = "black"),
              axis.text.x = ggplot2::element_text(size = 12, angle = 90, hjust = 1, face = "bold", colour = "black"),
              axis.title = ggplot2::element_text(size = 15, face = "bold", colour = "black"),
              panel.grid = ggplot2::element_blank()
            )
          
          # save individual (DIAGONAL) plot -> output_dir_histo
          indiv_file <- file.path(
            output_dir_histo,
            paste0(meta_col, i, "_density.png")
          )
          ggplot2::ggsave(indiv_file, plot = plot_matrix[[i, j]],
                          width = 5, height = 5, units = "in", dpi = 150, device = "png")
          
          # ---------- LOWER TRIANGLE ----------
        } else {
          
          weight_data <- feature_matrix[, c("Feature", "View",
                                            paste0("Factor", i),
                                            paste0("Factor", j))]
          colnames(weight_data) <- c("Feature", "View", "FactorX", "FactorY")
          
          ellipse_df <- weight_data |>
            dplyr::group_by(View) |>
            dplyr::group_modify(~{
              X <- as.matrix(.x[, c("FactorX", "FactorY")])
              if (nrow(X) < 3) return(data.frame(x = numeric(0), y = numeric(0)))
              
              mu <- colMeans(X)
              S  <- stats::cov(X, use = "pairwise.complete.obs") + diag(ridge, 2)
              eig <- eigen(S, symmetric = TRUE)
              
              axes <- sqrt(pmax(eig$values, min_axis))
              theta <- seq(0, 2*pi, length.out = 200)
              circle <- cbind(cos(theta), sin(theta))
              
              A <- eig$vectors %*% diag(axes) %*% t(eig$vectors)
              pts <- t(mu + ellipse_scale * padding * t(circle %*% A))
              
              data.frame(x = pts[, 1], y = pts[, 2])
            }) |>
            dplyr::ungroup()
          
          plot_matrix[[i, j]] <- ggplot(
            weight_data,
            aes(x = FactorX, y = FactorY, color = View)
          ) +
            geom_point(size = 2, alpha = 0.7) +
            geom_polygon(
              data = ellipse_df,
              aes(x = x, y = y, color = View, fill = View),
              inherit.aes = FALSE,
              alpha = 0.5,
              show.legend = FALSE
            ) +
            scale_color_manual(values = view.colors) +
            scale_fill_manual(values = view.colors) +
            theme_classic() +
            labs(title = paste("Feature Weights: Factor", i, "vs", j),
                 x = paste("Factor", i),
                 y = paste("Factor", j)) +
            ggplot2::theme(
              plot.title = ggplot2::element_text(face = "bold", size = 14, colour = "black"),
              axis.text.y = ggplot2::element_text(size = 12, colour = "black"),
              axis.text.x = ggplot2::element_text(size = 12, angle = 90, hjust = 1, face = "bold", colour = "black"),
              axis.title = ggplot2::element_text(size = 15, face = "bold", colour = "black"),
              panel.grid = ggplot2::element_blank()
            )
          
          # save individual (LOWER) plot -> output_dir_sample
          indiv_file <- file.path(
            output_dir_sample,
            paste0(meta_col, "_", i, "_vs_", j, "_weights.png")
          )
          ggplot2::ggsave(indiv_file, plot = plot_matrix[[i, j]],
                          width = 5, height = 5, units = "in", dpi = 150, device = "png")
        }
      }
    }
    
    out_file <- file.path(
      output_dir,
      paste0(meta_col, "_Factor_Pair_PCA_Matrix.png")
    )
    
    message("Saving PCA plot matrix for ", meta_col, " to ", out_file)
    
    png(out_file, width = 30, height = 20, units = "in", res = 300)
    gridExtra::grid.arrange(
      grobs = as.vector(plot_matrix),
      nrow = length(factor_numbers),
      ncol = length(factor_numbers)
    )
    tryCatch(while (dev.cur() > 1) dev.off(), error = function(e) {})
  }
  
  message("All PCA matrices saved successfully.")
}



########################################
###FACTOR COMPARISONS
########################################

generate_mofa_factor_comparison_plot <- function(
    MOFAobject.trained,
    group.colors,
    view.colors,
    metadata,
    output_dir = "Results/Covariate_factor_comparisons"
) {
  message("Creating MOFA factor comparison plots (one per grouping variable)...")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # additional output dirs (individual plots)
  output_dir_data   <- paste0(output_dir, "/individual_violin_plots")
  output_dir_sample <- paste0(output_dir, "/individual_data_loading_plots")
  
  if (!dir.exists(output_dir_data))   dir.create(output_dir_data, recursive = TRUE)
  if (!dir.exists(output_dir_sample)) dir.create(output_dir_sample, recursive = TRUE)
  
  # Factors present (assumes same number across views)
  factor_numbers <- seq_len(ncol(get_weights(MOFAobject.trained)[[1]]))
  
  .blank_plot <- function(txt = "No data") {
    ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = txt)
  }
  
  for (meta_col in names(group.colors)) {
    message("Processing grouping variable: ", meta_col)
    
    plot_matrix <- matrix(vector("list", length(factor_numbers) * 2L),
                          nrow = 2, ncol = length(factor_numbers))
    
    for (j in factor_numbers) {
      message(sprintf("Processing factor %d...", j))
      
      Zj <- as.data.frame(get_factors(MOFAobject.trained)$group1[, j, drop = FALSE])
      colnames(Zj) <- "FactorValue"
      Zj$sample <- rownames(Zj)
      
      group_df <- metadata[, c("sample", meta_col), drop = FALSE]
      factor_values <- merge(Zj, group_df, by = "sample", sort = FALSE)
      
      factor_values$Group <- factor(
        factor_values[[meta_col]],
        levels = names(group.colors[[meta_col]])
      )
      factor_values <- factor_values[
        !is.na(factor_values$FactorValue) & !is.na(factor_values$Group),
        , drop = FALSE
      ]
      factor_values$Group <- droplevels(factor_values$Group)
      
      ## Row 1: Violin + jitter (UPPER) -> output_dir_data
      if (nrow(factor_values) > 0L && nlevels(factor_values$Group) >= 1L) {
        p1 <- ggplot2::ggplot(
          factor_values,
          ggplot2::aes(x = Group, y = FactorValue, fill = Group)
        ) +
          ggplot2::geom_violin(trim = FALSE, alpha = 0.6) +
          ggplot2::geom_jitter(shape = 16, position = ggplot2::position_jitter(0.2),
                               alpha = 0.5, size = 1) +
          ggplot2::scale_fill_manual(values = group.colors[[meta_col]]) +
          ggplot2::theme_classic() +
          ggplot2::labs(
            title = paste("Factor", j, "-", meta_col),
            x = meta_col, y = "Factor Value"
          ) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 14, colour = "black"),
            axis.text  = ggplot2::element_text(size = 8, colour = "black"),
            axis.title = ggplot2::element_text(size = 10, face = "bold", colour = "black"),
            legend.text  = ggplot2::element_text(size = 8, colour = "black"),
            legend.title = ggplot2::element_text(size = 10, colour = "black"),
            panel.grid = ggplot2::element_blank()
          )
      } else {
        p1 <- .blank_plot(paste0("Factor ", j, " — no groups/data"))
      }
      plot_matrix[[1, j]] <- p1
      
      indiv_file_p1 <- file.path(
        output_dir_data,
        paste0(meta_col, "_Factor", j, "_violin.png")
      )
      ggplot2::ggsave(indiv_file_p1, plot = p1,
                      width = 5, height = 4, units = "in", dpi = 150, device = "png")
      
      ## Row 2: Mean vs |weight| (LOWER) -> output_dir_sample
      views <- views_names(MOFAobject.trained)
      scatter_list <- vector("list", length(views))
      
      for (vi in seq_along(views)) {
        view <- views[vi]
        
        if (nrow(factor_values) == 0L) next
        grp <- factor_values$Group
        if (nlevels(grp) == 0L) next
        
        raw_data_view <- MOFAobject.trained@data[[view]]$group1
        
        idx <- match(factor_values$sample, colnames(raw_data_view))
        if (any(is.na(idx))) next
        sub_raw <- raw_data_view[, idx, drop = FALSE]
        
        w_view <- MOFAobject.trained@expectations$W[[view]][, j]
        feats  <- intersect(rownames(sub_raw), names(w_view))
        if (length(feats) == 0L) next
        sub_raw <- sub_raw[feats, , drop = FALSE]
        w_view  <- w_view[feats]
        
        if (nlevels(grp) >= 2L) {
          G <- stats::model.matrix(~ -1 + grp)
          grp_names <- colnames(G)
        } else {
          G <- matrix(1, nrow = length(grp), ncol = 1)
          grp_names <- if (nlevels(grp) == 1L) levels(grp) else "Group"
          colnames(G) <- grp_names
        }
        
        M <- t(sub_raw)
        counts <- colSums(G)
        means <- (t(G) %*% M) / pmax(counts, 1)
        rownames(means) <- grp_names
        
        mean_df <- as.data.frame(t(means))
        mean_df$Feature <- rownames(mean_df)
        long <- reshape2::melt(
          mean_df, id.vars = "Feature",
          variable.name = "Group", value.name = "FeatureMean"
        )
        long$Weight <- abs(w_view[long$Feature])
        long$view   <- view
        long$Group  <- factor(long$Group, levels = grp_names)
        
        scatter_list[[vi]] <- long
      }
      
      scatter_list <- Filter(function(x) !is.null(x) && nrow(x) > 0, scatter_list)
      
      if (length(scatter_list) == 0L) {
        p2 <- .blank_plot(paste0("Factor ", j, " — no view/feature data"))
      } else {
        scatter_data <- do.call(rbind, scatter_list)
        p2 <- ggplot2::ggplot(
          scatter_data,
          ggplot2::aes(x = FeatureMean, y = Weight, colour = view)
        ) +
          ggplot2::geom_point(alpha = 0.8, size = 2) +
          ggplot2::theme_classic() +
          ggplot2::scale_color_manual(values = view.colors) +
          ggplot2::facet_wrap(~ Group, scales = "free_x", nrow = 1) +
          ggplot2::labs(
            title = "Mean Feature Values",
            y = paste("Factor", j, "Weight")
          ) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 14, colour = "black"),
            axis.text  = ggplot2::element_text(size = 8, colour = "black"),
            axis.title = ggplot2::element_text(size = 10, face = "bold", colour = "black"),
            legend.text  = ggplot2::element_text(size = 8, colour = "black"),
            legend.title = ggplot2::element_text(size = 10, colour = "black"),
            panel.grid = ggplot2::element_blank()
          )
      }
      plot_matrix[[2, j]] <- p2
      
      indiv_file_p2 <- file.path(
        output_dir_sample,
        paste0(meta_col, "_Factor", j, "_mean_vs_weight.png")
      )
      ggplot2::ggsave(indiv_file_p2, plot = p2,
                      width = 6, height = 4, units = "in", dpi = 150, device = "png")
    }
    
    message("Assembling grid for ", meta_col, "...")
    plot_list <- as.vector(plot_matrix)
    
    grid_plot <- gridExtra::arrangeGrob(
      grobs = plot_list,
      nrow = length(factor_numbers),
      ncol = 2,
      top = grid::textGrob(paste("FACTOR vs", meta_col),
                           gp = grid::gpar(fontsize = 14, fontface = "bold")),
      padding = grid::unit(30, "pt")
    )
    
    out_file <- file.path(output_dir, paste0(meta_col, "_Factor_Comparison.png"))
    message("Saving factor comparison figure for ", meta_col, " to ", out_file)
    
    dynamic_width <- length(group.colors[[meta_col]]) * 2 + 6
    grDevices::png(out_file, width = dynamic_width, height = 20, units = "in", res = 300)
    grid::grid.draw(grid_plot)
    grDevices::dev.off()
  }
}


########################################
# FACTOR CORRELATION PLOTS
########################################

generate_mofa_factor_correlation_plot <- function(MOFAobject.trained,
                                                  view.colors,
                                                  sample_metadata,
                                                  output_dir = "Results/Corrilation_and_Variance_plots") {
  message("Generating MOFA factor correlation plots...")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 0) Keep the object intact; do not drop 'group' in samples_metadata
  md <- samples_metadata(MOFAobject.trained)
  
  # 1) Align metadata to sample order (defensive but useful)
  samps <- as.character(unlist(samples_names(MOFAobject.trained)))
  md <- md[match(samps, md$sample), , drop = FALSE]
  
  # 2) Choose covariate columns from metadata, explicitly excluding 'sample' and 'group'
  covar_cols <- setdiff(colnames(md), c("sample", "group"))
  
  # (optional) Pre-coerce selected columns to numeric to avoid the “non-numeric… converting” warning
  md[ covar_cols ] <- lapply(md[ covar_cols ], function(x) {
    y <- if (is.numeric(x) || is.integer(x)) x else as.numeric(factor(x))
    # optional: impute missing to 0 for your design
    y[is.na(y)] <- 0
    y
  })
  
  # (optional) Drop columns with zero variance across all samples (avoids sd=0 from cor())
  covar_cols <- covar_cols[ vapply(md[covar_cols], function(x) { v <- var(x, na.rm = TRUE); is.finite(v) && v > 0 }, logical(1)) ]
  
  core_covar_cols <- unique(c("category", as.character(sample_metadata$category)))
  
  # ---- Plot 1: Factor vs Factor ----
  output_file <- file.path(output_dir, "Factor_vs_Sample_corrilation.png")
  png(filename = output_file, width = 12, height = 6, units = "in", res = 300)
  
  message("Rendering sample-to-factor correlation heatmap...")
  correlate_factors_with_covariates(
    MOFAobject.trained,
    covariates = core_covar_cols,   # <- character vector of cols from samples_metadata
    groups     = "group1",     # set explicitly if single-group model
    plot       = "r",
    transpose = TRUE,
    mar = c(2, 2, 2, 2),
    tl.cex = 1.2,
    cl.ratio = 0.25
  )
  
  tryCatch(while (dev.cur() > 1) dev.off(), error = function(e) {})
  message(sprintf("Saved Factor vs Sample plot to %s", output_file))
  
  other_covar_cols <- setdiff(covar_cols, core_covar_cols)
  
  if(length(other_covar_cols) >= 1){
    # ---- Plot 2: Factor vs Covariate ----
    output_file <- file.path(output_dir, "Factor_vs_Covariate_corrilation.png")
    png(filename = output_file, width = 12, height = 6, units = "in", res = 300)
    
    message("Rendering covariate-to-factor correlation heatmap...")
    correlate_factors_with_covariates(
      MOFAobject.trained,
      covariates = other_covar_cols,   # <- character vector of cols from samples_metadata
      groups     = "group1",     # set explicitly if single-group model
      plot       = "r",
      mar = c(2, 2, 2, 2),
      tl.cex = 1.2,
      cl.ratio = 0.25
    )
    
    tryCatch(while (dev.cur() > 1) dev.off(), error = function(e) {})
    message(sprintf("Saved Factor vs Covariates plot to %s", output_file))
    
    # ---- Plot 3: Factor vs Covatiate ----
    output_file <- file.path(output_dir, "Factor_vs_Covariate_corrilation_pval.png")
    png(filename = output_file, width = 12, height = 6, units = "in", res = 300)
    
    message("Rendering covariate-to-factor correlation Pval heatmap...")
    correlate_factors_with_covariates(
      MOFAobject.trained,
      covariates = other_covar_cols,   # <- character vector of cols from samples_metadata
      groups     = "group1",     # set explicitly if single-group model
      plot       = "log_p",
      transpose = TRUE,
      mar = c(2, 2, 2, 2),
      tl.cex = 1.2,
      cl.ratio = 0.25
    )
    
    tryCatch(while (dev.cur() > 1) dev.off(), error = function(e) {})
    message(sprintf("Saved Factor vs Covariates plot to %s", output_file))
    
    # ---- Table Factor vs Covatiate pvals ----
    dataframe <- correlate_factors_with_covariates(
      MOFAobject.trained,
      covariates = covar_cols,   # <- character vector of cols from samples_metadata
      groups     = "group1",     # set explicitly if single-group model
      plot       = "log_p",
      return_data = TRUE,
      alpha = 1,
      transpose = TRUE,
      mar = c(2, 2, 2, 2),
      tl.cex = 1.2,
      cl.ratio = 0.25
    )
    
    pvals <- 10^(-dataframe)
    output_file <- file.path(output_dir, "Factor_vs_Covariate_corrilation_pval.csv")
    write.csv(x = pvals,file = output_file)
    
  }
  
  # ---- Plot 4: Factor vs Factor ----
  output_file2 <- file.path(output_dir, "Factor_vs_Factor_corrilation.png")
  png(filename = output_file2, width = 8, height = 6, units = "in", res = 300)
  message("Rendering factor-to-factor correlation heatmap...")
  plot_factor_cor(
    MOFAobject.trained,
    title = "",
    mar = c(2, 2, 2, 2),
    tl.cex = 1.2,
    cl.ratio = 0.25
  )
  mtext("Factor vs Factor", side = 1, line = 4, font = 2)
  tryCatch(while (dev.cur() > 1) dev.off(), error = function(e) {})
  message(sprintf("Saved Factor vs Factor plot to %s", output_file2))
}

########################################
# PER VIEW PLOTS
########################################

# Helper to ensure output directories exist
ensure_output_dirs <- function(factor_name, view_name) {
  output_dir <- file.path("Results/Per_Factor_Analysis", factor_name, view_name)
  if (!dir.exists(output_dir))
    dir.create(output_dir, recursive = TRUE)
  return(output_dir)
}

save_weight_tables <- function(MOFAobject.trained,
                               group.colors,
                               outdir = "Results/Significance_testing/tables") {
  views <- names(MOFAobject.trained@expectations$W)
  
  # Create output root
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  for (meta_col in names(group.colors)) {
    combined_rows <- list()
    
    for (view_name in views) {
      stats_key <- paste(view_name, meta_col, sep = "_")
      stats_df  <- MOFAobject.trained@expectations$Statistics[[stats_key]]
      
      if (is.null(stats_df)) {
        message(sprintf("No statistics available for %s/%s, skipping.", view_name, meta_col))
        next
      }
      
      # Core matrices for this view
      W_view   <- MOFAobject.trained@expectations$W[[view_name]]
      adj_W    <- MOFAobject.trained@expectations$adjusted_W[[view_name]]
      raw_mat  <- MOFAobject.trained@data[[view_name]][["group1"]]
      
      feat_ids <- rownames(W_view)
      if (is.null(feat_ids)) stop("W has no rownames for view: ", view_name)
      
      # Human-friendly feature names
      feature_names <- gsub(paste0("_", view_name, "$"), "", feat_ids)
      
      # Base table
      base_df <- data.frame(
        Feature = feature_names,
        View    = view_name,
        stringsAsFactors = FALSE,
        row.names = feat_ids
      )
      
      # Stats
      stats_sub <- stats_df[feat_ids, , drop = FALSE]
      
      # Identify p-value columns
      p_cols <- grep("Pval$", colnames(stats_sub), ignore.case = TRUE, value = TRUE)
      
      
      # Convert to numeric + filter rows
      keep <- rep(TRUE, nrow(stats_sub))
      if (length(p_cols) > 0L) {
        for (pc in p_cols) {
          stats_sub[[pc]] <- suppressWarnings(as.numeric(stats_sub[[pc]]))
        }
        keep <- apply(stats_sub[, p_cols, drop = FALSE], 1, function(r) {
          all(!is.na(r)) && all(r <= 0.05)
        })
      }
      if (!any(keep)) {
        message(sprintf("All rows filtered for %s/%s.", view_name, meta_col))
        next
      }
      
      # Apply filter
      base_df   <- base_df[keep, , drop = FALSE]
      stats_sub <- stats_sub[keep, , drop = FALSE]
      
      # Add stats
      base_df <- cbind(base_df, stats_sub)
      
      # Add weights
      factor_names <- colnames(W_view)
      for (fn in factor_names) {
        base_df[[paste0("Weight_", fn)]] <- W_view[rownames(base_df), fn]
      }
      if (!is.null(adj_W)) {
        for (fn in factor_names) {
          base_df[[paste0("AdjWeight_", fn)]] <- adj_W[rownames(base_df), fn]
        }
      }
      
      # Raw data last
      raw_df <- as.data.frame(raw_mat[rownames(base_df), , drop = FALSE], check.names = FALSE)
      row_df <- cbind(base_df, raw_df)
      
      combined_rows[[view_name]] <- row_df
      message(sprintf("processed %d features for %s/%s", nrow(row_df), view_name, meta_col))
    }
    
    if (length(combined_rows) == 0) {
      message(sprintf("No rows left for %s after filtering.", meta_col))
      next
    }
    
    # Combine across views
    combined_df <- do.call(rbind, combined_rows)
    
    # File per meta_col in a single folder
    out_file <- file.path(outdir, sprintf("Weights_AllViews_AllFactors__%s.csv", meta_col))
    write.csv(combined_df, out_file, row.names = FALSE)
    message(sprintf("Wrote %s (n=%d rows, %d cols)", out_file, nrow(combined_df), ncol(combined_df)))
  }
}

########################################
# REGRESSION TESTING
########################################

plot_scatter_pval_vs_weight <- function(view.colors,
                                        indir = "Results/Significance_testing/tables",
                                        outdir = "Results/Significance_testing/regression_testing") {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  # ---- additional output dirs (individual plots) ----
  output_dir_r2      <- paste0(outdir, "/individual_R2_barplots")
  output_dir_scatter <- paste0(outdir, "/individual_R2_scatter_plots")
  
  if (!dir.exists(output_dir_r2)) dir.create(output_dir_r2, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir_scatter)) dir.create(output_dir_scatter, recursive = TRUE, showWarnings = FALSE)
  
  files <- list.files(indir, pattern = "\\.csv$", full.names = TRUE)
  if (length(files) == 0L) stop("No CSV files found in ", indir)
  
  for (file in files) {
    message("Processing file: ", basename(file))
    df <- utils::read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
    if (!all(c("Feature","View") %in% names(df))) next
    
    # P-values and weights
    p_cols <- grep("Pval$", names(df), value = TRUE)
    w_prefix <- "^Weight_"
    weight_cols <- grep(w_prefix, names(df), value = TRUE)
    if (length(p_cols) == 0L || length(weight_cols) == 0L) next
    
    for (pc in p_cols) df[[pc]] <- suppressWarnings(as.numeric(df[[pc]]))
    
    # Melt weights to long format
    long <- reshape2::melt(
      df,
      id.vars = c("Feature","View", p_cols),
      measure.vars = weight_cols,
      variable.name = "FactorCol",
      value.name = "Weight"
    )
    long$Factor <- sub(w_prefix, "", long$FactorCol)
    long$Weight <- abs(long$Weight)
    
    for (pc in p_cols) {
      dat <- long[, c("Feature","View","Factor","Weight", pc)]
      names(dat)[5] <- "Pval"
      dat <- dat[!is.na(dat$Pval), , drop = FALSE]
      if (nrow(dat) == 0L) next
      
      dat$Pval <- pmax(pmin(dat$Pval, 1), .Machine$double.xmin)
      dat$neglog10P <- -log10(dat$Pval)
      
      # All factors
      factor_numbers <- sort(unique(dat$Factor))
      
      plot_list <- list()
      
      # ---- Naming base (used for both combined + individual outputs) ----
      base <- tools::file_path_sans_ext(basename(file))
      base <- sub("^Weights_AllViews_AllFactors__", "", base)
      safe_base <- gsub("[^A-Za-z0-9_\\-]+","_", base)
      safe_pc <- gsub("[^A-Za-z0-9_\\-]+","_", pc)
      
      for (f in factor_numbers) {
        dsub <- dat[dat$Factor == f, , drop = FALSE]
        
        # ---- R² barplot ----
        rsq_tbl <- do.call(rbind, lapply(split(dsub, dsub$View), function(dv) {
          if (nrow(dv) < 3) return(NULL)
          fit <- lm(neglog10P ~ Weight, data = dv)
          r2 <- summary(fit)$r.squared
          data.frame(View = unique(dv$View), R2 = r2, stringsAsFactors = FALSE)
        }))
        
        if (!is.null(rsq_tbl) && nrow(rsq_tbl) > 0) {
          p_r2 <- ggplot2::ggplot(rsq_tbl, ggplot2::aes(x = View, y = R2, fill = View)) +
            ggplot2::geom_col() +
            ggplot2::scale_fill_manual(values = view.colors) +
            ggplot2::theme_classic() +
            ggplot2::ylim(0, 1) +
            ggplot2::labs(
              title = paste("Factor", f, "R² by View"),
              x = "View", y = expression(R^2)
            ) +
            ggplot2::theme(
              plot.title = ggplot2::element_text(face = "bold", size = 12),
              axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
            ) 
        } else {
          p_r2 <- ggplot2::ggplot() + ggplot2::theme_void() +
            ggplot2::labs(title = paste("Factor", f, "— no R² data"))
        }
        
        # ---- Scatter ----
        if (nrow(dsub) >= 3) {
          p_scatter <- ggplot2::ggplot(dsub, ggplot2::aes(x = Weight, y = neglog10P, colour = View)) +
            ggplot2::geom_point(alpha = 0.5, size = 1.3) +
            ggplot2::geom_smooth(method = "lm", se = FALSE, size = 0.6) +
            ggplot2::scale_color_manual(values = view.colors) +
            ggplot2::theme_classic() +
            ggplot2::labs(
              title = paste("Factor", f, "Weight vs –log10(P)"),
              x = "Weight", y = expression(-log[10](P)), colour = "View"
            ) +
            ggplot2::theme(
              plot.title = ggplot2::element_text(face = "bold", size = 12)
            )
        } else {
          p_scatter <- ggplot2::ggplot() + ggplot2::theme_void() +
            ggplot2::labs(title = paste("Factor", f, "— no scatter data"))
        }
        
        # Store both plots in correct order: col1=bar, col2=scatter
        plot_list <- append(plot_list, list(p_r2, p_scatter))
        
        # ---- Save individual plots ----
        out_r2 <- file.path(output_dir_r2,
                            paste0("R2__", safe_base, "_", safe_pc, "_Factor", f, ".png"))
        out_sc <- file.path(output_dir_scatter,
                            paste0("Scatter__", safe_base, "_", safe_pc, "_Factor", f, ".png"))
        
        ggplot2::ggsave(out_r2, plot = p_r2, width = 5, height = 4, units = "in", dpi = 150)
        ggplot2::ggsave(out_sc, plot = p_scatter, width = 6, height = 4, units = "in", dpi = 150)
      }
      
      # ---- Arrange grid: 2 columns, one row per factor ----
      final_plot <- gridExtra::arrangeGrob(
        grobs = plot_list,
        ncol = 2,
        nrow = length(factor_numbers),
        top = grid::textGrob(
          paste0("R² and Scatter by Factor — ", pc),
          gp = grid::gpar(fontsize = 14, fontface = "bold")
        )
      )
      
      # ---- Naming ----
      outfile <- file.path(outdir, paste0("R2_scatter__", safe_base, ".png"))
      
      grDevices::png(outfile, width = 16, height = 4 * length(factor_numbers), units = "in", res = 300)
      grid::grid.draw(final_plot)
      grDevices::dev.off()
      
      message("Saved: ", outfile)
    }
  }
}

########################################
# VIOLIN PLOTS
########################################

plot_violin_pval_vs_weight <- function(view.colors,
                                       indir = "Results/Significance_testing/tables",
                                       outdir = "Results/Significance_testing/violin_plots") {
  
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  files <- list.files(indir, pattern = "\\.csv$", full.names = TRUE)
  if (length(files) == 0L) stop("No CSV files found in ", indir)
  
  for (file in files) {
    message("Processing violin for: ", basename(file))
    df <- utils::read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
    
    if (!all(c("Feature","View") %in% names(df))) next
    
    # Hard-coded params
    p_cols <- grep("Pval$", names(df), value = TRUE)
    w_prefix <- "^Weight_"
    weight_cols <- grep(w_prefix, names(df), value = TRUE)
    
    if (length(p_cols) == 0L || length(weight_cols) == 0L) next
    for (pcol in p_cols) df[[pcol]] <- suppressWarnings(as.numeric(df[[pcol]]))
    
    # Melt into long form
    long <- reshape2::melt(
      df,
      id.vars = c("Feature","View", p_cols),
      measure.vars = weight_cols,
      variable.name = "FactorCol",
      value.name = "Weight"
    )
    long$Factor <- sub(w_prefix, "", long$FactorCol)
    long$Weight <- abs(long$Weight)
    
    for (pcol in p_cols) {
      dat <- long[, c("Feature","View","Factor","Weight", pcol)]
      names(dat)[5] <- "Pval"
      dat <- dat[!is.na(dat$Pval), , drop = FALSE]
      if (nrow(dat) == 0L) next
      
      dat$Pval <- pmax(pmin(dat$Pval, 1), .Machine$double.xmin)
      
      # More bins
      dat$Pbin <- cut(
        dat$Pval,
        breaks = c(0, 1e-5, 0.01, 1),
        include.lowest = TRUE, right = TRUE,
        labels = c("≤1e-3","(1e-3,0.01]","(0.01,1]")
      )
      
      # Violin plot faceted in grid: rows=Factor, cols=View
      p_violin <- ggplot2::ggplot(dat, ggplot2::aes(x = Pbin, y = Weight, fill = View)) +
        ggplot2::geom_violin(trim = FALSE, alpha = 0.7) +
        ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.5) +
        ggplot2::scale_fill_manual(values = view.colors) +
        ggplot2::facet_grid(Factor ~ View, scales = "free_y") +
        ggplot2::theme_classic() +
        ggplot2::labs(
          title = paste0("Weight distribution by P-value bin — ", pcol),
          x = "P-value bin", y = "Weight"
        ) +
        ggplot2::theme(
          legend.position = "none",
          strip.text = ggplot2::element_text(face = "bold")
        )
      
      # ---- File naming (only base filename, strip prefix) ----
      base <- tools::file_path_sans_ext(basename(file))
      base <- sub("^Weights_AllViews_AllFactors__", "", base)
      outfile <- file.path(outdir, paste0("Violin__", base, ".png"))
      
      # dynamic size
      nfac <- length(unique(dat$Factor))
      nview <- length(unique(dat$View))
      ggplot2::ggsave(outfile, plot = p_violin,
                      width = 3 + nview * 2, height = 2 + nfac * 1.5, dpi = 300)
      message("Saved: ", outfile)
    }
  }
}

########################################
# TOP LEVEL HEATMAP
########################################
plot_top_weights_heatmap_faceted_by_view <- function(MOFAobject.trained,
                                                     out_file = "./Results/Per_Factor_Analysis/TopWeights_Heatmap_FacetedByView.png",
                                                     top_n = 5,
                                                     factors = NULL,
                                                     use_abs = FALSE,
                                                     strip_view_suffix = TRUE,
                                                     prefix_view = TRUE,
                                                     gap_mm = 2) {
  stopifnot(!is.null(MOFAobject.trained))
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  
  if (!requireNamespace("MOFA2", quietly = TRUE)) stop("Package 'MOFA2' is required.")
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) stop("Package 'ComplexHeatmap' is required.")
  if (!requireNamespace("circlize", quietly = TRUE)) stop("Package 'circlize' is required.")
  if (!requireNamespace("grid", quietly = TRUE)) stop("Package 'grid' is required.")
  
  W_list <- MOFA2::get_weights(MOFAobject.trained)
  if (!is.list(W_list) || length(W_list) == 0) stop("No weights found in MOFA object.")
  
  # Factors to include
  if (is.null(factors)) {
    factors <- colnames(W_list[[1]])
  } else {
    factors <- as.character(factors)
  }
  if (length(factors) == 0) stop("No factors specified/found.")
  
  mats <- list()
  row_view <- character()
  
  for (view_name in names(W_list)) {
    W <- W_list[[view_name]]
    facs <- intersect(factors, colnames(W))
    if (length(facs) == 0) next
    
    feats_raw <- rownames(W)
    feats <- if (strip_view_suffix) gsub(paste0("_", view_name), "", feats_raw) else feats_raw
    
    # union of top features across factors (within view)
    top_feats <- unique(unlist(lapply(facs, function(f) {
      ord <- order(abs(W[, f]), decreasing = TRUE)
      feats[ord][seq_len(min(top_n, length(ord)))]
    })))
    
    mat <- W[match(top_feats, feats), facs, drop = FALSE]
    rownames(mat) <- top_feats
    if (use_abs) mat <- abs(mat)
    
    # avoid name collisions across views (internal rownames)
    if (prefix_view) rownames(mat) <- paste0(view_name, "::", rownames(mat))
    
    mats[[view_name]] <- mat
    row_view <- c(row_view, rep(view_name, nrow(mat)))
  }
  
  if (length(mats) == 0) stop("No data assembled for plotting.")
  
  # stack views (rows) into one matrix
  mat_all <- do.call(rbind, mats)
  row_split <- factor(row_view, levels = names(mats))
  
  # display labels: strip the view prefix for plotting only
  labels_row <- rownames(mat_all)
  if (prefix_view) labels_row <- sub("^[^:]+::", "", labels_row)
  
  # signed weights colour mapping centred at 0
  w_min <- min(mat_all, na.rm = TRUE)
  w_max <- max(mat_all, na.rm = TRUE)
  col_fun <- circlize::colorRamp2(c(w_min, 0, w_max), c("blue", "yellow", "red"))
  
  # row annotation (view) - optional; keep if you want a colour strip
  ha <- ComplexHeatmap::rowAnnotation(
    View = row_split,
    show_annotation_name = FALSE
  )
  
  ht <- ComplexHeatmap::Heatmap(
    mat_all,
    name = "Weight",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = row_split,
    row_gap = grid::unit(gap_mm, "mm"),
    left_annotation = ha,
    
    show_row_names = TRUE,
    row_labels = labels_row,
    # <- clean row labels shown on plot
    row_names_gp = grid::gpar(fontsize = 5), # <- smaller feature labels
    
    row_title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
    row_title_rot = 0,
    row_title_side = "left",
    
    column_names_rot = 90,
    column_names_gp = grid::gpar(fontsize = 9),
    
    heatmap_legend_param = list(title = "Weight")
  )
  
  # sensible device size: scale height by total rows
  n_rows <- nrow(mat_all)
  n_cols <- ncol(mat_all)
  
  grDevices::png(
    filename = out_file,
    width  = max(1200, 120 * n_cols),
    height = max(900, 12 * n_rows),
    res = 150
  )
  ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  grDevices::dev.off()
  
  invisible(list(matrix = mat_all, heatmap = ht))
}

########################################
# TOP WEIGHTS BAR PLOTS
########################################

plot_top_weights_view <- function(factor_name) {
  views <- names(MOFAobject.trained@expectations$W)
  for (view_name in views) {
    message(
      paste(
        "Generating top weights plot for view:",
        view_name,
        "and factor:",
        factor_name
      )
    )
    output_dir <- ensure_output_dirs(factor_name, view_name)
    
    top_weights_plot <- plot_top_weights(
      MOFAobject.trained.cleanNames,
      view = view_name,
      factors = factor_name,
      nfeatures = 25,
      abs = TRUE,
      scale = TRUE,
      sign = "all"
    ) +
      theme(
        plot.title = element_text(face = "bold", size = 14, colour = "black"),
        axis.text = element_text(size = 10, face = "bold", colour = "black"),
        axis.title = element_text(size = 12, face = "bold", colour = "black"),
        legend.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 12, colour = "black"),
        panel.grid = element_blank()
      )
    ggsave(
      plot = top_weights_plot,
      filename = file.path(output_dir, "Top_Weights.png"),
      width = 15,
      height = 15
    )
  }
}

########################################
# SCATTER PLOT PER VIEW
########################################

plot_scatter_per_view <- function(factor_name) {
  views <- names(MOFAobject.trained@expectations$W)
  
  for (view_name in views) {
    for (meta_col in names(group.colors)) {
      message(sprintf(
        "Generating scatter plot for view: %s, factor: %s, grouping: %s",
        view_name, factor_name, meta_col
      ))
      
      if(nrow(as.data.frame(MOFAobject.trained@expectations$W[view_name])) >= 2){
        
        output_dir <- ensure_output_dirs(factor_name, paste(view_name,"/scatterplots",sep=""))
        
        scatter_plot <- plot_data_scatter(
          MOFAobject.trained,
          view = view_name,
          factor = factor_name,
          features = 20,
          dot_size = 3,
          color_by = meta_col,
          legend = TRUE
        ) +
          scale_fill_manual(values = group.colors[[meta_col]]) +
          scale_color_manual(values = group.colors[[meta_col]]) +
          theme(
            plot.title = element_text(face = "bold", size = 14, colour = "black"),
            axis.text = element_text(size = 8, colour = "black"),
            axis.title = element_text(size = 10, face = "bold", colour = "black"),
            legend.text = element_text(size = 8, colour = "black"),
            legend.title = element_text(size = 10, colour = "black"),
            panel.grid = element_blank()
          )
        
        # Output file: prefix with grouping variable name
        out_file <- file.path(output_dir,
                              paste0(meta_col, "_Scatter.png"))
        
        ggsave(
          plot = scatter_plot,
          filename = out_file,
          width = 10,
          height = 5,
          units = "in"
        )
      } else {
        print("Not enough features to plot a scatter")
      }
    }
  }
}

########################################
# TOP WEIGHTS PER FACTOR
########################################

plot_top_weights_per_factor <- function(factor_name, top_n = 10) {
  message(paste("Generating top weights plot for factor:", factor_name))
  
  # Collect top features across all views
  top_features_list <- list()
  
  for (view_name in names(get_weights(MOFAobject.trained))) {
    weight_df <- data.frame(
      Feature = gsub(paste0("_", view_name), "", 
                     rownames(MOFAobject.trained.cleanNames@expectations$W[[view_name]])),
      View = view_name,
      Weight = MOFAobject.trained@expectations$W[[view_name]][, factor_name]
    )
    
    # Sort and select top features
    weight_df <- weight_df %>% arrange(desc(abs(Weight)))
    weight_df <- weight_df[0 < weight_df$Weight | 0 > weight_df$Weight,]
    top_features_list[[view_name]] <- head(weight_df, top_n)
  }
  
  # Combine and retain top features across all views
  top_features_df <- bind_rows(top_features_list)
  top_weights_df <- top_features_df %>%
    arrange(desc(abs(Weight))) 
  
  # Construct plot
  top_weights_plot <- ggplot(top_weights_df,
                             aes(x = reorder(Feature, -abs(Weight)),
                                 y = Weight,
                                 fill = View)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top Weights for", factor_name),
         x = "Feature",
         y = "Weight") +
    theme_classic() +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      axis.text.y = element_text(size = 7)
    )
  
  ggsave(
    plot = top_weights_plot,
    filename = file.path("./Results/Per_Factor_Analysis",factor_name, "Top-weights-bar.png"),
    width = 15,
    height = length(group.colors)/1.5 + 10
  )
}

########################################
# INDIVIDUAL HEATMAP 
########################################

plot_heatmap_per_view <- function(factor_name, top_features = 25) {
  views <- names(MOFAobject.trained.cleanNames@expectations$W)
  
  for (view_name in views) {
    message(sprintf("Generating heatmap for view: %s | factor: %s | all covariates",
                    view_name, factor_name))
    
    # Retrieve the underlying data matrix for this view
    data_view <- as.data.frame(MOFAobject.trained.cleanNames@data[[view_name]]$group1)
    
    # --- Skip empty / all-NA / zero-variance views ---
    if (all(is.na(data_view))) {
      message(sprintf("  Skipping view %s (all values NA)", view_name))
      next
    }
    # remove rows with all NAs
    data_view <- data_view[rowSums(!is.na(data_view)) > 0, , drop = FALSE]
    if (nrow(data_view) < 2 || ncol(data_view) < 2) {
      message(sprintf("  Skipping view %s (insufficient data for clustering)", view_name))
      next
    }
    # check variance
    if (all(apply(data_view, 1, sd, na.rm = TRUE) == 0)) {
      message(sprintf("  Skipping view %s (zero variance across features)", view_name))
      next
    }
    
    # Proceed with heatmap generation
    output_dir <- ensure_output_dirs(factor_name, view_name)
    
    heatmap_plot <- tryCatch(
      plot_data_heatmap(
        MOFAobject.trained.cleanNames,
        factor = factor_name,
        view = view_name,
        denoise = TRUE,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        show_colnames = TRUE,
        show_rownames = TRUE,
        annotation_samples = names(group.colors)[1],
        features = top_features,
        annotation_colors = group.colors[1],
        annotation_legend = TRUE,
        scale = "row"
      ),
      error = function(e) {
        message(sprintf("  Skipping view %s (heatmap failed: %s)", view_name, e$message))
        return(NULL)
      }
    )
    
    if (is.null(heatmap_plot)) next
    
    ggsave(
      plot = heatmap_plot,
      filename = file.path(output_dir, "Heatmap.png"),
      width = length(metadata$category)/4 + 5,
      height = length(group.colors)/4 + 5
    )
    
    tryCatch(while (dev.cur() > 1) dev.off(), error = function(e) {})
  }
}

########################################
# WEIGHT TABLES
########################################

create_weight_outputs <- function(factor_name,view_metadata) {
  # Ensure output directory exists
  output_dir <- ensure_output_dirs(factor_name, "")
  
  view_groups <- unlist(as.list(na.omit(unique(view_metadata$View_Group))))
  
  # Collect top features and weights across all views for this factor
  top_features_list <- list()
  all_weights_list <- list()
  for (view_name in names(get_weights(MOFAobject.trained))) {
    message(
      paste(
        "Generating top features lists for view:",
        view_name,
        "and factor:",
        factor_name
      )
    )
    weight_df <- data.frame(
      Feature = gsub(
        paste('_', view_name, sep = ""),
        "",
        rownames(MOFAobject.trained.cleanNames@expectations$W[[view_name]])
      ),
      View = view_name,
      Weight = MOFAobject.trained@expectations$W[[view_name]][, factor_name],
      AdjWeight = MOFAobject.trained@expectations$adjusted_W[[view_name]][, factor_name]
    )
    # Sort by absolute weight and select top 10 features
    weight_df <- weight_df %>% arrange(desc(abs(Weight)))
    top_features_list[[view_name]] <- head(weight_df, 10)
    all_weights_list[[view_name]] <- weight_df
  }
  # Merge and save all weights to a CSV file
  all_weights_df <- bind_rows(all_weights_list)
  write.csv(all_weights_df,
            paste0(output_dir, "/All_Weights.csv"),
            row.names = FALSE)
  
  # For each defined view group, combine top features and generate a bar plot
  for (group_name in names(view_groups)) {
    view_set <- view_groups[[group_name]]
    # Combine top features from all views in this group
    filtered_top_features <- top_features_list[view_set]
    top_features_df <- bind_rows(filtered_top_features)
    # Save top features table for this group
    write.csv(
      top_features_df,
      file = paste0(output_dir, "/", group_name, "-Top_Features.csv"),
      row.names = FALSE
    )
    # Create bar plot of top feature weights for this group
    top_weights_df <- top_features_df %>% arrange(desc(abs(Weight)))
    top_weights_df$Feature <- substr(top_weights_df$Feature, 1, 40)
    top_weights_plot <- ggplot(top_weights_df, aes(
      x = reorder(Feature, -abs(Weight)),
      y = Weight,
      fill = View
    )) +
      geom_bar(stat = "identity") +
      coord_flip() +
      labs(
        title = paste("Top Weights for", factor_name, "-", group_name),
        x = "Feature",
        y = "Weight"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(face = "bold", size = 14, colour = "black"),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text = element_text(size = 8, colour = "black"),
        axis.title = element_text(size = 10, face = "bold", colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        legend.title = element_text(size = 10, colour = "black"),
        panel.grid = element_blank()
      )
    ggsave(
      plot = top_weights_plot,
      filename = paste0(output_dir, "/", group_name, "-Top_Weights.png"),
      width = 10,
      height = 8
    )
  }
  return(invisible(NULL))
}

########################################
# Helper function to get GSEA input vectors for a given factor
########################################

get_gsea_input <- function(factor_name) {
  # Extract KEGG weights for the current factor across views
  kegg_weights <- MOFAobject.trained@expectations$KEGG_weights
  # Filter out any with less than 2 entries (if applicable)
  kegg_weights_filtered <- Filter(function(x)
    nrow(x) > 1, kegg_weights)
  
  # Extract named weight vectors for the specified factor from each available matrix/data frame
  extracted_list <- lapply(kegg_weights_filtered, function(df) {
    if ((is.matrix(df) ||
         is.data.frame(df)) && factor_name %in% colnames(df)) {
      vec <- df[, factor_name]
      names(vec) <- rownames(df)
      return(vec)
    } else {
      return(NULL)
    }
  })
  # Remove any NULL entries
  extracted_list <- extracted_list[!sapply(extracted_list, is.null)]
  
  # Deduplicate weights by taking mean for duplicated names, and sort in decreasing order
  dedup_vector <- function(x) {
    v <- tapply(x, names(x), mean)
    v <- as.vector(v)
    names(v) <- names(x <- tapply(x, names(x), mean))
    v[order(v, decreasing = TRUE)]
  }
  # Apply deduplication to all extracted vectors
  gsea_input <- lapply(extracted_list, dedup_vector)
  return(gsea_input)
}

########################################
# Perform KEGG enrichment analysis for a given factor
########################################

perform_kegg_enrichment <- function(factor_name,
                                    MOFAobject.trained,
                                    view_metadata,
                                    TERM2GENE_list,
                                    TERM2NAME_list,
                                    outdir = "Results/Per_Factor_Analysis") {
  
  output_dir <- file.path(outdir, factor_name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  message("Generating FSEA for factor: ", factor_name)
  
  # ---- Build GSEA input ----
  kegg_weights <- MOFAobject.trained@expectations$KEGG_weights
  kegg_weights_filtered <- Filter(function(x) nrow(x) > 1, kegg_weights)
  
  get_gsea_input <- function(kegg_weights_filtered, factor_name) {
    
    # helper: deduplicate by name using mean, keep numeric, sort decreasing
    dedup_vector <- function(x) {
      keep <- !is.na(names(x)) & names(x) != ""
      x <- x[keep]
      
      v <- tapply(x, names(x), mean, na.rm = TRUE)
      v <- sort(v, decreasing = TRUE)
      
      as.numeric(v) |> `names<-`(names(v))
    }
    
    extracted_list <- lapply(kegg_weights_filtered, function(df) {
      
      if (!is.matrix(df) && !is.data.frame(df)) {
        return(NULL)
      }
      
      if (!(factor_name %in% colnames(df))) {
        return(NULL)
      }
      
      vec <- df[, factor_name, drop = TRUE]
      names(vec) <- rownames(df)
      
      dedup_vector(vec)
    })
    
    extracted_list[!vapply(extracted_list, is.null, logical(1))]
  }
  
  gsea_input <- get_gsea_input(
    kegg_weights_filtered = kegg_weights_filtered,
    factor_name = factor_name
  )
  
  # ---- Loop once over all inputs ----
  res_list <- list()
  
  for (nm in names(gsea_input)) {
    
    print(nm)
    # split name: e.g. "metabolomics_hvg"
    parts <- strsplit(nm, "_")[[1]]
    view <- paste(parts[-length(parts)], collapse = "_")
    sp   <- parts[length(parts)]
    
    if (view %in% subset(view_metadata, Type == "Metabolic" & Active == "Y")$View) {
      TERM2GENE <- TERM2GENE_list[[sp]]$metabolite
      TERM2NAME <- TERM2NAME_list[[sp]]$metabolite
    } else if (view %in% subset(view_metadata, Type == "RNAseq" & Active == "Y")$View) {
      TERM2GENE <- TERM2GENE_list[[sp]]$gene
      TERM2NAME <- TERM2NAME_list[[sp]]$gene
    } else {
      next
    }
    
    ids_gsea  <- names(gsea_input[[nm]])
    ids_sets <- unique(TERM2GENE[, 2])
    
    overlap <- intersect(ids_gsea, ids_sets)
    
    if (length(overlap) < 10) {
      message(
        sprintf(
          "Skipping %s: only %d overlapping IDs",
          nm, length(overlap)
        )
      )
      next
    }
    
    gsea_res <- clusterProfiler::GSEA(
      geneList   = gsea_input[[nm]],
      TERM2GENE  = TERM2GENE,
      TERM2NAME  = TERM2NAME,
      pvalueCutoff = 0.05,
      pAdjustMethod = "none",
      minGSSize = 10,
      maxGSSize = 500,
      nPermSimple = 10000,
      eps = 0
    )
    
    if (!is.null(gsea_res) && nrow(gsea_res@result) > 0) {
      df <- gsea_res@result
      df$Cluster <- nm
      res_list[[nm]] <- df
    }
  }
  
  # ---- Stitch back into compareClusterResult ----
  if (length(res_list) > 0) {
    cc_df <- do.call(rbind, res_list)
    
    # Move Cluster to the first column
    if ("Cluster" %in% colnames(cc_df)) {
      cc_df <- cc_df[, c("Cluster", setdiff(colnames(cc_df), "Cluster"))]
    }
    
    combined_result <- new("compareClusterResult",
                           compareClusterResult = cc_df,
                           geneClusters = gsea_input)
    
    output_file <- file.path(output_dir, "KEGG_Enrichment_Weighted.csv")
    utils::write.csv(combined_result, output_file, row.names = FALSE)
    message("Saved enrichment table: ", output_file)
    
    # Plot
    KEGG_Enrichment_plot <- enrichplot::dotplot(
      combined_result,
      color = "pvalue",
      font.size = 8,
      showCategory = 5,
      title = paste0(factor_name, ": FSEA"),
      split = ".sign"
    ) +
      facet_grid(. ~ .sign) +
      theme(
        plot.title.position = "panel",              # key: keep title over panels, not legend
        plot.title = element_text(
          face = "bold", size = 14, colour = "black",
          margin = margin(b = 6)
        ),
        legend.key.height = unit(0.8, "cm"),
        legend.margin = margin(t = 10, r = 5, b = 5, l = 5),
        axis.text = element_text(size = 8, colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title = element_text(size = 10, face = "bold", colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        legend.title = element_text(size = 10, colour = "black"),
        panel.grid = element_blank()
      )
    
    outfile <- file.path(output_dir, "KEGG_Enrichment.png")
    ggsave(plot = KEGG_Enrichment_plot, 
           filename = outfile,
           width = length(unique(combined_result@compareClusterResult$Cluster))/2 + 5, 
           height = length(unique(combined_result@compareClusterResult$ID))/4 + 3
           )
    
    message("Saved enrichment plot: ", outfile)
  } else {
    message("No significant GSEA results available for factor: ", factor_name)
  }
  
  return(invisible(NULL))
}

########################################
# PLOT ENRICHMENT NETWORK
########################################

plot_enrichment_network <- function(factor_name,
                                    gene_annotation,
                                    metabolite_annotation,
                                    view.colors) {
  # Ensure output directory exists
  output_dir <- ensure_output_dirs(factor_name, "")
  
  message(paste("Generating Igraph cnplot for factor:", factor_name))
  
  # ---- Read enrichment results ----
  enrichment_file <- file.path(output_dir, "KEGG_Enrichment_Weighted.csv")
  if (!file.exists(enrichment_file)) {
    message("No enrichment results file found for factor ", factor_name)
    return(invisible(NULL))
  }
  enrichment_df <- read.csv(enrichment_file, stringsAsFactors = FALSE)
  
  # ---- KEGG ID → name mapping ----
  kegg_to_name <- c()
  if (is.list(metabolite_annotation)) {
    for (sp in names(metabolite_annotation)) {
      df <- metabolite_annotation[[sp]]
      if (is.data.frame(df) && nrow(df) > 0 && "KEGG_ID" %in% colnames(df)) {
        map <- df[!is.na(df$KEGG_ID), ]
        kegg_to_name <- c(
          kegg_to_name,
          setNames(map$Metabolite, gsub("^cpd:", "", map$KEGG_ID))
        )
      }
    }
  }
  if (is.list(gene_annotation)) {
    for (sp in names(gene_annotation)) {
      df <- gene_annotation[[sp]]
      if (is.data.frame(df) && nrow(df) > 0 && "KEGG_ID" %in% colnames(df)) {
        map <- df[!is.na(df$KEGG_ID), ]
        kegg_to_name <- c(
          kegg_to_name,
          setNames(map$gene, gsub("^\\w+:", "", map$KEGG_ID))
        )
      }
    }
  }
  
  # Replace KEGG IDs with readable names
  replace_kegg_with_name <- function(KEGG_IDs, mapping) {
    ids <- unlist(strsplit(KEGG_IDs, "/", fixed = TRUE))
    names_vec <- mapping[ids]
    names_vec[is.na(names_vec)] <- ids[is.na(names_vec)]
    paste(names_vec, collapse = "/")
  }
  annotated_enrichment_df <- enrichment_df
  annotated_enrichment_df$core_enrichment <- sapply(
    enrichment_df$core_enrichment,
    replace_kegg_with_name,
    mapping = kegg_to_name
  )
  
  # ---- Build edge list ----
  edge_list <- list()
  for (i in seq_len(nrow(annotated_enrichment_df))) {
    cluster <- annotated_enrichment_df$Cluster[i]
    term <- annotated_enrichment_df$Description[i]
    compounds_str <- annotated_enrichment_df$core_enrichment[i]
    pvalue <- annotated_enrichment_df$pvalue[i]
    enrichmentScore <- annotated_enrichment_df$enrichmentScore[i]
    if (is.na(compounds_str) || compounds_str == "") next
    compounds <- strsplit(compounds_str, "/", fixed = TRUE)[[1]]
    edge_list[[i]] <- data.frame(
      Cluster = rep(cluster, length(compounds)),
      term = rep(term, length(compounds)),
      pvalue = rep(pvalue, length(compounds)),
      enrichmentScore = rep(enrichmentScore, length(compounds)),
      compound = compounds,
      stringsAsFactors = FALSE
    )
  }
  if (length(edge_list) == 0) {
    message("No enriched terms to visualize for factor ", factor_name)
    return(invisible(NULL))
  }
  edges_all <- do.call(rbind, edge_list)
  edges_all <- edges_all[edges_all$pvalue < 0.05, ]
  
  # ---- Normalize cluster names (strip _species for coloring) ----
  edges_all <- edges_all %>%
    mutate(
      Cluster_clean = sub("_[^_]+$", "", Cluster),
      compound_label = compound,
      compound_clustered = paste0(compound, " (", Cluster, ")")
    )
  
  # ---- Select top 100 edges ----
  edges <- edges_all %>%
    dplyr::select(from = term,
                  to = compound_clustered,
                  enrichmentScore,
                  Cluster = Cluster_clean,
                  pvalue) %>%
    dplyr::arrange(desc(abs(enrichmentScore))) %>%
    head(100)
  
  if (nrow(edges) == 0) {
    message("No significant enriched connections to plot for factor ", factor_name)
    return(invisible(NULL))
  }
  
  # ---- Build igraph ----
  g <- graph_from_data_frame(edges, directed = FALSE)
  terms <- unique(edges$from)
  V(g)$type <- ifelse(V(g)$name %in% terms, "term", "compound")
  
  # Aggregate compound attributes
  compound_attrs <- edges_all %>%
    dplyr::select(node = compound_clustered,
                  label = compound_label,
                  Cluster = Cluster_clean,
                  pvalue) %>%
    dplyr::group_by(node) %>%
    dplyr::summarise(
      pvalue = mean(pvalue, na.rm = TRUE),
      Cluster = dplyr::first(as.character(Cluster)),
      label = dplyr::first(label),
      .groups = "drop"
    )
  
  # Initialize attributes
  V(g)$pvalue <- NA_real_
  V(g)$Cluster <- NA_character_
  V(g)$label <- V(g)$name
  
  # Assign attributes
  compound_idx <- match(V(g)$name, compound_attrs$node)
  valid_idx <- !is.na(compound_idx)
  V(g)$pvalue[valid_idx] <- compound_attrs$pvalue[compound_idx[valid_idx]]
  V(g)$Cluster[valid_idx] <- compound_attrs$Cluster[compound_idx[valid_idx]]
  V(g)$label[valid_idx] <- compound_attrs$label[compound_idx[valid_idx]]
  V(g)$label <- gsub(" - ", "\n", V(g)$label)
  
  # ---- Node sizes and colors ----
  V(g)$size <- -log10(V(g)$pvalue) * 3
  V(g)$size[is.na(V(g)$size)] <- 10
  V(g)$size[V(g)$size > 10] <- 10
  V(g)$color <- ifelse(
    V(g)$type == "term",
    "grey95",
    view.colors[V(g)$Cluster]
  )
  V(g)$color[is.na(V(g)$color)] <- "grey96"
  
  V(g)$label.cex <- ifelse(V(g)$type == "term", 0.4, 0.4)
  V(g)$label.color <- ifelse(V(g)$type == "term", "#006400", "black")
  
  # Edge width scaled
  E(g)$width <- rescale(abs(E(g)$enrichmentScore), to = c(1, 5))
  
  # ---- Plot ----
  png(
    filename = file.path(output_dir, "Top_100_edges_KEGG_cnet.png"),
    width = 1500, height = 1500, res = 300
  )
  par(family = "Arial", mar = c(2, 2, 5, 2))
  set.seed(42)
  coords <- layout_with_fr(g)
  
  plot(
    g,
    layout = coords,
    vertex.label = V(g)$label,
    vertex.label.cex = V(g)$label.cex,
    vertex.label.color = V(g)$label.color,
    vertex.label.family = "Arial",
    vertex.label.dist = 0.5,
    edge.color = "grey80"
  )
  title(main = "CnetPlot of Top 100 Significantly Enriched (p < 0.05) Features", cex.main = 1)
  
  # ---- Legends ----
  legend("bottomleft", legend = names(view.colors), col = view.colors,
         pch = 21, pt.bg = view.colors, pt.cex = 1, cex = 0.5, bty = "n",
         title = "Dataset View")
  
  legend_pvals <- c(0.05, 0.01, 0.001)
  legend_sizes <- -log10(legend_pvals) * 3
  pt_cex_legend <- legend_sizes / 4
  legend("topright",
         legend = c("    p = 0.05", "", "    p = 0.01", "", "    p = 0.001"),
         pt.cex = c(pt_cex_legend[1], NA, pt_cex_legend[2], NA, pt_cex_legend[3]),
         pch = 21, pt.bg = "white", col = "black", cex = 0.5, bty = "n",
         title = "Enrichment Pval")
  
  legend("topleft",
         legend = c("Low enrichment", "High enrichment"),
         lwd = c(1, 5), col = "grey80", cex = 0.5, bty = "n",
         title = "Enrichment Score")
  
  tryCatch(while (dev.cur() > 1) dev.off(), error = function(e) {})
  return(invisible(NULL))
}


########################################
# GENERATE PATHVIEW 
########################################

generate_simple_pathway_graph <- function(factor_name, view_metadata, MOFAobject.trained,
                                          TERM2GENE_list = NULL) {
  enrichment_file <- file.path("Results/Per_Factor_Analysis", factor_name,
                               "KEGG_Enrichment_Weighted.csv")
  if (!file.exists(enrichment_file)) {
    message("No enrichment results for ", factor_name)
    return(invisible(NULL))
  }
  
  enrichment_df <- read.csv(enrichment_file, stringsAsFactors = FALSE)
  enrichment_df <- enrichment_df[enrichment_df$pvalue < 0.05, ]
  if (nrow(enrichment_df) == 0) {
    message("No significant pathways for ", factor_name)
    return(invisible(NULL))
  }
  
  weights <- MOFAobject.trained@expectations$KEGG_weights
  
  for (i in seq_len(nrow(enrichment_df))) {
    cluster <- enrichment_df$Cluster[i]   # e.g. "LCMS_Aphid_api"
    parts   <- strsplit(cluster, "_")[[1]]
    view    <- paste(parts[-length(parts)], collapse = "_")
    sp      <- parts[length(parts)]
    clean_id <- gsub(sp, "", enrichment_df$ID[i])
    raw_id   <- enrichment_df$ID[i]
    
    if (!cluster %in% names(weights)) {
      message("No weights for cluster ", cluster, " → skipping")
      next
    }
    
    mat <- weights[[cluster]]
    if (is.null(mat) || !(factor_name %in% colnames(mat))) {
      message("No factor weights for ", cluster, " in ", factor_name)
      next
    }
    
    vec <- mat[, factor_name]
    names(vec) <- rownames(mat)
    
    # Output dir
    output_dir <- file.path("Results/Per_Factor_Analysis", factor_name, view, "KEGG_pathway_plots")
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    message("Running Pathview for pathway: ", sp, clean_id, " (cluster ", cluster, ")")
    
    pv <- tryCatch({
      # 1. Load pathway graph
      
      xml_dir <- file.path("Annotations/kegg_annotation_cache", sp)
      if (!dir.exists(xml_dir)) dir.create(xml_dir, recursive = TRUE)
      graph <- ggkegg::pathway(raw_id, use_cache = FALSE,directory = xml_dir)
      
      # 2. Attach weights
      graph <- graph %>%
        activate(nodes) %>%
        mutate(weight = vec[graphics_name],
               degree = centrality_degree(mode = "all"))
      
      # 3. Extract pathway title (fallback to raw_id if missing)
      title_row <- as_tibble(graph, active = "nodes") %>%
        filter(type == "map", grepl("^TITLE:", graphics_name))
      
      if (nrow(title_row) == 0) {
        pathway_title <- raw_id
      } else {
        pathway_title <- sub("^TITLE:", "", title_row$graphics_name)
        
        # Replace any character that is not A–Z, a–z, 0–9, or underscore
        pathway_title <- gsub("[^A-Za-z0-9_]", "_", pathway_title)
        
        # Collapse multiple underscores to a single one
        pathway_title <- gsub("_+", "_", pathway_title)
        
        # Trim underscores from start/end
        pathway_title <- gsub("^_|_$", "", pathway_title)
      }
      
      # 4. Plot with KEGG background + weights overlay
      p <- graph |>
        ggraph(layout = "manual", x = x, y = y) +
        geom_node_rect(aes(fill = weight,
                           filter = type %in% c("compound", "gene")),
                       colour = "black", size = 0.2, na.rm = TRUE) +
        overlay_raw_map() +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                             midpoint = 0, name = "Weight") +
        theme_void()
      
      # 5. Save output
      output_path <- file.path(output_dir, paste0(pathway_title, ".png"))
      ggsave(output_path, plot = p, width = 10, height = 8, dpi = 300)
      message("Saved ggkegg pathway plot: ", output_path)
      
    }, error = function(e) {
      message("ggkegg failed for ", raw_id, " (", sp, "): ", conditionMessage(e))
    })
    
    
  }
}

########################################
# RUN FELLA
########################################

run_fella_enrichment <- function(factor_name, view_metadata, fella_list, TERM2GENE_list, MOFAobject.trained) {
  
  # Derive metabolomic sheets from view_metadata
  metabolomic_views <- subset(view_metadata, Active == "Y" & Type == "Metabolic", select = View)$View
  if (length(metabolomic_views) == 0) {
    warning("No active metabolomic views found in view_metadata.")
    return(NULL)
  }
  
  for (metabolomic_view in metabolomic_views) {
    message(sprintf("View: %s | Factor: %s", metabolomic_view, factor_name))
    
    # --- 1. Find KEGG species for this view ---
    sp_string <- view_metadata$Kegg_Species[match(metabolomic_view, view_metadata$View)]
    sp_list <- unlist(strsplit(sp_string, ","))
    
    for (sp in sp_list) {
      message(sprintf("species: %s", sp))
      
      if (is.na(sp) || !sp %in% names(fella_list)) {
        message("  Skipped: no FELLA object for species ", sp)
        next
      }
      
      background <- TERM2GENE_list[[sp]][["metabolite"]][["KEGG_ID"]]
      
      # --- 2. Extract weights for this factor/view ---
      kegg_weights <- MOFAobject.trained@expectations$KEGG_weights[[paste(metabolomic_view, sp, sep = "_")]]
      if (is.null(kegg_weights) || !factor_name %in% colnames(kegg_weights)) {
        message("  Skipped: no data or factor missing for ", sp)
        next
      }
      
      vec <- kegg_weights[, factor_name]
      names(vec) <- rownames(kegg_weights)
      vec <- vec[!is.na(vec)]
      if (length(vec) == 0L) {
        message("  Skipped: all weights NA.")
        next
      }
      
      # Deduplicate KEGG IDs
      vec_dedup <- tapply(vec, names(vec), mean, na.rm = TRUE)
      vec_dedup <- sort(vec_dedup[!is.na(vec_dedup)], decreasing = TRUE)
      
      input_ids <- colnames(fella_list[[sp]]@diffusion@matrix)
      valid_ids <- intersect(names(vec_dedup), input_ids)
      if (length(valid_ids) == 0L) {
        message("  Skipped: no matching KEGG IDs in FELLA DB for ", sp)
        next
      }
      
      # --- 3. Run FELLA diffusion ---
      fella.user <- FELLA::defineCompounds(
        compounds = valid_ids,
        data = fella_list[[sp]],
        compoundsBackground = background
      )
      fella.user <- FELLA::enrich(
        compounds = valid_ids, 
        method = "diffusion", 
        approx = "normality", 
        data = fella_list[[sp]])
      
      # Results
      fella.table <- FELLA::generateResultsTable(
        object = fella.user, data = fella_list[[sp]], method = "diffusion", nlimit = 50
      )
      
      output_dir <- ensure_output_dirs(factor_name, metabolomic_view)
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      write.csv(fella.table,
                file.path(output_dir, paste0("FELLA_Enrichment_", sp, ".csv")),
                row.names = FALSE)
      
      # --- 4. Generate and save FELLA enrichment graph ---
      fella_graph <- FELLA::generateResultsGraph(
        object = fella.user,
        data = fella_list[[sp]],
        method = "diffusion",
        nlimit = 100,
        threshold = 0.05
      )
      
      if (igraph::vcount(fella_graph) > 0) {
        output_path <- file.path(output_dir, paste0("FELLA_enrichment_graph_top50_", sp, ".png"))
        
        node_ids <- igraph::V(fella_graph)$name
        
        # Determine type by KEGG slot membership
        get_node_type <- function(id) {
          if (id %in% names(fella_list[[sp]]@keggdata@id$compound)) return("compound")
          if (id %in% names(fella_list[[sp]]@keggdata@id$reaction)) return("reaction")
          if (id %in% names(fella_list[[sp]]@keggdata@id$enzyme))   return("enzyme")
          if (id %in% names(fella_list[[sp]]@keggdata@id$pathway))  return("pathway")
          if (id %in% names(fella_list[[sp]]@keggdata@id$module))   return("module")
          return("unknown")
        }
        node_types <- vapply(node_ids, get_node_type, character(1))
        
        # Human-readable labels
        node_labels <- vapply(node_ids, function(id) {
          if (id %in% names(fella_list[[sp]]@keggdata@id2name)) {
            fella_list[[sp]]@keggdata@id2name[[id]][1]
          } else {
            id
          }
        }, character(1))
        
        # Map weights to node_ids, default 0 if missing
        node_weights <- sapply(node_ids, function(id) {
          if (id %in% names(vec_dedup)) {
            return(vec_dedup[[id]])
          } else {
            return(0)
          }
        })
        
        # Add attrs
        fella_graph <- igraph::set_vertex_attr(fella_graph, "type", value = node_types)
        fella_graph <- igraph::set_vertex_attr(fella_graph, "label", value = node_labels)
        fella_graph <- igraph::set_vertex_attr(fella_graph, "weight", value = node_weights)
        
        # Plot with ggraph
        p <- ggraph(fella_graph, layout = "fr") +
          geom_edge_link(alpha = 0.2, colour = "grey70") +
          geom_node_point(aes(shape = type, color = weight), size = 5, alpha = 0.9, na.rm = TRUE) +
          geom_node_text(aes(label = label), size = 2.5, repel = TRUE) +
          scale_shape_manual(values = c(compound = 15, reaction = 16, enzyme = 17, pathway = 18, module = 3, unknown = 1),
                             name = "Category") +
          scale_color_gradient2(low = "blue", mid = "grey80", high = "red", midpoint = 0, limits = c(-0.5, 0.5), name = "Weight") +
          theme_classic() +
          theme(
            axis.line = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank()
          ) +
          ggtitle(sprintf("FELLA diffusion — %s | %s (%s)", factor_name, metabolomic_view, sp))
        
        ggsave(output_path, plot = p, width = 8, height = 8, dpi = 300)
        message("Saved FELLA enrichment graph for ", factor_name,
                " (", metabolomic_view, ", ", sp, ")")
      }
    }
  }
}

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

print("Go Have a party, we are done here...")
