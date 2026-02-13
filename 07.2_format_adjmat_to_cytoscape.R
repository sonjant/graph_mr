# Format Airwave CF MR multiomics results for Cytoscape

# Clean environment
rm(list=ls())

# Paths
if (Sys.getenv("GRAPH_MR_ROOT") == "") stop("Set GRAPH_MR_ROOT to the pipeline root directory (see README).")
setwd(Sys.getenv("GRAPH_MR_ROOT"))
source("config/config.R")

# Set number of cores
ncores = 4

# Load packages
library("tidyverse")
library("vroom") 
library("readxl")
library("igraph")

# Arguments: method_type, dataset
args <- commandArgs(trailingOnly = TRUE)
method_type <- args[1]
dataset     <- if (length(args) >= 2L) args[2] else Sys.getenv("GRAPH_MR_DATASET")
if (!nzchar(dataset)) stop("Provide dataset name or set GRAPH_MR_DATASET.")
base_dir <- paste0(method_type, "_", dataset)

# Read the result adjacency matrices (unadjusted and adjusted)
# Unadjusted models
bonferroni_matrix_unadjusted <- vroom(paste0(results_dir, base_dir, "/matrices_chem_names/fstat10-bonferroni-totalfx-unadjusted-matrix.csv"), delim = ",", show_col_types = F) %>%
  as.data.frame()

nomsig_matrix_unadjusted <- vroom(paste0(results_dir, base_dir, "/matrices_chem_names/fstat10-nominal-totalfx-unadjusted-matrix.csv"), show_col_types = F) %>%
  as.data.frame()

# Adjusted models
bonferroni_matrix_adjusted <- vroom(paste0(results_dir, base_dir, "/matrices_chem_names/fstat10-bonferroni-totalfx-agesex-matrix.csv"), delim = ",", show_col_types = F) %>%
  as.data.frame()

nomsig_matrix_adjusted <- vroom(paste0(results_dir, base_dir, "/matrices_chem_names/fstat10-nominal-totalfx-agesex-matrix.csv"), show_col_types = F) %>%
  as.data.frame()

# Optional metadata for node colours (nodes in graph are display names after 06.1)
meta_path <- paste0(data_dir, base_dir, "/preprocessed/exposure_metadata.tsv")
node_metadata <- if (file.exists(meta_path)) {
  m <- vroom(meta_path, show_col_types = FALSE)
  if ("display_name" %in% names(m)) m <- dplyr::rename(m, Exposure = display_name)
  else if ("exposure_id" %in% names(m)) m <- dplyr::rename(m, Exposure = exposure_id)
  m
} else {
  data.frame(Exposure = character(), stringsAsFactors = FALSE)
}
pathway_col <- intersect(c("pathway", "class", "SUPER_PATHWAY", "Protein class"), names(node_metadata))[1]

# Set and create the output file path for the igraph objects
output_path <- paste0(results_dir, base_dir, "/cytoscape/")
dir.create(paste0(output_path), showWarnings = F, recursive = T)


# Create igraph .gml output of total_fx results for Cytoscape----
writeLines(c(" ", " ", "Started script 07.2 - Format adj. matrix for cytoscape", " ", " "))

# Function to create and write graph from adjacency matrix
create_and_write_graph <- function(adj_matrix, output_filename, input_type, output_path, suffix = "") {
  
  # Create an edge list 
  edges <- which(!is.na(as.matrix(adj_matrix[ , -which(names(adj_matrix) == "Exposure")])), arr.ind = TRUE)
  
  # Extract parent and child names for each non-NA edge
  edges <- data.frame(from = adj_matrix$Exposure[edges[, 1]],
                     to = colnames(adj_matrix)[edges[, 2] + 1],  # Adjust column index to skip Exposure column
                     color = "black")
  
  # Get unique nodes from the edge list
  nodes <- unique(c(edges$from, edges$to))
  
  # Assign colors to nodes based on input_type
  if (input_type == 'metabolomics') {
    # Get unique SUPER_PATHWAY categories
    categories <- unique(metabolomics_dictionary$SUPER_PATHWAY[!is.na(metabolomics_dictionary$SUPER_PATHWAY)])
    categories <- sort(as.character(categories))
    
    # Generate colors for each category
    n_categories <- length(categories)
    all_possible_colors <- rainbow(n_categories)
    all_colors <- setNames(all_possible_colors, categories)
    
    # Assign colors to nodes based on their SUPER_PATHWAY
    node_colors <- data.frame(
      Exposure = nodes,
      color = sapply(nodes, function(node) {
        pathway <- metabolomics_dictionary$SUPER_PATHWAY[metabolomics_dictionary$Exposure == node]
        if (length(pathway) > 0 && !is.na(pathway[1])) {
          all_colors[as.character(pathway[1])]
        } else {
          "#999999"  # Default gray color if no pathway found
        }
      })
    )
    
  } else if (input_type == 'proteomics') {
    # Get unique Protein class categories
    categories <- unique(protein_dictionary$`Protein class`[!is.na(protein_dictionary$`Protein class`)])
    categories <- sort(as.character(categories))
    
    # Generate colors for each category
    n_categories <- length(categories)
    all_possible_colors <- rainbow(n_categories)
    all_colors <- setNames(all_possible_colors, categories)
    
    # Assign colors to nodes based on their Protein class
    node_colors <- data.frame(
      Exposure = nodes,
      color = sapply(nodes, function(node) {
        protein_class <- protein_dictionary$`Protein class`[protein_dictionary$Exposure == node]
        if (length(protein_class) > 0 && !is.na(protein_class[1])) {
          all_colors[as.character(protein_class[1])]
        } else {
          "#999999"  # Default gray color if no class found
        }
      })
    )
    
  } else if (input_type == 'tripleintersection' || input_type == 'phenomolecular') {
    # Set the node categories and corresponding colours for source types
    all_possible_categories <- c("Protein", "Metabolite", "Phenotype")
    all_possible_colors <- c("#ff9933", "#0000ee", "#ff6699")
    all_colors <- setNames(all_possible_colors, all_possible_categories)
    
    # Assign colors to nodes based on which dictionary contains the Exposure
    node_colors <- data.frame(
      Exposure = nodes,
      color = sapply(nodes, function(node) {
        if (node %in% protein_dictionary$Exposure) {
          all_colors["Protein"]
        } else if (node %in% metabolomics_dictionary$Exposure) {
          all_colors["Metabolite"]
        } else if (node %in% phenomic_dictionary$Exposure) {
          all_colors["Phenotype"]
        } else {
          "#999999"  # Default gray color if not found in any dictionary
        }
      })
    )
  } else {
    # Default case: use gray for all nodes
    node_colors <- data.frame(
      Exposure = nodes,
      color = "#999999"
    )
  }
  
  # Create the graph and assign attributes
  g <- graph_from_data_frame(edges, directed = TRUE, vertices = node_colors)
  
  # Set edge colors
  E(g)$color <- edges$color
  
  # Set node colors
  V(g)$color <- V(g)$color
  
  # Write the graph to a .gml file with appropriate suffix
  if (suffix == "") {
    write_graph(g, file = paste0(output_path, output_filename), format = "gml")
  } else {
    # Insert suffix before .gml extension
    output_filename_suffixed <- gsub("\\.gml$", paste0("_", suffix, ".gml"), output_filename)
    write_graph(g, file = paste0(output_path, output_filename_suffixed), format = "gml")
  }
  
  return(g)
}

# Create and write graphs for unadjusted models
writeLines("Creating graphs for unadjusted models...")
bonferroni_graph_unadjusted <- create_and_write_graph(
  adj_matrix = bonferroni_matrix_unadjusted,
  output_filename = paste0("bonferroni_", base_dir, ".gml"),
  input_type = dataset,
  output_path = output_path,
  suffix = "unadjusted"
)

nomsig_graph_unadjusted <- create_and_write_graph(
  adj_matrix = nomsig_matrix_unadjusted,
  output_filename = paste0("nomsig_", base_dir, ".gml"),
  input_type = dataset,
  output_path = output_path,
  suffix = "unadjusted"
)

# Create and write graphs for adjusted models
writeLines("Creating graphs for adjusted models...")
bonferroni_graph_adjusted <- create_and_write_graph(
  adj_matrix = bonferroni_matrix_adjusted,
  output_filename = paste0("bonferroni_", base_dir, ".gml"),
  input_type = dataset,
  output_path = output_path,
  suffix = "agesex"
)

nomsig_graph_adjusted <- create_and_write_graph(
  adj_matrix = nomsig_matrix_adjusted,
  output_filename = paste0("nomsig_", base_dir, ".gml"),
  input_type = dataset,
  output_path = output_path,
  suffix = "agesex"
)


writeLines(c(" ", " ", "Finished script 07.2 - Formatted adj. matrix for cytoscape", " ", " "))


