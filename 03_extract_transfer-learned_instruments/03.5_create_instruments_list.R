# Purpose: Create a list of instrumental variables (IVs) for the dosage extraction.

# Clean the environment
rm(list=ls())

# Paths
if (Sys.getenv("GRAPH_MR_ROOT") == "") stop("Set GRAPH_MR_ROOT to the pipeline root directory (see README).")
setwd(Sys.getenv("GRAPH_MR_ROOT"))
source("config/config.R")

# Load packages
library("tidyverse")
library("vroom")
library("parallel")
library("gtools")

# Specify a base_dir based on the command line inputs
args <- commandArgs(trailingOnly = TRUE)
method_type <- args[1]
input_type <- args[2]
base_dir <- paste0(method_type, "_", input_type)

# Set the number of cores you require
ncores <- 2

# Create an input List of file paths to all exposure IVs (same as the output path in 03.4)
input_list <- list.files(path = paste0(results_dir, base_dir, "/instruments"), pattern = "complete_cases_instruments.tsv", full.names = T, recursive = T)

# Optional: protein annotation (set GRAPH_MR_ANNOTATION_DIR or use data/annotation/OlinkProteomics/)
olink_csv <- paste0(annotation_dir, "OlinkProteomics/olink-explore-3072-assay-list-2024-03-19.csv")
protein_names <- if (file.exists(olink_csv)) {
  vroom(olink_csv) %>% dplyr::rename(UniProt = "UniProt ID", protein_name = "Protein name") %>% dplyr::select(UniProt, protein_name)
} else {
  data.frame(UniProt = character(), protein_name = character())
}
protein38_path <- paste0(annotation_dir, "OlinkProteomics/proteinID_hg38_coordinates.txt")
protein38coords <- if (file.exists(protein38_path)) {
  vroom(protein38_path, col_names = c("UniProtID", "hg38chr", "hg38start", "hg38end")) %>% distinct(.)
} else {
  data.frame(UniProtID = character(), hg38chr = integer(), hg38start = integer(), hg38end = integer())
}
bed_path <- paste0(annotation_dir, "OlinkProteomics/protein_coords_hg19_hg38.bed")
protein_19_38_coords <- if (file.exists(bed_path)) {
  vroom(bed_path, col_names = c("hg19chr", "hg19start", "hg19end", "hg38chr_start_end", "rm")) %>% 
  dplyr::select("hg19chr", "hg19start", "hg19end", "hg38chr_start_end") %>% 
  dplyr::mutate(hg19chr = str_remove(hg19chr, "^chr")) %>% 
  tidyr::separate(hg38chr_start_end, into = c("hg38chr", "pos"), sep = ":") %>%
  tidyr::separate(pos, into = c("hg38start", "hg38end"), sep = "-") %>%
  dplyr::mutate(hg38chr = str_remove(hg38chr, "^chr"),
                hg38chr = as.integer(hg38chr),
                hg38start = as.integer(hg38start),
                hg38end = as.integer(hg38end),
                hg19chr = as.integer(hg19chr)) %>%
  distinct(.)
} else {
  data.frame(hg19chr = integer(), hg19start = integer(), hg19end = integer(), hg38chr = integer(), hg38start = integer(), hg38end = integer())
}

# Merge protein coordinates and for any duplicates, pick the longest coordinate duplicate
merged_protein_coords <- dplyr::inner_join(protein_19_38_coords, protein38coords, by = c("hg38chr", "hg38start", "hg38end")) %>% 
  dplyr::mutate(hg19_length = hg19end - hg19start) %>%      
  dplyr::group_by(UniProtID) %>%
  dplyr::slice_max(order_by = hg19_length, n = 1, with_ties = F) %>%  
  dplyr::ungroup()

# Set the output file path to the directory containing the instrumental variables, which was the "output" path in script 02.5
output <- paste0(results_dir, base_dir, "/instruments/")


# Start script ---- ---- ---- ---- ---- ----
writeLines(c(" ", " ", "Started script 03.5 - Make instruments list.", " ", " "))

# Create a list of IVs for each exposure 
create_iv_lists <- function(input_list) {
  mclapply(input_list, mc.cores = ncores, FUN = function(exposure_i) {
    # Extract the name of the input exposure 
    exposure_name <- basename(dirname(exposure_i))
    
    # Read IVs
    exposure_i_ivs <- vroom(exposure_i, delim = "\t", col_types = c("CHROM" = "i", "POS" = "i"))
    
    # Check that IVs exist
    if (nrow(exposure_i_ivs) > 0) {
      
      # Remove variants that contain >1 base
      exposure_i_ivs <- exposure_i_ivs %>%
        dplyr::filter(nchar(REF) == 1 & nchar(ALT) == 1)
      
      # Check again if IVs still exist after filtering
      if (nrow(exposure_i_ivs) > 0) {
        
        # Case 1: Protein with coordinates available — filter for cis region
        if (exposure_name %in% protein_names$UniProt && exposure_name %in% merged_protein_coords$UniProtID) {
          coords <- merged_protein_coords %>%
            dplyr::filter(UniProtID == exposure_name) %>%
            dplyr::mutate(hg19start = hg19start - 500000,
                          hg19end = hg19end + 500000)
          
          exposure_i_ivs <- exposure_i_ivs %>%
            dplyr::filter(CHROM == coords$hg19chr,
                          POS >= coords$hg19start,
                          POS <= coords$hg19end) %>%
            dplyr::mutate(list = paste0(CHROM, "\t", POS)) %>%
            dplyr::select(list)
          
          if (nrow(exposure_i_ivs) > 0) {
            write.table(exposure_i_ivs,
                        file = paste0(output, exposure_name, "/complete_cases_instruments_list.txt"),
                        col.names = FALSE, row.names = FALSE, quote = FALSE, append = FALSE)
            cat("Exposure:", exposure_name, "list saved\n")
          } else {
            cat("No SNPs in cis window for exposure:", exposure_name, "\n")
          }
          
          # Case 2: Protein but no coordinates — skip
        } else if (exposure_name %in% protein_names$UniProt) {
          cat("Skipping exposure with no coordinates:", exposure_name, "\n")
          
          # Case 3: Not a protein — save all variants
        } else {
          exposure_i_ivs <- exposure_i_ivs %>%
            dplyr::mutate(list = paste0(CHROM, "\t", POS)) %>%
            dplyr::select(list)
          
          write.table(exposure_i_ivs,
                      file = paste0(output, exposure_name, "/complete_cases_instruments_list.txt"),
                      col.names = FALSE, row.names = FALSE, quote = FALSE, append = FALSE)
          cat("Exposure:", exposure_name, "list saved\n")
        }
        
      } else {
        cat("No valid SNPs (after filtering REF/ALT) for exposure:", exposure_name, "\n")
      }
      
    } else {
      cat("No variants found for exposure:", exposure_name, "\n")
    }
    
    return(NULL) # always NULL for mclapply
  })
}


iv_lists <- create_iv_lists(input_list)

writeLines(c(" ", " ", "Finished script 03.5 - Made instruments list.", " ", " "))
