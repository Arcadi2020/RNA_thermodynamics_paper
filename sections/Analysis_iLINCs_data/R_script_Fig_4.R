# Script: Figure 4 — iLINCS cohort filtering, pathway NES heatmaps, volcano+boxplots, and per-sample MFE violins
#
# Purpose
#   Reproduce Fig. 4 from iLINCS transcriptomic perturbation data at 3.33 µM and 10 µM.
#   The script builds “growth-arrest / no-mitogenic-signaling” cohorts via GSEA filters,
#   summarizes cell line–drug coverage, and visualizes ΔG_fold effects (Normalized MFE)
#   at the cohort and per-sample levels.
#
# Data loaded
#   • Precomputed R environments:
#       - Analysis_MFE_iLINCS_all_datasets_3.33uM_cancer_lines_environment.RData
#       - Analysis_MFE_iLINCS_all_datasets_10uM_cancer_lines_environment_selected.RData
#   • GSEA outputs (per condition):
#       - Results_GSEA/combined_kk2_results.csv   # KEGG gseKEGG NES
#       - Results_GSEA/combined_gsea_results.csv  # proliferation signatures NES
#   • Per-sample MFE summaries:
#       - Results_all_cancer_cell_lines_3.33uM.csv
#       - Results_all_cancer_cell_lines_10uM.csv
#   • Metadata / lookups:
#       - cellinfo_beta.csv                       # cell line → primary_disease
#       - iLINCs_complete_data_selected_samples.csv  # Drug → TARGET, TARGET_PATHWAY
#       - Extended_Table_ANNOTATED_iLINCs/Table_ANNOTATED_iLINCs.csv  # ensembl_id, Normalized_MFE, p_value, z_score



library(readr)
library(dplyr)
library(biomaRt)
library(ggplot2)
library(ggpmisc)
library(ggrepel)
library(forcats)  
library(stringr)
library(pheatmap)
library(patchwork)
library(tidyr)
library(gridExtra)



# Load the environment of the script output
load("path_to/Analysis_MFE_iLINCS_all_datasets_3.33uM_cancer_lines_environment.RData")




# Fig. 4A ####
# Load the list object
# This is the total number of iLINCs samples that has been downloaded
directory <- "path_to/Output_script_cancer_cell_lines_3.33uM"
result_list <- readRDS(file.path(directory, "iLINCs_data_unfiltered.rds"))
# 1,609 samples



# Determine the number of unique cell lines and drugs
# Get the names of the list
sample_names <- names(result_list)

# Extract cell line names (everything before the first "_")
cell_lines <- sub("_.*", "", sample_names)
unique_cell_lines <- unique(cell_lines)
length(unique_cell_lines)
# 37 cells lines

# Extract drug names (between the first and second "_")
drug_names <- sapply(strsplit(sample_names, "_"), `[`, 2)
unique_drugs <- unique(drug_names)
length(unique_drugs)

# Optional: print out the unique names if you want to inspect them
head(unique_cell_lines)
head(unique_drugs)



# Filtering only cancer cell lines
directory <- "path_to/Output_script_cancer_cell_lines_3.33uM"
result_list_2 <- readRDS(file.path(directory, "iLINCs_filtered_samples.rds"))
# 774 samples



# Create a data frame with the names of the cell lines and drugs
iLINCs_ell_lines <- data.frame(
  iLINCs_cell_line = sapply(names(result_list_2), function(x) {
    # Extract part before the first _
    sub("_.*", "", x)
  }),
  Anti_cancer_agents = sapply(names(result_list_2), function(x) {
    # Extract part between the first and second _
    sub("^[^_]*_([^_]*).*", "\\1", x)
  }),
  original_ID = names(result_list_2),  # Add the original names as a new column
  stringsAsFactors = FALSE
)

# Remove row names
rownames(iLINCs_ell_lines) <- NULL



# How many unique cancer cell lines are present in the dataset
Unique_cancer_cell_lines <- unique(iLINCs_ell_lines$iLINCs_cell_line)
print(length(Unique_cancer_cell_lines))
# 20 cancer cell lines

# How many unique anti-cancer agents are present in the dataset
Unique_anti_cancer_agents <- unique(iLINCs_ell_lines$Anti_cancer_agents)
print(length(Unique_anti_cancer_agents))
# 40 anti-cancer agents













# The samples are subjected to GSEA. The output is NES for specific mitogenic signaling pathways and specific cell proliferation signatures 
# In the first filtering step samples with negative NES values for the cell cycle KEGG pathway (hsa04110) are seleted (totally 423).

# Open the combined_kk2 generated using the gseKEGG function of the clusterProfiler package
GSEA_combined_kk2 <- read.delim("path_to/Results_GSEA/combined_kk2_results.csv",
                                sep = ",", header = T, stringsAsFactors = F)

# Open the combined_gsea generated using the GSEA function of the clusterProfiler to assess the proliferation signatures
GSEA_combined_iLINCs <- read.delim("path_to/Results_GSEA/combined_gsea_results.csv",
                                   sep = ",", header = T, stringsAsFactors = F)

# Join the GSEA_combined_kk2 data frame with the results of the proliferation signatures GSEA
iLINCs_complete_GSEA <- dplyr::bind_rows(GSEA_combined_iLINCs, GSEA_combined_kk2)

# Identify Sample with NES values for all three proliferation signatures
required_ids <- c("hsa04110", "Signt_1", "Signt_2", "Signt_3")

# Get Sample names that have all three IDs
valid_samples <- iLINCs_complete_GSEA %>%
  filter(ID %in% required_ids) %>%
  group_by(Sample) %>%
  filter(all(required_ids %in% ID)) %>%
  pull(Sample) %>%
  unique()

# Filter the original data frame to include only the valid samples
iLINCs_complete_GSEA <- iLINCs_complete_GSEA %>%
  filter(Sample %in% valid_samples)



# Selected KEGG signaling pathways and proliferation signatures
Selected_KEGG_pathways <- c("hsa04150", "hsa04062", "hsa04668", "hsa04657", "hsa04064", "hsa04110", "hsa03030", "hsa03410", "hsa03430", "hsa03420",
                            "hsa03320", "hsa03460", "hsa03440", "hsa04014", "hsa04066", "hsa04218", "hsa04512", "hsa04024", "hsa04071", "hsa04070",
                            "hsa04068", "hsa04310", "hsa03082", "hsa04015", "hsa04625", "hsa04371", "hsa04390", "hsa04010", "hsa04622", "hsa04350",
                            "hsa04152", "hsa04926", "hsa04630", "hsa04072", "hsa04115", "hsa04151", "hsa04910", "hsa04330", "hsa04012", "hsa04022",
                            "hsa04020", "hsa04370", "hsa04210", "hsa03450", "Signt_1", "Signt_2", "Signt_3")

GSEA_combined_kk2_filtered <- dplyr::filter(iLINCs_complete_GSEA, ID %in% Selected_KEGG_pathways)
Samples_3.33uM_N1 <- unique(GSEA_combined_kk2_filtered$Sample)
# 57 samples





# The second filtering step retains only samples which have also negative NES values for the cell proliferatin signatures (totally 407)
# Identify all samples with a significant down-regulation of the cell cycle and negative NES values of the proliferation signatures
# Identify the samples where ID is "hsa04110", "Signt_1", "Signt_2", or "Signt_3" and NES >= 0
iLINCs_proliferating <- GSEA_combined_kk2_filtered[(GSEA_combined_kk2_filtered$ID %in% c("hsa04110", "Signt_1", "Signt_2", "Signt_3") & GSEA_combined_kk2_filtered$NES >= 0), ]
iLINCs_proliferating <- dplyr::distinct(iLINCs_proliferating, Sample, .keep_all = TRUE)

# Remove these samples from the complete dataset
GSEA_combined_kk2_filtered_cell_cycle_arrest <- GSEA_combined_kk2_filtered[!(GSEA_combined_kk2_filtered$Sample %in% iLINCs_proliferating$Sample), ]
Samples_3.33uM_N2 <- unique(GSEA_combined_kk2_filtered_cell_cycle_arrest$Sample)
# 50 samples






# The third filtering step removes all samples which have positive NES values for any mitogenic signaling pathway (totally 216)
# Identify the samples with a significant down-regulation of cell cycle and mitogenic signaling pathways
# Check all the pathways that may foster cell proliferation
ids_to_check <- c(
  "hsa04150", "hsa04062", "hsa04668", "hsa04657", "hsa04064", "hsa03320", 
  "hsa04014", "hsa04066", "hsa04068", "hsa04310", "hsa04015", "hsa04070",
  "hsa04625", "hsa04371", "hsa04390", "hsa04010", "hsa04622", "hsa04350", 
  "hsa04152", "hsa04926", "hsa04630", "hsa04072", "hsa04151", "hsa04512",
  "hsa04910", "hsa04330", "hsa04012", "hsa04022", "hsa04020", "hsa03030",
  "hsa03420", "hsa04071"
)

# Identify samples with NES > 0 for the specified IDs
samples_to_exclude <- GSEA_combined_kk2_filtered_cell_cycle_arrest %>%
  dplyr::filter(ID %in% ids_to_check, NES >= 0) %>%
  dplyr::select(Sample) %>%
  dplyr::distinct()

# Remove these samples from the original data frame
iLINCS_cell_cycle_arrest_mitogenic_signaling <- dplyr::filter(GSEA_combined_kk2_filtered_cell_cycle_arrest, 
                                                              !Sample %in% samples_to_exclude$Sample)

Samples_3.33uM_N3 <- unique(iLINCS_cell_cycle_arrest_mitogenic_signaling$Sample)
# 24 samples













# Select samples with a significant down-regulation of proliferation signatures and  mitogenic signaling pathways

ids_to_check <- c(
  "hsa04150", "hsa04062", "hsa04668", "hsa04657", "hsa04064", "hsa03320", 
  "hsa04014", "hsa04066", "hsa04068", "hsa04310", "hsa04015", "hsa04070",
  "hsa04625", "hsa04371", "hsa04390", "hsa04010", "hsa04622", "hsa04350", 
  "hsa04152", "hsa04926", "hsa04630", "hsa04072", "hsa04151", "hsa04512",
  "hsa04910", "hsa04330", "hsa04012", "hsa04022", "hsa04020", "hsa03030",
  "hsa03420", "hsa04071"
)

# Identify samples with NES > 0 for the specified IDs
samples_to_exclude <- GSEA_combined_kk2_filtered_cell_cycle_arrest %>%
  dplyr::filter(ID %in% ids_to_check, NES >= 0) %>%
  dplyr::select(Sample) %>%
  dplyr::distinct()

# Remove these samples from the original data frame
iLINCS_cell_cycle_arrest_mitogenic_signaling <- dplyr::filter(GSEA_combined_kk2_filtered_cell_cycle_arrest, 
                                                              !Sample %in% samples_to_exclude$Sample)



# Load library
# Load KEGGREST
library(KEGGREST)


# 1. IDs to check
ids_to_check <- c(
  "hsa04150", "hsa04062", "hsa04668", "hsa04657", "hsa04064", "hsa03320", 
  "hsa04014", "hsa04066", "hsa04068", "hsa04310", "hsa04015", "hsa04070",
  "hsa04625", "hsa04371", "hsa04390", "hsa04010", "hsa04622", "hsa04350", 
  "hsa04152", "hsa04926", "hsa04630", "hsa04072", "hsa04151", "hsa04512",
  "hsa04910", "hsa04330", "hsa04012", "hsa04022", "hsa04020", "hsa03030",
  "hsa03420", "hsa04071", "Signt_1", "Signt_2", "Signt_3"
)

# 2. Exclude proliferating samples
samples_to_exclude <- GSEA_combined_kk2_filtered_cell_cycle_arrest %>%
  dplyr::filter(ID %in% ids_to_check, NES >= 0) %>%
  dplyr::select(Sample) %>%
  dplyr::distinct()

iLINCS_cell_cycle_arrest_mitogenic_signaling <- dplyr::filter(
  GSEA_combined_kk2_filtered_cell_cycle_arrest, 
  !Sample %in% samples_to_exclude$Sample
)

# 3. Retrieve KEGG pathway names via KEGGREST
library(KEGGREST)

# Separate IDs (exclude Signt_* for KEGG query)
kegg_ids <- ids_to_check[!grepl("Signt", ids_to_check)]

kegg_info <- lapply(kegg_ids, function(id) {
  info <- keggGet(id)
  data.frame(ID = id, Name = info[[1]]$NAME, stringsAsFactors = FALSE)
})
kegg_info_df <- do.call(rbind, kegg_info)

# Add Signt_* manually
extra_signatures <- data.frame(
  ID = c("Signt_1", "Signt_2", "Signt_3"),
  Name = c("Signt_1", "Signt_2", "Signt_3"),
  stringsAsFactors = FALSE
)
kegg_info_df <- rbind(kegg_info_df, extra_signatures)

# Clean names (remove suffix)
kegg_info_df$Name <- sub(" - Homo sapiens \\(human\\)$", "", kegg_info_df$Name)
kegg_info_df$Name <- sub(" signaling pathway$", " signaling", kegg_info_df$Name)

# 4. Explicit order from your figure
ordered_names <- c(
  "AMPK signaling",
  "Apelin signaling",
  "Apoptosis",
  "ATP-chromatin remodel.",
  "Base excision repair",
  "C-type lectin receptor signaling",
  "Cellular senescence",
  "cGMP-PKG signaling",
  "Chemokine signaling",
  "DNA replication",
  "ECM-receptor interaction",
  "ErbB signaling",
  "Fanconi anemia",
  "FoxO signaling",
  "HIF-1 signaling",
  "Hippo signaling",
  "Homologous recombination",
  "IL-17 signaling",
  "JAK-STAT signaling",
  "MAPK signaling",
  "Mismatch repair",
  "mTOR signaling",
  "NF-kappa B signaling",
  "Notch signaling",
  "Nucleotide excision repair",
  "p53 signaling",
  "PI3K-Akt signaling",
  "PPAR signaling",
  "Rap1 signaling",
  "Relaxin signaling",
  "RIG-I-like receptor signaling",
  "TGF-beta signaling",
  "TNF signaling",
  "Wnt signaling",
  "Cell cycle",
  "Signt_1",
  "Signt_2",
  "Signt_3"
)

# 5. Generate heatmap
generate_heatmap <- function(kk2_df, kegg_info_df, ordered_names) {
  # Clean KEGG names
  kegg_info_df$Name <- sub(" - Homo sapiens \\(human\\)$", "", kegg_info_df$Name)
  kegg_info_df$Name <- sub(" signaling pathway$", " signaling", kegg_info_df$Name)
  
  # Map IDs to cleaned pathway names
  kk2_df <- kk2_df %>%
    dplyr::mutate(Description = kegg_info_df$Name[match(ID, kegg_info_df$ID)]) %>%
    dplyr::filter(!is.na(Description))
  
  # Filter to only pathways in ordered_names
  kk2_df <- kk2_df %>% dplyr::filter(Description %in% ordered_names)
  
  # Factor with explicit order
  kk2_df$Description <- factor(kk2_df$Description, levels = ordered_names)
  
  # Remove duplicates and reshape
  kk2_df <- kk2_df %>%
    dplyr::select(Description, Sample, NES) %>%
    dplyr::distinct(Description, Sample, .keep_all = TRUE) %>%
    tidyr::spread(Sample, NES)
  
  # Replace NA with 0
  kk2_df[is.na(kk2_df)] <- 0
  
  # Convert Description to character for rownames
  kk2_df$Description <- as.character(kk2_df$Description)
  rownames(kk2_df) <- kk2_df$Description
  kk2_df <- kk2_df[, -1, drop = FALSE]
  
  # ADD MISSING PATHWAYS WITH ZEROS
  missing_paths <- setdiff(ordered_names, rownames(kk2_df))
  if (length(missing_paths) > 0) {
    zero_matrix <- matrix(0, nrow = length(missing_paths), ncol = ncol(kk2_df))
    rownames(zero_matrix) <- missing_paths
    colnames(zero_matrix) <- colnames(kk2_df)
    kk2_df <- rbind(kk2_df, zero_matrix)
  }
  
  # Reorder rows to match ordered_names
  kk2_df <- kk2_df[ordered_names, , drop = FALSE]
  
  # Generate heatmap
  pheatmap(as.matrix(kk2_df),
           color = colorRampPalette(c("blue", "white"))(50),
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           display_numbers = TRUE,
           number_format = "%.2g",
           fontsize_row = 11,
           fontsize_col = 2,
           fontsize_number = 0,
           border_color = "black",
           lwd = 0.3)
}



# 6. Plot heatmap
Heatmap_cell_cycle <- generate_heatmap(iLINCS_cell_cycle_arrest_mitogenic_signaling,
                                       kegg_info_df,
                                       ordered_names)
print(Heatmap_cell_cycle)


Samples_3.33uM_N3 <- unique(iLINCS_cell_cycle_arrest_mitogenic_signaling$Sample)
# 24 samples







# Check the total number of cell lines and drug combinations in this final dataset
Cell_line_drug_combinations <- iLINCS_cell_cycle_arrest_mitogenic_signaling %>%
  dplyr::mutate(Sample = gsub("(^[^_]+_[^_]+)_.*", "\\1", Sample)) %>%
  dplyr::distinct(Sample)
print(nrow(Cell_line_drug_combinations))
# 24 combinations

# Write the final dataset to a file
iLINCs_Cell_line_drug_combinations_final <-  Cell_line_drug_combinations %>%
  separate(Sample, into = c("Cell_line", "Drug"), sep = "_", remove = FALSE)
output_directory <- "path_to/Fig_4/Other_files"
write.csv(iLINCs_Cell_line_drug_combinations_final, file = paste0(output_directory, "/iLINCs_Cell_line_drug_combinations_3.33uM.csv"), row.names = FALSE)


# How many different cell lines are present in the final dataset
Cell_lines <- iLINCS_cell_cycle_arrest_mitogenic_signaling %>%
  dplyr::mutate(Sample = gsub("(^[^_]+)_.*", "\\1", Sample)) %>%
  dplyr::distinct(Sample)
print(nrow(Cell_lines))
# 7 human cancer cell lines

# How many drugs are present in the final dataset
Drugs <- iLINCS_cell_cycle_arrest_mitogenic_signaling %>%
  dplyr::mutate(Sample = gsub("^[^_]+_([^_]+)_.*", "\\1", Sample)) %>%
  dplyr::distinct(Sample)
print(nrow(Drugs))
# 17 different drugs










# Load the complete iLINCs 3.33 uM data
all_results_2_complete <- read.csv("path_to/Results_all_cancer_cell_lines_3.33uM.csv",
                                   header = TRUE, sep = ",", stringsAsFactors = FALSE)


# Add a new column called 'Log10_MFE_change' containing the log10 of the ratio between the mean MFE of the DOWN-regulated genes and the mean MFE of the UP-regulated genes
library(dplyr)
library(tidyr)

iLINCS_analysis <- all_results_2_complete %>%
  mutate(
    Log10_MFE_change = log10(MeanMFE_DOWN / MeanMFE_UP)
  ) %>%
  relocate(MeanMFE_UP, .after = MeanMFE_DOWN) %>%
  relocate(Log10_MFE_change, .after = MeanMFE_UP) %>%
  # Create temporary split columns but keep SampleName
  separate(SampleName, into = c("cell_line", "Drug", "rest"), sep = "_", extra = "merge", remove = FALSE) %>%
  dplyr::select(-rest)






# Check the total number of cell lines and drug combinations in this final dataset
Cell_line_drug_combinations <- iLINCS_analysis %>%
  dplyr::distinct(SampleName)
print(nrow(Cell_line_drug_combinations))
# 24 combinations


# How many different cell lines are present in the final dataset
Cell_lines <- iLINCS_analysis %>%
  dplyr::distinct(cell_line)
print(nrow(Cell_lines))
# 7 human cancer cell lines

# How many drugs are present in the final dataset
Drugs <- iLINCS_analysis %>%
  dplyr::distinct(Drug)
print(nrow(Drugs))
# 16 different drugs









# Add information on cell lines and drugs and group the samples accordingly

# Open the iLINCs cell_info file
iLINCs_cell_info <- read.delim("path_to/cellinfo_beta.csv",
                               sep = ",", header = T, stringsAsFactors = F) %>%
  dplyr::rename(iLINCs_cell_line = cell_iname)


iLINCS_analysis <- dplyr::rename(iLINCS_analysis, iLINCs_cell_line = cell_line)

# Merge the iLINCS_analysis data frame with iLINCs_cell_info on the iLINCs_cell_line and cell_iname columns
iLINCs_cell_info_complete <- dplyr::left_join(iLINCS_analysis, iLINCs_cell_info, by = c("iLINCs_cell_line", "primary_disease")) %>%
  dplyr::select(iLINCs_cell_line, Drug, primary_disease)
  # dplyr::mutate(primary_disease = ifelse(iLINCs_cell_line == "P1A82", "teratocarcinoma", primary_disease))




# Add the drug targets
# In the original data frame, some of the drug targets, were missing. They have been added together with the relative signaling pathways. The file has been save as .csv file and loaded here below.
iLINCs_drug_info_complete <- read.delim("path_to/iLINCs_complete_data_selected_samples.csv",
                                        sep = ",", header = T, stringsAsFactors = F) %>%
  dplyr::select("Anti_cancer_agents", "TARGET", "TARGET_PATHWAY") %>%
  dplyr::rename(Drug = Anti_cancer_agents) %>%
  dplyr::distinct(Drug, .keep_all = TRUE)



# Join the data frame to generate the complete plot
IlINCs_complete_df <- iLINCs_cell_info_complete %>%
  dplyr::left_join(iLINCs_drug_info_complete, by = "Drug") %>%
  dplyr::mutate(TARGET_PATHWAY = str_replace(TARGET_PATHWAY, 'ATR-Chk1 DNA damage response', 'DNA damage response'))





# Plot the data
# Create a new column with primary_disease for grouping iLINCs_cell_lines
IlINCs_complete_df <- IlINCs_complete_df %>%
  mutate(primary_disease_label = factor(primary_disease, levels = unique(primary_disease)),
         TARGET_PATHWAY_label = factor(TARGET_PATHWAY, levels = unique(TARGET_PATHWAY)))


# Total n. of cell lines:
Total_cell_lines <- unique(IlINCs_complete_df$iLINCs_cell_line)
print(length(Total_cell_lines))

# Total n. of anti-cancer agents:
Total_anti_cancer_agents <- unique(IlINCs_complete_df$Drug)
print(length(Total_anti_cancer_agents))


# Create the tile plot
Tile_plot_iLINCs_cell_lines_drugs_grouped <- ggplot(IlINCs_complete_df, aes(x = Drug, y = iLINCs_cell_line)) +
  geom_tile(fill = "blue", color = "white") +  # Create tile plot with colored tiles
  theme_bw() +
  labs(x = "Anti-cancer Agents", y = "iLINCs Cell Lines", title = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.ticks.x = element_blank(),  # Remove x-axis ticks
        axis.ticks.y = element_blank(),  # Remove y-axis ticks
        panel.grid = element_blank(),  # Remove grid lines
        strip.text.x = element_text(size = 12, angle = 90, hjust = 0, vjust = 0.5),  # Align TARGET_PATHWAY labels
        strip.text.y = element_text(size = 12, angle = 0, hjust = 0, vjust = 0.5),  # Align primary_disease labels
        strip.background = element_blank(),  # Remove facet background
        panel.spacing = unit(0.15, "lines")) +  # Reduce space between groups
  facet_grid(primary_disease ~ TARGET_PATHWAY, scales = "free", space = "free")  # Group by primary_disease and TARGET_PATHWAY

# Print the plot
print(Tile_plot_iLINCs_cell_lines_drugs_grouped)




# How many different cell lines are present in the final datset:
iLINCS_cell_lines_3.33uM_unique <- unique(IlINCs_complete_df$iLINCs_cell_line)
length(iLINCS_cell_lines_3.33uM_unique)

# How many cancer types does the dataset include:
iLINCS_cancer_types_3.33uM_unique <- unique(IlINCs_complete_df$primary_disease)
length(iLINCS_cancer_types_3.33uM_unique)

# How many distinct drugs does the dataset contain:
iLINCS_drugs_3.33uM_unique <- unique(IlINCs_complete_df$Drug)
length(iLINCS_drugs_3.33uM_unique)















# Fig. 4B ####

# iLINCs volcano plot

# iLINCs volcano plot all samples
# Generate a volcano plot with different colors for selected samples
iLINCS_analysis_plot_membrane_receptors <- ggplot(iLINCS_analysis, aes(x = Log10_MFE_change, y = -log10(wilcox_test_all_UP_vs_DOWN_pValue))) +
  theme_classic() +
  geom_point(aes(
    shape = ifelse(Log10_MFE_change > 0 & -log10(wilcox_test_all_UP_vs_DOWN_pValue) > 1.301, 16, 21),
    fill  = ifelse(Log10_MFE_change > 0 & -log10(wilcox_test_all_UP_vs_DOWN_pValue) > 1.301, "blue", "white"),
    color = ifelse(Log10_MFE_change > 0 & -log10(wilcox_test_all_UP_vs_DOWN_pValue) > 1.301, "blue", "black")
  ),
  size = 4, stroke = 0.8) +
  scale_shape_identity() +  # Use identity for shape mapping
  scale_fill_identity() +   # Use identity for fill mapping
  scale_color_identity() +  # Use identity for color mapping
  scale_alpha_identity() +  # Use identity for alpha mapping
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Add vertical line at x = 0
  geom_hline(yintercept = 1.301, linetype = "dashed", color = "black") +  # Add horizontal line at y = 1.301
  labs(x = "Log10_MFE_change", y = "-Log10_ks_test_all_UP_vs_DOWN_pValue", title = "") +
  theme(plot.title = element_text(color = "black", size = 14, face = "bold"),
        legend.position = "none") +
  xlim(-0.06, 0.06) +
  geom_text_repel(aes(label = ifelse(-log10(wilcox_test_all_UP_vs_DOWN_pValue) > 1.301, paste(iLINCs_cell_line, Drug, sep = "_"), "")), 
                  size = 3, box.padding = 0.5, point.padding = 0.5, nudge_x = 0.004, nudge_y = 0.001, max.overlaps = 15)

print(iLINCS_analysis_plot_membrane_receptors)









# Plot with selected names for the manuscript
iLINCS_analysis_plot_membrane_receptors <- ggplot(iLINCS_analysis, aes(
  x = iLINCS_analysis$Log10_MFE_change,
  y = -log10(iLINCS_analysis$wilcox_test_all_UP_vs_DOWN_pValue))) +
  theme_classic() +
  geom_point(aes(
    fill  = ifelse(Log10_MFE_change > 0 & -log10(wilcox_test_all_UP_vs_DOWN_pValue) > 1.301, "blue", "lightgrey"),
    color = ifelse(Log10_MFE_change > 0 & -log10(wilcox_test_all_UP_vs_DOWN_pValue) > 1.301, "black", "black"),
    alpha = ifelse(Log10_MFE_change > 0 & -log10(wilcox_test_all_UP_vs_DOWN_pValue) > 1.301, 1, 0.5) # conditional transparency
  ),
  size = 4,
  shape = 21,
  stroke = 0.6) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_alpha_identity() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1.301, linetype = "dashed", color = "black") +
  labs(x = expression(Log[10]~normalized~Delta*G[fold]~"change"), 
       y = expression(-Log[10]~p.value),
       title = "") +
  theme(plot.title = element_text(color = "black", size = 14, face = "bold"),
        legend.position = "none") +
  ylim(0, 8) +
  xlim(-0.06, 0.06) +
  geom_text_repel(aes(
    label = ifelse(-log10(wilcox_test_all_UP_vs_DOWN_pValue) > 1.301 & 
                     paste(iLINCs_cell_line, Drug, sep = "_") %in% 
                     c("HCC515_neratinib", "HT29_trametinib", "PC3_palbociclib"),
                   paste(iLINCs_cell_line, Drug, sep = "_"), "")), 
    size = 3, box.padding = 1, point.padding = 0.5,
    nudge_x = 0.004, nudge_y = 0.001, max.overlaps = Inf, 
    segment.color = 'grey50', segment.size = 0.5)

print(iLINCS_analysis_plot_membrane_receptors)

# Save the plot
ggsave("iLINCS_analysis_volcano_plot.pdf", plot = iLINCS_analysis_plot_membrane_receptors, width = 5, height = 8, path = "path_to/Desktop")









# Classify samples into quadrants
library(dplyr)
library(ggplot2)
library(ggrepel)

thr <- 1.301

iLINCS_analysis <- iLINCS_analysis %>%
  mutate(
    Quadrant = case_when(
      Log10_MFE_change > 0 & -log10(wilcox_test_all_UP_vs_DOWN_pValue) >  thr ~ "UP_regulated_high_MFE_TX",
      Log10_MFE_change > 0 & -log10(wilcox_test_all_UP_vs_DOWN_pValue) <= thr ~ "Q2",
      Log10_MFE_change < 0 & -log10(wilcox_test_all_UP_vs_DOWN_pValue) >  thr ~ "DOWN_regulated_high_MFE_TX",
      Log10_MFE_change < 0 & -log10(wilcox_test_all_UP_vs_DOWN_pValue) <= thr ~ "Q4",
      TRUE ~ "Other"
    )
  )

# Counts & percentages
n_total <- nrow(iLINCS_analysis)

# UP (either significant or not)
N_total_iLINCs_UP <- iLINCS_analysis %>%
  filter(Quadrant %in% c("UP_regulated_high_MFE_TX", "Q2")) %>%
  nrow()
Total_iLINCs_UP <- 100 * N_total_iLINCs_UP / n_total

# UP significant only
N_iLINCs_UP_pval <- iLINCS_analysis %>%
  filter(Quadrant == "UP_regulated_high_MFE_TX") %>%
  nrow()
iLINCs_pval_UP <- 100 * N_iLINCs_UP_pval / n_total
Percentage_iLINCs_UP <- 100 * N_iLINCs_UP_pval / N_total_iLINCs_UP

# DOWN (either significant or not)
N_total_iLINCs_DOWN <- iLINCS_analysis %>%
  filter(Quadrant %in% c("DOWN_regulated_high_MFE_TX", "Q4")) %>%
  nrow()
Total_iLINCs_DOWN <- 100 * N_total_iLINCs_DOWN / n_total

# DOWN significant only
N_iLINCs_DOWN_pval <- iLINCS_analysis %>%
  filter(Quadrant == "DOWN_regulated_high_MFE_TX") %>%
  nrow()
iLINCs_pval_DOWN <- 100 * N_iLINCs_DOWN_pval / n_total
Percentage_iLINCs_DOWN <- 100 * N_iLINCs_DOWN_pval / N_total_iLINCs_DOWN

# Print the results
summary_df <- data.frame(
  Category = c("UP (sig+ns)", "UP (sig only)", "DOWN (sig+ns)", "DOWN (sig only)"),
  Description = c(
    "All UP-regulated samples (includes both significant and non-significant)",
    "Significant UP-regulated samples (Log10_MFE_change > 0 & p < 0.05)",
    "All DOWN-regulated samples (includes both significant and non-significant)",
    "Significant DOWN-regulated samples (Log10_MFE_change < 0 & p < 0.05)"
  ),
  N = c(
    N_total_iLINCs_UP,
    N_iLINCs_UP_pval,
    N_total_iLINCs_DOWN,
    N_iLINCs_DOWN_pval
  ),
  Percent = c(
    round(Total_iLINCs_UP, 2),
    round(iLINCs_pval_UP, 2),
    round(Total_iLINCs_DOWN, 2),
    round(iLINCs_pval_DOWN, 2)
  )
)

print(summary_df)













# Plot to show the distribution of Log10_MFE_change
# Generate the horizontal box plot with reduced box width and adjusted x-axis limits
iLINCS_analysis_boxplot_horizontal <- ggplot(iLINCS_analysis, aes(y = "", x = Log10_MFE_change)) +
  theme_classic() +
  geom_boxplot(fill = "lightgrey", color = "black", width = 0.8) +
  #geom_jitter(height = 0.18, color = "black", size = 3, alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.8, color = "black") +  # Add vertical line at x = 0
  labs(y = "", x = "", title = "") +
  xlim(-0.06, 0.06) +  # Adjust x-axis limits here
  theme(plot.title = element_text(color = "black", size = 14, face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

print(iLINCS_analysis_boxplot_horizontal)

# Save the plot
ggsave("iLINCS_analysis_boxplot_horizontal.pdf", plot = iLINCS_analysis_boxplot_horizontal, width = 5, height = 1.5, path = "path_to/Desktop")


# Null hypothesis (H0) the mean of the distribution is 0.
# Alternative hypothesis (HA) the mean of the distribution is greater than 0.
# Perform a One-Sample t-Test:

test_result <- t.test(iLINCS_analysis$Log10_MFE_change, mu = 0, alternative = "greater")
print(test_result)










# Combine the two plots into a single figure
iLINCS_analysis_boxplot_horizontal <- ggplot(iLINCS_analysis, aes(y = "", x = Log10_MFE_change)) +
  theme_classic() +
  geom_boxplot(fill = "lightgrey", color = "black", width = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.8, color = "black") +  # Add vertical line at x = 0
  labs(y = "", x = "", title = "") +
  xlim(-0.06, 0.06) +  # Align with scatter plot
  theme(plot.title = element_text(color = "black", size = 14, face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),  # Remove x-axis labels
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())  # Remove x-axis ticks

# Combine the two plots with a reduced boxplot height (70% reduction)
combined_plot <- iLINCS_analysis_boxplot_horizontal / iLINCS_analysis_plot_membrane_receptors + 
  plot_layout(heights = c(0.07, 0.93))  # 20% for boxplot, 80% for scatter plot

# Display the combined plot
print(combined_plot)

# Save the combined plot
ggsave("Combined_iLINCS_analysis_3.33_uM.tiff", plot = combined_plot, width = 5, height = 8.5, path = "path_to/Desktop", dpi = 300)




# Null hypothesis (H0) the mean of the distribution is 0.
# Alternative hypothesis (HA) the mean of the distribution is greater than 0.
# Perform a One-Sample t-Test:

test_result <- t.test(iLINCS_analysis$Log10_MFE_change, mu = 0, alternative = "greater")
print(test_result)


wilcox_result <- wilcox.test(iLINCS_analysis$Log10_MFE_change, mu = 0, alternative = "greater", exact = FALSE)
print(wilcox_result)



















# Fig. 4C ####

# Violin plots samples 3.33 uM

# Load the complete data generated in STEP 7 the R script:
# path_to/Analysis_MFE_iLINCS_all_datasets_10uM_cancer_lines_AC.R

# Load the list object from the specified directory (STEP 7 of the original R script)
directory <- "path_to/Output_script_cancer_cell_lines_3.33uM"
result_list_2 <- readRDS(file.path(directory, "iLINCs_filtered_samples.rds"))




# Select the samples of interest:
all_results_2_complete <- read.csv("path_to/Results_all_cancer_cell_lines_3.33uM.csv",
                                   header = TRUE, sep = ",", stringsAsFactors = FALSE)

# Edit the names of the samples to match the iLINCs format
all_results_2_complete$iLINCs_cell_line_drug <- sub("^(([^_]*_[^_]*)_).*", "\\1", all_results_2_complete$SampleName)

# Remove the trailing underscore
all_results_2_complete$iLINCs_cell_line_drug <- gsub("_$", "", all_results_2_complete$iLINCs_cell_line_drug)

iLINCs_all_results_2_complete_filtered <- dplyr::filter(all_results_2_complete,  iLINCs_cell_line_drug %in% c("HCC515_neratinib", "HT29_trametinib", "PC3_palbociclib")) %>%
  dplyr::filter(wilcox_test_all_UP_vs_DOWN_pValue < 0.05) %>% # Filter for significant samples
  dplyr::mutate(SampleName = str_remove(SampleName, "_MFE\\.pdf$"))


# Filter the result_list_2 list object to only include the selected samples
iLINCs_selected_samples <- result_list_2[grep(paste(iLINCs_all_results_2_complete_filtered$SampleName, collapse = "|"), names(result_list_2))]




# Open the Table_ANNOTATED_iLINCs
Table_ANNOTATED_iLINCs <- read.delim("path_to/Table_ANNOTATED_iLINCs.csv",
                                     sep = ",", header = T, stringsAsFactors = F)
Table_ANNOTATED_iLINCs <- dplyr::rename(Table_ANNOTATED_iLINCs, ensembl_id = ensembl_gene_id)



# Initialize plot_list as an empty list
plot_list <- list()

# Initialize all_results as an empty data frame with appropriate column names
all_results <- data.frame(SampleName = character(),
                          Net_MFE_change = numeric(),
                          MeanMFE_UP = numeric(),
                          MeanMFE_DOWN = numeric(),
                          TTEST_pValue = numeric(),
                          ks_test_all_UP_vs_DOWN_pValue = numeric(),
                          wilcox_test_all_UP_vs_DOWN_pValue = numeric(),
                          stringsAsFactors = FALSE)

# Function to analyze the data and store plots
process_result_list <- function(Sample_name) { 
  message("Processing ", Sample_name)
  df <- iLINCs_selected_samples[[Sample_name]] # Access individual data frames
  
  # Modify title_plot_1 to remove text after the second underscore
  split_title <- strsplit(Sample_name, "_")[[1]]
  if (length(split_title) > 2) {
    title_plot_1 <- paste(split_title[1:2], collapse = "_")
  } else {
    title_plot_1 <- Sample_name  # If there are less than 2 underscores, keep the original
  }
  
  title_plot_1 <- paste0(title_plot_1)  # Optionally add the rest of the title
  
  # Join the data with Table_ANNOTATED_iLINCs
  Sample_name_TX <- dplyr::right_join(df, Table_ANNOTATED_iLINCs, by = "ensembl_id")
  
  # Filter only statistically significant z-scores
  Sample_name_TX_filtered <- dplyr::filter(Sample_name_TX, p_value < 0.05)
  
  # Calculate means for up-regulated and down-regulated genes
  MFE_Sample_name_UP <- dplyr::filter(Sample_name_TX_filtered, z_score > 0)
  Mean_MFE_Sample_name_UP <- mean(MFE_Sample_name_UP$Normalized_MFE, na.rm = TRUE)
  
  MFE_Sample_name_DOWN <- dplyr::filter(Sample_name_TX_filtered, z_score < 0)
  Mean_MFE_Sample_name_DOWN <- mean(MFE_Sample_name_DOWN$Normalized_MFE, na.rm = TRUE)
  
  # Perform statistical tests
  TTEST <- t.test(MFE_Sample_name_DOWN$Normalized_MFE, MFE_Sample_name_UP$Normalized_MFE, paired = FALSE)
  ks_test_all_UP_vs_DOWN <- ks.test(MFE_Sample_name_DOWN$Normalized_MFE, MFE_Sample_name_UP$Normalized_MFE)
  # Perform Wilcoxon test (Mann–Whitney U test)
  wilcox_test_all_UP_vs_DOWN <- wilcox.test(
    MFE_Sample_name_DOWN$Normalized_MFE,
    MFE_Sample_name_UP$Normalized_MFE,
    alternative = "two.sided",
    paired = FALSE,
    exact = FALSE  # safer for large sample sizes
  )
  
  # Calculate the net MFE change
  Sample_name_TX_net_MFE_change <- Sample_name_TX_filtered %>% 
    dplyr::mutate(Net_MFE_change = (Normalized_MFE * z_score) * (-1))
  Net_MFE_change <- sum(Sample_name_TX_net_MFE_change$Net_MFE_change, na.rm = TRUE)
  
  # Generate violin plot with a white box plot inside
  MFE_Sample_name <- ggplot(Sample_name_TX_filtered, aes(x = factor(z_score > 0), y = Normalized_MFE)) +
    geom_violin(fill = "lightgrey", alpha = 0.7) +
    geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = NA) +  # Box plot on top
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 20),
      axis.title.y = element_blank(),              # Keep this if you want no label, or remove it to add a label
      axis.text.x = element_blank(),
      # axis.text.y = element_blank(),             # REMOVE this line to show y-axis values
      plot.title = element_blank()
    ) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    labs(y = "Normalized MFE", x = NULL) +          # Now we add a label to y-axis if you want it
    scale_x_discrete(labels = c("z_score < 0", "z_score > 0"))
  
  
  # Save the plot
  ggsave(
    filename = paste0(Sample_name, "_MFE_plot.pdf"),  # Create a unique filename
    plot = MFE_Sample_name,
    width = 7,
    height = 5,
    path = "path_to/Desktop"  # Specify the path to save the plot
  )
  
  # Store the plot in the list
  plot_list[[Sample_name]] <<- MFE_Sample_name
  
  # Create a data frame with the results for the current sample
  sample_results <- data.frame(
    SampleName = Sample_name,  # Use Sample_name here
    Net_MFE_change = Net_MFE_change,
    MeanMFE_UP = Mean_MFE_Sample_name_UP,
    MeanMFE_DOWN = Mean_MFE_Sample_name_DOWN,
    TTEST_pValue = TTEST$p.value,
    ks_test_all_UP_vs_DOWN_pValue = ks_test_all_UP_vs_DOWN$p.value,
    wilcox_test_all_UP_vs_DOWN_pValue = wilcox_test_all_UP_vs_DOWN$p.value
  )
  
  # Append the sample results to the all_results data frame
  all_results <<- rbind(all_results, sample_results)
}

# Run the function for each iLINCS data frame
Data_result_list <- lapply(names(iLINCs_selected_samples), process_result_list)

# Arrange and display all violin plots in a single figure with 3 rows and increased vertical spacing
library(patchwork)

# Option 1 — simplest: use plot_spacer() between plots
combined_plot <- plot_list[[1]] / 
  plot_spacer() / 
  plot_list[[2]] / 
  plot_spacer() / 
  plot_list[[3]] +
  plot_layout(ncol = 1, heights = c(1, 0.1, 1, 0.1, 1)) & 
  theme(plot.margin = margin(20, 10, 20, 10))

# Save combined plot
ggsave(
  filename = "combined_violin_plots.pdf",
  plot = combined_plot,
  width = 5,
  height = 16,  # Increased height to better accommodate spacing
  path = "path_to/Desktop"
)

print(all_results)


















# ANALYSIS iLINCs DATA 10 µM ####
# The data have been generated in the R script which analyses the iLINCs datasets:
# path_to/Analysis_MFE_iLINCS_all_datasets_10uM_cancer_lines_AC.R


# Load the environment of the script output
load("path_to/Analysis_MFE_iLINCS_all_datasets_10uM_cancer_lines_environment_selected.RData")


library(readr)
library(dplyr)
library(biomaRt)
library(ggplot2)
library(ggpmisc)
library(ggrepel)
library(forcats)  
library(stringr)
library(pheatmap)
library(patchwork)
library(tidyr)
library(gridExtra)




# Fig. 4D ####

# Open the combined_kk2 generated using the gseKEGG function of the clusterProfiler package
# These data show the samples that have ≥1 pathway with a significant NES value
GSEA_combined_kk2 <- read.delim("path_to/combined_kk2_results.csv",
                                sep = ",", header = T, stringsAsFactors = F)

# Open the combined_gsea generated using the GSEA function of the clusterProfiler to assess the proliferation signatures
# These data show the samples that have ≥1 proliferation signature with a significant NES value 
GSEA_combined_iLINCs <- read.delim("path_to/combined_gsea_results.csv",
                                   sep = ",", header = T, stringsAsFactors = F)

# Join the GSEA_combined_kk2 data frame with the results of the proliferation signatures GSEA
iLINCs_complete_GSEA <- dplyr::bind_rows(GSEA_combined_iLINCs, GSEA_combined_kk2)

# Identify Sample with significant NES values for all three proliferation signatures and for the cell cycle KEGG pathway
required_ids <- c("hsa04110", "Signt_1", "Signt_2", "Signt_3")

# Get Sample names that have all three IDs
valid_samples <- iLINCs_complete_GSEA %>%
  filter(ID %in% required_ids) %>%
  group_by(Sample) %>%
  filter(all(required_ids %in% ID)) %>%
  pull(Sample) %>%
  unique()





# Filter the original data frame to include only the valid samples
iLINCs_complete_GSEA <- iLINCs_complete_GSEA %>%
  filter(Sample %in% valid_samples)







# Selected KEGG signaling pathways and proliferation signatures
Selected_KEGG_pathways <- c("hsa04150", "hsa04062", "hsa04668", "hsa04657", "hsa04064", "hsa04110", "hsa03030", "hsa03410", "hsa03430", "hsa03420",
                            "hsa03320", "hsa03460", "hsa03440", "hsa04014", "hsa04066", "hsa04218", "hsa04512", "hsa04024", "hsa04071", "hsa04070",
                            "hsa04068", "hsa04310", "hsa03082", "hsa04015", "hsa04625", "hsa04371", "hsa04390", "hsa04010", "hsa04622", "hsa04350",
                            "hsa04152", "hsa04926", "hsa04630", "hsa04072", "hsa04115", "hsa04151", "hsa04910", "hsa04330", "hsa04012", "hsa04022",
                            "hsa04020", "hsa04370", "hsa04210", "hsa03450", "Signt_1", "Signt_2", "Signt_3")

GSEA_combined_kk2_filtered <- dplyr::filter(iLINCs_complete_GSEA, ID %in% Selected_KEGG_pathways)
Samples_N1 <- unique(GSEA_combined_kk2_filtered$Sample)
# Totally 259 samples that have significant NES values in each of the three proliferation signatures, the cell cycle KEGG pathway and any other KEGG pathways





# The second filtering step retains only samples which have also negative NES values for the cell proliferation signatures (totally 407)
# Identify all samples with a significant down-regulation of the cell cycle and negative NES values of the proliferation signatures
# Identify the samples where ID is "hsa04110", "Signt_1", "Signt_2", or "Signt_3" and NES >= 0
iLINCs_proliferating <- GSEA_combined_kk2_filtered[(GSEA_combined_kk2_filtered$ID %in% c("hsa04110", "Signt_1", "Signt_2", "Signt_3") & GSEA_combined_kk2_filtered$NES >= 0), ]
iLINCs_proliferating <- dplyr::distinct(iLINCs_proliferating, Sample, .keep_all = TRUE)

# Remove these samples from the complete dataset
GSEA_combined_kk2_filtered_cell_cycle_arrest <- GSEA_combined_kk2_filtered[!(GSEA_combined_kk2_filtered$Sample %in% iLINCs_proliferating$Sample), ]
Samples_N2 <- unique(GSEA_combined_kk2_filtered_cell_cycle_arrest$Sample)
# 244 samples









# The third filtering step removes all samples which have positive NES values for any mitogenic signaling pathway (totally 216)
# Identify the samples with a significant down-regulation of cell cycle and mitogenic signaling pathways
# Check all the pathways that may foster cell proliferation
ids_to_check <- c(
  "hsa04150", "hsa04062", "hsa04668", "hsa04657", "hsa04064", "hsa03320", 
  "hsa04014", "hsa04066", "hsa04068", "hsa04310", "hsa04015", "hsa04070",
  "hsa04625", "hsa04371", "hsa04390", "hsa04010", "hsa04622", "hsa04350", 
  "hsa04152", "hsa04926", "hsa04630", "hsa04072", "hsa04151", "hsa04512",
  "hsa04910", "hsa04330", "hsa04012", "hsa04022", "hsa04020", "hsa03030",
  "hsa03420", "hsa04071"
)

# Identify samples with NES > 0 for the specified IDs
samples_to_exclude <- GSEA_combined_kk2_filtered_cell_cycle_arrest %>%
  dplyr::filter(ID %in% ids_to_check, NES >= 0) %>%
  dplyr::select(Sample) %>%
  dplyr::distinct()

# Remove these samples from the original data frame
iLINCS_cell_cycle_arrest_mitogenic_signaling <- dplyr::filter(GSEA_combined_kk2_filtered_cell_cycle_arrest, 
                                                              !Sample %in% samples_to_exclude$Sample)

Samples_N3 <- unique(iLINCS_cell_cycle_arrest_mitogenic_signaling$Sample)
# 147 samples








# Check the total number of cell lines and drug combinations in this final dataset
Cell_line_drug_combinations <- iLINCS_cell_cycle_arrest_mitogenic_signaling %>%
  mutate(Sample = gsub("(^[^_]+_[^_]+)_.*", "\\1", Sample)) %>%
  dplyr::distinct(Sample)
# 125 cell line and drug combinations

# Write the final dataset to a file
iLINCs_Cell_line_drug_combinations_final <-  Cell_line_drug_combinations %>%
  separate(Sample, into = c("Cell_line", "Drug"), sep = "_", remove = FALSE)
output_directory <- "path_to/Fig_3"
write.csv(iLINCs_Cell_line_drug_combinations_final, file = paste0(output_directory, "/iLINCs_Cell_line_drug_combinations.csv"), row.names = FALSE)


# How many different cell lines are present in the final dataset
Cell_lines <- iLINCS_cell_cycle_arrest_mitogenic_signaling %>%
  mutate(Sample = gsub("(^[^_]+)_.*", "\\1", Sample)) %>%
  dplyr::distinct(Sample)
# 36 human cancer cell lines

# How many drugs are present in the final dataset
Drugs <- iLINCS_cell_cycle_arrest_mitogenic_signaling %>%
  mutate(Sample = gsub("^[^_]+_([^_]+)_.*", "\\1", Sample)) %>%
  dplyr::distinct(Sample)
# 33 different drugs
























# Select samples with a significant down-regulation of proliferation signatures and  mitogenic signaling pathways

ids_to_check <- c(
  "hsa04150", "hsa04062", "hsa04668", "hsa04657", "hsa04064", "hsa03320", 
  "hsa04014", "hsa04066", "hsa04068", "hsa04310", "hsa04015", "hsa04070",
  "hsa04625", "hsa04371", "hsa04390", "hsa04010", "hsa04622", "hsa04350", 
  "hsa04152", "hsa04926", "hsa04630", "hsa04072", "hsa04151", "hsa04512",
  "hsa04910", "hsa04330", "hsa04012", "hsa04022", "hsa04020", "hsa03030",
  "hsa03420", "hsa04071"
)

# Identify samples with NES > 0 for the specified IDs
samples_to_exclude <- GSEA_combined_kk2_filtered_cell_cycle_arrest %>%
  dplyr::filter(ID %in% ids_to_check, NES >= 0) %>%
  dplyr::select(Sample) %>%
  dplyr::distinct()

# Remove these samples from the original data frame
iLINCS_cell_cycle_arrest_mitogenic_signaling <- dplyr::filter(
  GSEA_combined_kk2_filtered_cell_cycle_arrest, 
  !Sample %in% samples_to_exclude$Sample
)

# Define the heatmap generation function with row reordering
generate_heatmap <- function(kk2_df) {
  
  # Prepare data for the heatmap
  kk2_df <- kk2_df %>%
    dplyr::select(Description, Sample, NES) %>%
    tidyr::spread(Sample, NES)
  
  # Replace NA values with 0 for better visualization
  kk2_df[is.na(kk2_df)] <- 0
  
  # Define desired order for Description: all others first, Signt_1, Signt_2, Signt_3 last
  desired_order <- c(
    setdiff(kk2_df$Description, c("Signt_1", "Signt_2", "Signt_3")),
    "Signt_1", "Signt_2", "Signt_3"
  )
  
  # Set factor and reorder rows accordingly
  kk2_df$Description <- factor(kk2_df$Description, levels = desired_order)
  kk2_df <- kk2_df %>% dplyr::arrange(Description)
  
  # Set row names and remove Description column
  rownames(kk2_df) <- kk2_df$Description
  kk2_df <- kk2_df[,-1]
  
  # Generate heatmap
  pheatmap(as.matrix(kk2_df), 
           color = colorRampPalette(c("blue", "white"))(50), 
           cluster_rows = FALSE, 
           cluster_cols = FALSE, 
           display_numbers = TRUE, 
           number_format = "%.2g", 
           main = "",
           fontsize_row = 11, 
           fontsize_col = 2,
           fontsize_number = 0,
           border_color = "black",
           labels_col = NA)
}


# Generate the heatmap using the combined_kk2_df
Heatmap_cell_cycle <- generate_heatmap(iLINCS_cell_cycle_arrest_mitogenic_signaling)
print(Heatmap_cell_cycle)


# Check the number of samples that show cell cycle arrest and no mitogenic signaling
Samples_N3 <- unique(iLINCS_cell_cycle_arrest_mitogenic_signaling$Sample)
# 143 samples

# Save the plot
ggsave("Heatmap_cell_cycle.pdf", plot = Heatmap_cell_cycle, width = 12, height = 14, path = "path_to/Desktop")
























# Fig. 4E ####

# Load the complete iLINCs data
all_results_2_complete <- read.csv("path_to/Results_all_cancer_cell_lines_10uM.csv",
                                   header = TRUE, sep = ",", stringsAsFactors = FALSE)


# Add a new column called 'Log10_MFE_change' containing the log10 of the ratio between the mean MFE of the DOWN-regulated genes and the mean MFE of the UP-regulated genes
iLINCS_analysis <- dplyr::mutate(all_results_2_complete, Log10_MFE_change = log10(MeanMFE_DOWN/MeanMFE_UP)) %>%
  dplyr::relocate(MeanMFE_UP, .after = MeanMFE_DOWN) %>%
  dplyr::relocate(Log10_MFE_change, .after = MeanMFE_UP)






# Add information on cell lines and drugs and group the samples accordingly

# Open the iLINCs cell_info file
iLINCs_cell_info <- read.delim("path_to/cellinfo_beta.csv",
                               sep = ",", header = T, stringsAsFactors = F) %>%
  dplyr::rename(iLINCs_cell_line = cell_iname)


iLINCS_analysis <- dplyr::rename(iLINCS_analysis, iLINCs_cell_line = cell_line) %>%
  dplyr::select(-primary_disease)

# Merge the iLINCS_analysis data frame with iLINCs_cell_info on the iLINCs_cell_line and cell_iname columns
iLINCs_cell_info_complete <- dplyr::left_join(iLINCS_analysis, iLINCs_cell_info, by = "iLINCs_cell_line") %>%
  dplyr::select(iLINCs_cell_line, Drug, primary_disease) %>%
  dplyr::mutate(primary_disease = ifelse(iLINCs_cell_line == "P1A82", "teratocarcinoma", primary_disease))




# Add the drug targets
# In the original data frame, some of the drug targets, were missing. They have been added together with the relative signaling pathways. The file has been save as .csv file and loaded here below.
iLINCs_drug_info_complete <- read.delim("path_to/iLINCs_complete_data_selected_samples.csv",
                                        sep = ",", header = T, stringsAsFactors = F) %>%
  dplyr::select("Anti_cancer_agents", "TARGET", "TARGET_PATHWAY") %>%
  dplyr::rename(Drug = Anti_cancer_agents)



# Join the data frame to generate the complete plot
IlINCs_complete_df <- iLINCs_cell_info_complete %>%
  dplyr::left_join(iLINCs_drug_info_complete, by = c("Drug"), relationship = "many-to-many") %>%
  dplyr::mutate(TARGET_PATHWAY = str_replace(TARGET_PATHWAY, 'ATR-Chk1 DNA damage response', 'DNA damage response')) %>%
  dplyr::distinct(iLINCs_cell_line, Drug, primary_disease, TARGET, TARGET_PATHWAY)




# Plot the data
# Create a new column with primary_disease for grouping iLINCs_cell_lines
IlINCs_complete_df <- IlINCs_complete_df %>%
  mutate(primary_disease_label = factor(primary_disease, levels = unique(primary_disease)),
         TARGET_PATHWAY_label = factor(TARGET_PATHWAY, levels = unique(TARGET_PATHWAY)))


# Total n. of cell lines:
Total_cell_lines <- unique(IlINCs_complete_df$iLINCs_cell_line) # 98
# Total n. of anti-cancer agents:
Total_anti_cancer_agents <- unique(IlINCs_complete_df$Drug) # 42







# iLINCs volcano plot all samples
# Generate a volcano plot with different colors for selected samples
iLINCS_analysis_plot_membrane_receptors <- ggplot(iLINCS_analysis, aes(x = Log10_MFE_change, y = -log10(wilcox_test_all_UP_vs_DOWN_pValue))) +
  theme_classic() +
  # Define points with different aesthetics based on conditions
  geom_point(aes(shape = ifelse(-log10(wilcox_test_all_UP_vs_DOWN_pValue) < 1.301, 21, 
                                ifelse(Log10_MFE_change > 0 & -log10(wilcox_test_all_UP_vs_DOWN_pValue) > 1.301, 16, 16)),
                 fill = ifelse(-log10(wilcox_test_all_UP_vs_DOWN_pValue) < 1.301, "white", 
                               ifelse(Log10_MFE_change > 0 & -log10(wilcox_test_all_UP_vs_DOWN_pValue) > 1.301, "blue", "black")),
                 color = ifelse(-log10(wilcox_test_all_UP_vs_DOWN_pValue) < 1.301, "black", 
                                ifelse(Log10_MFE_change > 0 & -log10(wilcox_test_all_UP_vs_DOWN_pValue) > 1.301, "blue", "black")),
                 alpha = ifelse(Log10_MFE_change < 0, 0.6, 1)),  # Alpha for shading
             size = 4, stroke = 0.8) +  # Stroke for empty circles' border
  scale_shape_identity() +  # Use identity for shape mapping
  scale_fill_identity() +   # Use identity for fill mapping
  scale_color_identity() +  # Use identity for color mapping
  scale_alpha_identity() +  # Use identity for alpha mapping
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Add vertical line at x = 0
  geom_hline(yintercept = 1.301, linetype = "dashed", color = "black") +  # Add horizontal line at y = 1.301
  labs(x = "Log10_MFE_change", y = "-Log10_ks_test_all_UP_vs_DOWN_pValue", title = "") +
  theme(plot.title = element_text(color = "black", size = 14, face = "bold"),
        legend.position = "none") +
  xlim(-0.05, 0.06) +
  geom_text_repel(aes(label = ifelse(-log10(wilcox_test_all_UP_vs_DOWN_pValue) > 1.301, paste(iLINCs_cell_line, Drug, sep = "_"), "")), 
                  size = 3, box.padding = 0.5, point.padding = 0.5, nudge_x = 0.004, nudge_y = 0.001, max.overlaps = 20)

print(iLINCS_analysis_plot_membrane_receptors)







iLINCS_analysis_plot_membrane_receptors <- ggplot(iLINCS_analysis, aes(x = Log10_MFE_change, y = -log10(iLINCS_analysis$wilcox_test_all_UP_vs_DOWN_pValue))) +
  theme_classic() +
  geom_point(aes(
    fill = ifelse(Log10_MFE_change > 0 & -log10(iLINCS_analysis$wilcox_test_all_UP_vs_DOWN_pValue) > 1.301, "blue", 
                  ifelse(Log10_MFE_change > 0 & -log10(iLINCS_analysis$wilcox_test_all_UP_vs_DOWN_pValue) < 1.301, "lightgrey",
                         ifelse(Log10_MFE_change < 0 & -log10(iLINCS_analysis$wilcox_test_all_UP_vs_DOWN_pValue) > 1.301, "grey50", 
                                ifelse(Log10_MFE_change < 0 & -log10(iLINCS_analysis$wilcox_test_all_UP_vs_DOWN_pValue) < 1.301, "lightgrey", "white"))))),
    color = "black",  # Border color for all points
    size = 4,  # Set a larger fixed size for all points
    shape = 21,  # Set a consistent shape for all points
    alpha = ifelse(-log10(iLINCS_analysis$wilcox_test_all_UP_vs_DOWN_pValue) < 1.301, 0.6, 1),
    stroke = 0.6) +  # Stroke for empty circles' border
  scale_fill_identity() +
  scale_color_identity() +
  scale_alpha_identity() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1.301, linetype = "dashed", color = "black") +
  labs(x = "Log10_MFE_change", y = "-Log10_ks_test_all_UP_vs_DOWN_pValue", title = "") +
  theme(plot.title = element_text(color = "black", size = 14, face = "bold"),
        legend.position = "none") +
  xlim(-0.05, 0.06) +
  geom_text_repel(aes(label = ifelse(-log10(iLINCS_analysis$wilcox_test_all_UP_vs_DOWN_pValue) > 1.301 & 
                                       paste(iLINCs_cell_line, Drug, sep = "_") %in% c("ES2_lapatinib", "PC3_torin-1", "A375_tozasertib"),
                                     paste(iLINCs_cell_line, Drug, sep = "_"), "")), 
                  size = 3, box.padding = 0.5, point.padding = 0.5, nudge_x = 0.004, nudge_y = 0.001, max.overlaps = Inf, 
                  segment.color = 'grey50', segment.size = 0.5)

print(iLINCS_analysis_plot_membrane_receptors)









# Total n. of cell lines with a UP-regulation of high-MFE transcripts either with significant or non-significant p values

iLINCS_analysis <- iLINCS_analysis %>%
  mutate(Quadrant = case_when(
    Log10_MFE_change > 0 & -log10(wilcox_test_all_UP_vs_DOWN_pValue) > 1.301 ~ "UP_regulated_high_MFE_TX",
    Log10_MFE_change > 0 & -log10(wilcox_test_all_UP_vs_DOWN_pValue) < 1.301 ~ "Q2",
    Log10_MFE_change < 0 & -log10(wilcox_test_all_UP_vs_DOWN_pValue) > 1.301 ~ "DOWN_regulated_high_MFE_TX",
    Log10_MFE_change < 0 & -log10(wilcox_test_all_UP_vs_DOWN_pValue) < 1.301 ~ "Q4",
    TRUE ~ "Other"
  ))

# Step 2: Count samples in each quadrant and filter for UP_regulated_high_MFE_TX and DOWN_regulated_high_MFE_TX
quadrant_counts <- iLINCS_analysis %>%
  filter(Quadrant %in% c("UP_regulated_high_MFE_TX", "DOWN_regulated_high_MFE_TX")) %>%
  group_by(Quadrant) %>%
  summarize(Count = n(), .groups = 'drop')


N_total_iLINCs_UP <- iLINCS_analysis %>%
  dplyr::filter(Quadrant %in% c("UP_regulated_high_MFE_TX", "Q2")) %>%
  nrow()
Total_iLINCs_UP <- (100/(nrow(iLINCS_analysis))*N_total_iLINCs_UP)
print(Total_iLINCs_UP)
# 62.5% of the samples corresponding to 91 cell lines 


# Cell lines with a significant UP-regulation of high MFE transcripts
N_iLINCs_UP_p.val <- dplyr::filter(iLINCS_analysis, Quadrant == "UP_regulated_high_MFE_TX") %>%
  nrow()
iLINCs_p.val_UP <- (100/(nrow(iLINCS_analysis))*N_iLINCs_UP_p.val)
print(iLINCs_p.val_UP)
# 27.2% of the total number of samples has a significant UP-regulation of high MFE transcripts

Perecentage_iLINCs_UP <- (100/(Total_iLINCs_UP)*iLINCs_p.val_UP)
print(Perecentage_iLINCs_UP)
# This corresponds to 43.4% of the samples with am up-regulation of high MFE transcripts





# Total n. of ell lines with a DOWN-regulation of high-MFE transcripts either with significant or non-significant p values
N_total_iLINCs_DOWN <- iLINCS_analysis %>%
  dplyr::filter(Quadrant %in% c("DOWN_regulated_high_MFE_TX", "Q4")) %>%
  nrow()
print(N_total_iLINCs_DOWN)
# 55 cell lines

Total_iLINCs_DOWN <- (100/(nrow(iLINCS_analysis))*N_total_iLINCs_DOWN)
print(Total_iLINCs_DOWN)
# Corresponding to 37.4% of the total number of samples


# Cell lines with a significant DOWN-regulation of high MFE transcripts
N_iLINCs_DOWN_p.val <- dplyr::filter(iLINCS_analysis, Quadrant == "DOWN_regulated_high_MFE_TX") %>%
  nrow()
print(N_iLINCs_DOWN_p.val)
# 13 cell lines with a significant DOWN-regulation of high MFE transcripts

Perecentage_iLINCs_DOWN <- (100/(N_total_iLINCs_DOWN)*N_iLINCs_DOWN_p.val)
print(Perecentage_iLINCs_DOWN)
# 23.6% of the samples with a DOWN-regulation of high MFE transcripts has significant p values

iLINCs_p.val_DOWN <- (100/(nrow(iLINCS_analysis))*N_iLINCs_DOWN_p.val)
print(iLINCs_p.val_DOWN)
# 8.57% of all the samples has a significant DOWN-regulation of high MFE transcripts











# Plot selected for the figure

library(ggrepel)

# Add the labeling of the samples with a significant up-regulation of the heat capacity

# Filter the subset of samples with a significant up-regulation of transcripts with increased heta capacity
iLINCS_heat_capacity_UP_p_val_plot <- dplyr::filter(iLINCS_analysis, Quadrant == "UP_regulated_high_MFE_TX") %>%
  dplyr::filter(Log10_Cp_37C_change > 0 & wilcox_test_all_UP_vs_DOWN_Cp_37C < 0.05)

# Create the plot
iLINCS_analysis_plot_membrane_receptors <- ggplot(iLINCS_analysis, aes(x = Log10_MFE_change, y = -log10(iLINCS_analysis$wilcox_test_all_UP_vs_DOWN_pValue))) +
  theme_classic() +
  
  # First layer
  geom_point(aes(
    fill = ifelse(Log10_MFE_change > 0 & -log10(wilcox_test_all_UP_vs_DOWN_pValue) > 1.301, "blue", 
                  ifelse(paste(iLINCs_cell_line, Drug, sep = "_") %in% 
                           paste(iLINCS_heat_capacity_UP_p_val_plot$iLINCs_cell_line, iLINCS_heat_capacity_UP_p_val_plot$Drug, sep = "_"), 
                         NA,
                         ifelse(Log10_MFE_change > 0 & -log10(wilcox_test_all_UP_vs_DOWN_pValue) < 1.301, "lightgrey",
                                ifelse(Log10_MFE_change < 0 & -log10(wilcox_test_all_UP_vs_DOWN_pValue) > 1.301, "lightgrey", 
                                       ifelse(Log10_MFE_change < 0 & -log10(wilcox_test_all_UP_vs_DOWN_pValue) < 1.301, "lightgrey", "white")))))),
    color = "black",  
    size = 4,  
    shape = 21,  
    alpha = ifelse(-log10(iLINCS_analysis$wilcox_test_all_UP_vs_DOWN_pValue) < 1.301, 0.6, 0.8),
    stroke = 0.6
  ) +
  
  # Second layer
  geom_point(data = iLINCS_heat_capacity_UP_p_val_plot, aes(
    x = Log10_MFE_change, 
    y = -log10(wilcox_test_all_UP_vs_DOWN_pValue), 
    fill = "blue"),  
    color = "black",  
    size = 4,  
    shape = 21,  
    stroke = 0.6
  ) +
  
  scale_fill_identity() +
  scale_color_identity() +
  scale_alpha_identity() +
  
  # Add lines
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1.301, linetype = "dashed", color = "black") +
  
  # Labels
  labs(x = expression(Log[10]~normalized~Delta*G[fold]~"change"), 
       y = expression(-Log[10]~p.value),
       title = "") +
  theme(plot.title = element_text(color = "black", size = 14, face = "bold"),
        legend.position = "none") +
  xlim(-0.06, 0.06) +
  ylim(0, 8) +
  
  # ---- label ONLY the top-right (blue) quadrant and only the three requested samples ----
geom_text_repel(
  aes(
    label = ifelse(
      Log10_MFE_change > 0 &
        -log10(wilcox_test_all_UP_vs_DOWN_pValue) > 1.301 &
        paste(iLINCs_cell_line, Drug, sep = "_") %in% 
        c("ES2_lapatinib", "PC3_torin-1", "A375_tozasertib"),
      paste(iLINCs_cell_line, Drug, sep = " "),
      ""
    )
  ),
  size = 3,
  box.padding = 1.0,    # extra space around the text
  point.padding = 1.1,  # keeps labels off nearby dots
  force = 6,            # stronger repulsion
  force_pull = 1,
  max.time = 5,
  min.segment.length = 0,
  nudge_x = 0.002,
  nudge_y = 0.05,
  max.overlaps = Inf,
  segment.color = "grey50",
  segment.size = 0.5,
  seed = 42
)

# Print the plot
print(iLINCS_analysis_plot_membrane_receptors)














# Plot to show the distribution of Log10_MFE_change
# Generate the horizontal box plot with reduced box width and adjusted x-axis limits
iLINCS_analysis_boxplot_horizontal <- ggplot(iLINCS_analysis, aes(y = "", x = Log10_MFE_change)) +
  theme_classic() +
  geom_boxplot(fill = "lightgrey", color = "black", width = 0.8) +
  #geom_jitter(height = 0.18, color = "black", size = 3, alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.8, color = "black") +  # Add vertical line at x = 0
  labs(y = "", x = "", title = "") +
  xlim(-0.06, 0.06) +  # Adjust x-axis limits here
  theme(plot.title = element_text(color = "black", size = 14, face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

print(iLINCS_analysis_boxplot_horizontal)

# Save the plot
ggsave("iLINCS_analysis_boxplot_horizontal.pdf", plot = iLINCS_analysis_boxplot_horizontal, width = 5, height = 1.5, path = "path_to/Desktop")


# Null hypothesis (H0) the mean of the distribution is 0.
# Alternative hypothesis (HA) the mean of the distribution is greater than 0.
# Perform a One-Sample t-Test:

test_result <- t.test(iLINCS_analysis$Log10_MFE_change, mu = 0, alternative = "greater")
print(test_result)






# Combine the two plots into a single figure
iLINCS_analysis_boxplot_horizontal <- ggplot(iLINCS_analysis, aes(y = "", x = Log10_MFE_change)) +
  theme_classic() +
  geom_boxplot(fill = "lightgrey", color = "black", width = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.8, color = "black") +  # Add vertical line at x = 0
  labs(y = "", x = "", title = "") +
  xlim(-0.06, 0.06) +  # Align with scatter plot
  theme(plot.title = element_text(color = "black", size = 14, face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),  # Remove x-axis labels
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())  # Remove x-axis ticks

# Combine the two plots with a reduced boxplot height (70% reduction)
combined_plot <- iLINCS_analysis_boxplot_horizontal / iLINCS_analysis_plot_membrane_receptors + 
  plot_layout(heights = c(0.07, 0.93))  # 20% for boxplot, 80% for scatter plot

# Display the combined plot
print(combined_plot)

# Save the combined plot
ggsave("combined_iLINCS_analysis_10_uM.tiff", plot = combined_plot, width = 5, height = 8.5, path = "path_to/Desktop", dpi = 300)

# Null hypothesis (H0): the mean of the distribution is 0.
# Alternative hypothesis (HA): the mean of the distribution is greater than 0.

# One-sample t-test (parametric)
t_test_result <- t.test(iLINCS_analysis$Log10_MFE_change, mu = 0, alternative = "greater")
print(t_test_result)

# Wilcoxon signed-rank test (non-parametric, recommended for skewed data)
# H0: median(Log10_MFE_change) = 0
# HA: median(Log10_MFE_change) > 0
wilcox_result <- wilcox.test(iLINCS_analysis$Log10_MFE_change, mu = 0, alternative = "greater", exact = FALSE)
print(wilcox_result)






# How many different cell lines are present in the final datset:
iLINCS_cell_lines_10uM_unique <- unique(iLINCS_analysis$iLINCs_cell_line)
length(iLINCS_cell_lines_10uM_unique)

# How many cancer types does the dataset include:
iLINCS_cancer_types_10uM_unique <- unique(iLINCS_analysis$subtype)
length(iLINCS_cancer_types_10uM_unique)

# How many distinct drugs does the dataset contain:
iLINCS_drugs_10uM_unique <- unique(iLINCS_analysis$Drug)
length(iLINCS_drugs_10uM_unique)


















# Fig. 4F ####

# Edit the names of the samples to match the iLINCs format
all_results_2_complete$iLINCs_cell_line_drug <- sub("^(([^_]*_[^_]*)_).*", "\\1", all_results_2_complete$SampleName)

# Remove the trailing underscore
all_results_2_complete$iLINCs_cell_line_drug <- gsub("_$", "", all_results_2_complete$iLINCs_cell_line_drug)

iLINCs_all_results_2_complete_filtered <- dplyr::filter(all_results_2_complete,  iLINCs_cell_line_drug %in% c("PC3_torin-1", "A375_tozasertib", "ES2_lapatinib")) %>%
  dplyr::filter(ks_test_all_UP_vs_DOWN_pValue < 0.05)  # Filter for significant samples


# Filter the result_list_2 list object to only include the selected samples
iLINCs_selected_samples <- result_list_2[grep(paste(iLINCs_all_results_2_complete_filtered$SampleName, collapse = "|"), names(result_list_2))]

# Open the Table_ANNOTATED_iLINCs
Table_ANNOTATED_iLINCs <- read.delim("path_to/Table_ANNOTATED_iLINCs.csv",
                                     sep = ",", header = T, stringsAsFactors = F)
Table_ANNOTATED_iLINCs <- dplyr::rename(Table_ANNOTATED_iLINCs, ensembl_id = ensembl_gene_id)




# Generate the data of the normalized_MFE values
# Initialize plot_list as an empty list
plot_list <- list()

# Initialize all_results as an empty data frame with appropriate column names
all_results <- data.frame(SampleName = character(),
                          Net_MFE_change = numeric(),
                          MeanMFE_UP = numeric(),
                          MeanMFE_DOWN = numeric(),
                          TTEST_pValue = numeric(),
                          ks_test_all_UP_vs_DOWN_pValue = numeric(),
                          wilcox_test_all_UP_vs_DOWN_pValue = numeric(),
                          stringsAsFactors = FALSE)

# Function to analyze the data and store plots
process_result_list <- function(Sample_name) { 
  message("Processing ", Sample_name)
  df <- iLINCs_selected_samples[[Sample_name]] # Access individual data frames
  
  # Modify title_plot_1 to remove text after the second underscore
  split_title <- strsplit(Sample_name, "_")[[1]]
  if (length(split_title) > 2) {
    title_plot_1 <- paste(split_title[1:2], collapse = "_")
  } else {
    title_plot_1 <- Sample_name  # If there are less than 2 underscores, keep the original
  }
  
  title_plot_1 <- paste0(title_plot_1)  # Optionally add the rest of the title
  
  # Join the data with Table_ANNOTATED_iLINCs
  Sample_name_TX <- dplyr::right_join(df, Table_ANNOTATED_iLINCs, by = "ensembl_id")
  
  # Filter only statistically significant z-scores
  Sample_name_TX_filtered <- dplyr::filter(Sample_name_TX, p_value < 0.05)
  
  # Calculate means for up-regulated and down-regulated genes
  MFE_Sample_name_UP <- dplyr::filter(Sample_name_TX_filtered, z_score > 0)
  Mean_MFE_Sample_name_UP <- mean(MFE_Sample_name_UP$Normalized_MFE, na.rm = TRUE)
  
  MFE_Sample_name_DOWN <- dplyr::filter(Sample_name_TX_filtered, z_score < 0)
  Mean_MFE_Sample_name_DOWN <- mean(MFE_Sample_name_DOWN$Normalized_MFE, na.rm = TRUE)
  
  # Perform statistical tests
  TTEST <- t.test(MFE_Sample_name_DOWN$Normalized_MFE, MFE_Sample_name_UP$Normalized_MFE, paired = FALSE)
  ks_test_all_UP_vs_DOWN <- ks.test(MFE_Sample_name_DOWN$Normalized_MFE, MFE_Sample_name_UP$Normalized_MFE)
  wilcox_test_all_UP_vs_DOWN <- wilcox.test(MFE_Sample_name_DOWN$Normalized_MFE,
                                            MFE_Sample_name_UP$Normalized_MFE)  # Added Wilcoxon test
  
  # Calculate the net MFE change
  Sample_name_TX_net_MFE_change <- Sample_name_TX_filtered %>% 
    dplyr::mutate(Net_MFE_change = (Normalized_MFE * z_score) * (-1))
  Net_MFE_change <- sum(Sample_name_TX_net_MFE_change$Net_MFE_change, na.rm = TRUE)
  
  # Generate violin plot with a white box plot inside
  MFE_Sample_name <- ggplot(Sample_name_TX_filtered, aes(x = factor(z_score > 0), y = Normalized_MFE)) +
    geom_violin(fill = "lightgrey", alpha = 0.7) +
    geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = NA) +  # Box plot on top
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 20),
      axis.text.x = element_blank(),  # Removes x-axis text
      plot.title = element_text(size = 10)
    ) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    labs(y = "", x = NULL, title = title_plot_1) +  # Use modified title
    scale_x_discrete(labels = c("z_score < 0", "z_score > 0"))
  
  # Save the plot individually for each sample
  ggsave(
    filename = paste0(Sample_name, "_MFE_plot.pdf"),  # Create a unique filename
    plot = MFE_Sample_name,
    width = 5,
    height = 4,
    path = "path_to/Desktop"  # Specify the path to save the plot
  )
  
  # Store the plot in the list
  plot_list[[Sample_name]] <<- MFE_Sample_name
  
  # Create a data frame with the results for the current sample
  sample_results <- data.frame(
    SampleName = Sample_name,  # Use Sample_name here
    Net_MFE_change = Net_MFE_change,
    MeanMFE_UP = Mean_MFE_Sample_name_UP,
    MeanMFE_DOWN = Mean_MFE_Sample_name_DOWN,
    TTEST_pValue = TTEST$p.value,
    ks_test_all_UP_vs_DOWN_pValue = ks_test_all_UP_vs_DOWN$p.value,
    wilcox_test_all_UP_vs_DOWN_pValue = wilcox_test_all_UP_vs_DOWN$p.value  # Added Wilcoxon result
  )
  
  # Append the sample results to the all_results data frame
  all_results <<- rbind(all_results, sample_results)
}

# Run the function for each iLINCs data frame
Data_result_list <- lapply(names(iLINCs_selected_samples), process_result_list)

# Convert to data frame to clearly see the data
# Bind rows and include a source column
Data_combined <- bind_rows(Data_result_list, .id = "Source") %>%
  dplyr::distinct(SampleName, .keep_all = TRUE)
print(Data_combined)















