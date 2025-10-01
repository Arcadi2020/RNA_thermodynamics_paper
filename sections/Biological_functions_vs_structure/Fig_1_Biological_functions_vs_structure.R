# Script: Figure 1 – RNA secondary-structure stability across the human transcriptome
# Purpose:
#   Reproduce all plots and summary statistics reported in Fig. 1 (A–F) and
#   Supplementary Figs. 1–2 of the manuscript. The analysis quantifies and visualizes
#   transcript RNA folding stability using length-normalized minimum free energy
#   (Normalized_MFE = MFE / Transcript_length; units: kcal·mol⁻¹·nt⁻¹), and relates
#   stability to RNA biotypes, tissue-specific expression, curated pathways, and
#   evolution. It also exports gene lists for external GSEA.
#
# High-level workflow:
#   1) Load the annotated transcript-level table with per-transcript MFE, heat capacity,
#      lengths, biotypes, and identifiers (Gencode_TX_norm_MFE_HC_table).
#   2) Fig. 1A: (i) Describe global Normalized_MFE distribution for all transcripts and
#      for protein-coding mRNAs; (ii) plot Normalized_MFE vs transcript ID.
#   3) Fig. 1B: Compare Normalized_MFE distributions between protein-coding RNAs, tRNAs,
#      and rRNAs (violin/density plots + KS/Wilcoxon tests).
#   4) Suppl. Fig. 1a–b: Identify protein-coding transcripts in the top (least
#      stable; 95th percentile) and bottom (most stable; 5th percentile) tails of the
#      Normalized_MFE distribution; export gene lists for GO GSEA; plot top terms.
#   5) Fig. 1C: Map Human Protein Atlas (HPA) tissue-specific genes onto the transcript
#      table and compare Normalized_MFE across tissues (box plots with means/95% CI).
#   6) Fig. 1D: Evaluate Normalized_MFE across curated signaling and DNA repair pathways
#      (per-pathway mean ± SD; box plots; reference line for cross-pathway mean).
#   7) **Suppl. Fig. 2:** Validate curated pathway results using MSigDB/GSEA gene sets.
#        - Import KEGG, Reactome, GO, and Hallmark gene sets (JSON format).
#        - Parse into R objects, extract gene symbols, and filter transcript table.
#        - Compute per-pathway mean ± SD of Normalized_MFE.
#        - Plot results (point/error bar plots + box plots) and compare against global mean.
#        - Observation: confirms consistency between curated and external gene sets for
#          key categories (DNA repair, signaling, EMT, stress response).
#   8) Fig. 1E: Contrast Error-prone vs High-fidelity DNA repair mechanisms (polymerases,
#      BER vs NHEJ) using Normalized_MFE distributions; test significance (KS/Wilcoxon).
#   9) Fig. 1F: Phylostratigraphic analysis of Normalized_MFE across evolutionary ages
#      (Trigos et al. PNAS phylostrata dataset), plotting means ± SD and box plots with
#      ANOVA, KS tests, and linear trend analysis.
#
# Required inputs (update paths as needed):
#   - Transcript metrics (CSV):
#       path_to/Gencode_TX_norm_MFE_HC_table.csv
#   - Tissue specificity (TSV from HPA “proteinatlas.tsv.zip”, unzipped to Proteinatlas.txt):
#       /Volumes/.../Proteinatlas.txt
#   - Curated pathway gene sets (R script defining character vectors of HGNC symbols):
#       .../R Scripts/All_gene_sets_AC.R
#   - GENCODE tRNA annotations (GTF) and GRCh38 genome FASTA for sequence extraction:
#       gencode.v47.tRNAs.gtf.gz (or v48) and GRCh38.p14.genome.fa.gz
#   - **MSigDB/GSEA gene sets (JSON):**
#       Downloaded from http://www.gsea-msigdb.org/gsea/downloads.jsp
#       (e.g., c2.cp.v7.5.1.json and pathway-specific JSONs in /GSEA_gene_sets/).
#   - Trigos et al. (PNAS 2017) phylostratigraphy table (CSV).
#
# Outputs:
#   - Human_protein_coding_TX.pdf                    # Fig. 1A scatter
#   - MFE_violin_tRNA_mRNA_rRNA_counts.pdf          # Fig. 1B
#   - GSEA_high_MFE.pdf / GSEA_low_MFE.pdf          # Suppl. Fig. 1
#   - High_MFE.csv / Low_MFE.csv                    # Suppl. Fig. 1
#   - MFE_tissues.pdf                                # Fig. 1C
#   - MFE_curated_signaling_pathways.pdf            # Fig. 1D
#   - Normalized_MFE_GSEA_signaling_pathways.pdf    # Suppl. Fig. 2
#   - Boxplot PDFs for selected GSEA categories      # Suppl. Fig. 2
#   - EP_HF_DNArepair_replication.pdf               # Fig. 1E
#   - MFE_phylostratigraphy.pdf                     # Fig. 1F
#
# Statistical notes:
#   - Distributions compared with KS and Wilcoxon tests.
#   - Values reported as mean ± SD (kcal·mol⁻¹·nt⁻¹).
#   - For Suppl. Fig. 2, each pathway’s mean ± SD compared against global mean line.
#
# Assumptions and conventions:
#   - Normalized_MFE = MFE / transcript length.
#   - External pathways use Ensembl gene IDs → HGNC mapping for joins.
#   - MSigDB JSON parsing requires `rjson` and `rrapply`.
#
# Figure mapping:
#   - Fig. 1A: Scatter/global distributions.
#   - Fig. 1B: Violin plots (protein-coding vs tRNA vs rRNA).
#   - Suppl. Fig. 1: Tail gene lists + GSEA plots.
#   - Fig. 1C: Tissue-specific box plots.
#   - Fig. 1D: Curated pathway box plots.
#   - **Suppl. Fig. 2: GSEA/KEGG/Reactome pathway means ± SD and box plots.**
#   - Fig. 1E: Error-prone vs high-fidelity DNA repair distributions.
#   - Fig. 1F: Phylostratigraphy vs Normalized_MFE (means ± SD, KS, ANOVA).
#






library(dplyr)
library(readr)
library(ggplot2)
library(scales)
library(gtools)
library(tibble)
library(reshape2)
library(tidyverse)
library(ggpubr)
library(tidyr)


# Load Table_ANNOTATED_MFE_heat_capacity ####
Gencode_TX_norm_MFE_HC_table <- read.delim("path_to/Gencode_TX_norm_MFE_HC_table.csv",
                                                         sep = ",", header = T, stringsAsFactors = F)
Table_ANNOTATED_gene_names <- Gencode_TX_norm_MFE_HC_table




# Analysis entire Human annotated transcriptome

# Generate a plot with the distribution of the Normalized_MFE values
ggplot(Table_ANNOTATED_gene_names, aes(x = Normalized_MFE)) +
  geom_histogram(binwidth = 0.001, color = "black") +
  labs(title = "Distribution of the Minimum Free Energy (Normalized_MFE) values of the Human transcriptome",
       x = "Normalized_MFE (kcal/mol)",
       y = "Frequency") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_reverse(limits = c(0, -0.8)) +
  # scale_x_continuous(breaks = seq(-10, 10, by = 0.1)) +
  # scale_y_continuous(breaks = seq(0, 10000, by = 1000)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(legend.position = "none")



# Plot Normalized_MFE vs transcript ID
ggplot(Table_ANNOTATED_gene_names, aes(x = reorder(ensembl_transcript_id_version, Normalized_MFE), y = Normalized_MFE)) +
  theme_classic() +
  geom_point(stat = "identity", colour = "black", size = 2) +
  # theme(axis.text.y = element_text(size=1)) +
  theme(axis.text.y = element_blank()) + # Remove y-axis labels
  theme(axis.title=element_text(size=20)) +
  xlab("Human annotated transcriptome") + 
  ylab("Normalized_MFE (kcal/mol)") +
  ggtitle("Distribution of the Human transcriptome Normalized_MFE") +
  coord_flip()







# Analysis of human protein-coding mRNAs

Table_ANNOTATED_gene_filtered <- dplyr::filter(Table_ANNOTATED_gene_names, transcript_biotype == "protein_coding")

# It is also possible to exclude lncRNAs, ribosomal, Clorf_IDs, MicroRNAs, mitochondrial external_gene_names and snRNAs from the analysis. Transfer RNAs are not present in the initial dataset.
# Table_ANNOTATED_gene_filtered <- Table_ANNOTATED_gene_names %>%
#   dplyr::filter(!is.na(external_gene_name) & external_gene_name != "") %>%
#   dplyr::filter(!startsWith(external_gene_name, "LINC")) %>%
#   dplyr::filter(!startsWith(external_gene_name, "C1orf")) %>%
#   dplyr::filter(!startsWith(external_gene_name, "MIR")) %>%
#   dplyr::filter(!startsWith(external_gene_name, "MTND")) %>%
#   dplyr::filter(!startsWith(external_gene_name, "RN7")) %>%
#   dplyr::filter(!startsWith(external_gene_name, "RNA5")) %>%
#   dplyr::filter(!startsWith(external_gene_name, "RNU1")) %>%
#   dplyr::filter(!startsWith(external_gene_name, "RNU4")) %>%
#   dplyr::filter(!startsWith(external_gene_name, "RNU6"))



# Generate a plot with the distribution of the Normalized_MFE values
ggplot(Table_ANNOTATED_gene_filtered, aes(x = Normalized_MFE)) +
  geom_histogram(binwidth = 0.001, color = "black") +
  labs(title = "Distribution of Normalized Minimum Free Energy Values for Human Protein-Coding mRNAs",
       x = "Normalized_MFE (kcal/mol)",
       y = "Frequency") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_reverse(limits = c(0, -0.8)) +
  # scale_x_continuous(breaks = seq(-10, 10, by = 0.1)) +
  # scale_y_continuous(breaks = seq(0, 10000, by = 1000)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(legend.position = "none")












# Fig.1A ####
# Plot the normalized MFE distribution of the Human protein-coding transcriptome (Normalized_MFE vs transcript ID).
Human_protein_coding_TX <- ggplot(Table_ANNOTATED_gene_filtered, aes(x = reorder(ensembl_transcript_id, Normalized_MFE), y = Normalized_MFE)) +
  theme_classic() +
  geom_point(stat = "identity", colour = "black", size = 1) +
  theme(axis.text.y = element_blank()) + # Remove y-axis labels
  theme(
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.ticks.y = element_blank(), # Remove y-axis tick marks
    axis.title.y = element_blank(), # Remove y-axis title
    axis.title = element_text(size = 20), 
    axis.text.x = element_text(size = 18, color = "black"),
    plot.margin = margin(t = 50, r = 10, b = 10, l = 10), # Extend top margin
    axis.line = element_line(linewidth = 0.2) # Reduce axis thickness
  ) +
  xlab("") + 
  ylab("") +
  ggtitle("") +
  coord_flip()

print(Human_protein_coding_TX)







# Fig.1B visualization secondary structures
# Visualize the RNA secondary structures of selected transcripts
# I have used the RNAfold webServer: http://rna.tbi.univie.ac.at//cgi-bin/RNAWebSuite/RNAfold.cgi



# Identify PROTEIN CODING transcripts with similar length but with the extremes of the Normalized_MFE distribution

# Define a function to get the required rows for each group
# This code generates a data frame which contains several groups of three transcripts with the same length but with the highest, mid and lowest Normalized_MFE values.
select_extremes <- function(df) {
  df <- df %>%
    arrange(Normalized_MFE) %>%
    mutate(row_index = row_number())
  
  mid_index <- ceiling(nrow(df) / 2)
  
  df %>%
    filter(row_index == 1 | row_index == mid_index | row_index == nrow(df)) %>%
    dplyr::select(-row_index)
}

# Apply the function to each group and combine the results
Result_selected_TX_MFE <- Table_ANNOTATED_gene_names %>%
  group_by(Transcript_length) %>%
  group_modify(~ select_extremes(.x)) %>%
  ungroup() %>%
  dplyr::filter(transcript_biotype == "protein_coding")

# View the result
Result_selected_TX_MFE_df <- as.data.frame(Result_selected_TX_MFE)

# Filter Transcript_length values that are present at least 3 times
# This will select transcripts that have more common lengths
# Result_selected_TX_MFE_df <- Result_selected_TX_MFE_df %>%
#   group_by(Transcript_length) %>%
#   filter(n() >= 3) %>%
#   ungroup()

# Look at transcriots with similar length but very different Normalized_MFE values
# Identify the average transcript length (this will broaden the spectrum of candidate transcripts)

# The mean transcript lengths in the adult cortex is 2,530 bp. These measurements were based on full-length transcript sequencing.
# Data have been taken from the paper:
# Full-length transcript sequencing of human and mouse cerebral cortex identifies widespread isoform diversity and alternative splicing, Szi Kay Leung et al. Cell Rep. 2021


# I select the transcripts with a length of 2,530 bp
Result_selected_TX_MFE_df_2530 <- Result_selected_TX_MFE_df %>%
  dplyr::filter(Transcript_length >= 2520 & Transcript_length <= 2540)


# The two transcripts with a similar length but very different Normalized_MFE values, identified by the code are:
# ENST00000413702.5
# ENST00000216241.14


# Get the fasta sequences of both transcripts
# Load the biomaRt package
library(biomaRt)

# Connect to the Ensembl BioMart database
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Define the transcript IDs
transcript_ids <- c("ENST00000413702.5", "ENST00000216241.14")

# Retrieve the FASTA sequences
fasta_sequences <- getSequence(
  id = transcript_ids,
  type = "ensembl_transcript_id_version",
  seqType = "cdna", # Use "cdna" for transcript sequences
  mart = ensembl
)



# Get the sequence of the transcripts
Max_target_index <- dplyr::filter(fasta_sequences, ensembl_transcript_id_version == "ENST00000413702.5") %>%
  dplyr::select(cdna) %>%
  dplyr::pull() # Paste this sequence into the RNAfold webServer as explained here above

Min_target_index <- dplyr::filter(fasta_sequences, ensembl_transcript_id_version == "ENST00000216241.14") %>%
  dplyr::select(cdna) %>%
  dplyr::pull() # Paste this sequence into the RNAfold webServer as explained here above


















# Supplementary Fig.1a ####
# Select the upper and lower quantiles of the distribution of the human protein coding RNAs
filtered_data <- na.omit(Table_ANNOTATED_gene_filtered$Normalized_MFE)

# Calculate the lower and upper quantiles
lower_quantile <- quantile(filtered_data, probs = 0.05)
upper_quantile <- quantile(filtered_data, probs = 0.95)

# Print the quantiles
print(paste("Lower Quantile (5th percentile):", lower_quantile))
print(paste("Upper Quantile (95th percentile):", upper_quantile))

# Subset the data to find the extremes
extremes <- filtered_data[filtered_data <= lower_quantile | filtered_data >= upper_quantile]
extremes_df <- as.data.frame(extremes)

# Filter the genes above the upper quantile. There are the transcripts with very unstable secondary structures (less negative MFE values)
High_MFE <- Table_ANNOTATED_gene_filtered %>%
  filter(Normalized_MFE > upper_quantile) %>%
  distinct(external_gene_name, .keep_all = TRUE)



# The gene list has been uploaded on https://geneontology.org/ for the Gene Set Enrichment Analysis (GSEA) in the  biological process database
# The output has been exported as a table (txt format) and then converted into csv with excel
# Read the CSV file
GSEA_High_MFE <- read.csv("path_to/High_MFE_genes.csv", header = TRUE, stringsAsFactors = FALSE)
# Remove the first 5 rows
GSEA_High_MFE <- GSEA_High_MFE[-c(1:5), ]
# Replace the headers with the first row
colnames(GSEA_High_MFE) <- GSEA_High_MFE[1, ]
GSEA_High_MFE <- GSEA_High_MFE[-1, ]  # Remove the first row now that it's used as headers
# Reset the row numbers (optional, for clean numbering)
rownames(GSEA_High_MFE) <- NULL





GSEA_High_MFE_filtered <- dplyr::filter(GSEA_High_MFE, `upload_1 (fold Enrichment)` > 1.6) %>%
  dplyr::slice(1:25)

# Use gsub to remove the pattern in the specified column
GSEA_High_MFE_filtered$`GO biological process complete` <- gsub(
  "\\s*\\(GO:\\d+\\)",  # Pattern to match (GO:XXXXX), including the space before it
  "",                   # Replacement: an empty string
  GSEA_High_MFE_filtered$`GO biological process complete`) # Target column


# Reduce the length of some of the GO terms
# GSEA_High_MFE_filtered$GO.biological.process.complete <- ifelse(GSEA_High_MFE_filtered$GO.biological.process.complete == "regulation of DNA damage response, signal transduction by p53 class mediator ", "regulation of DNA damage response ", GSEA_High_MFE_filtered$GO.biological.process.complete)
# GSEA_High_MFE_filtered$GO.biological.process.complete <- ifelse(GSEA_High_MFE_filtered$GO.biological.process.complete == "positive regulation of telomere maintenance via telomere lengthening ", "positive regulation of telomere maintenance ", GSEA_High_MFE_filtered$GO.biological.process.complete)

# Convert `upload_1 (fold Enrichment)` to numeric (force conversion and handle NA values)
GSEA_High_MFE_filtered$`upload_1 (fold Enrichment)` <- as.numeric(GSEA_High_MFE_filtered$`upload_1 (fold Enrichment)`)

# Check for NAs in `upload_1 (fold Enrichment)` and handle them (optional, e.g., remove rows with NAs)
GSEA_High_MFE_filtered <- GSEA_High_MFE_filtered %>%
  filter(!is.na(`upload_1 (fold Enrichment)`))

# Generate bar plot
GSEA_high_MFE <- ggplot(
  GSEA_High_MFE_filtered, 
  aes(x = -log10(`upload_1 (fold Enrichment)`), y = reorder(`GO biological process complete`, -log10(`upload_1 (fold Enrichment)`)))) +
  geom_bar(stat = "identity", fill = "steelblue") +  # Bar color
  scale_x_reverse() +  # Reverse the x-axis
  scale_y_discrete(limits = rev) +  # Reverse the y-axis order
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),  # Rotate x-axis labels
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 8)
  ) +
  xlab("Fold Enrichment (-log10)") +
  ylab("") +
  ggtitle("GSEA of the Most Unstable RNAs")

# Print the plot
print(GSEA_high_MFE)

# Save the combined plot
ggsave("GSEA_high_MFE.pdf", plot = GSEA_high_MFE, width = 7, height = 6, path = "/Users/z3514253/Desktop")










# Supplementary Fig.1b ####
# Filter the genes below the lower quantile. There are the transcripts with very stable secondary structures (highly negative MFE values)
Low_MFE <- Table_ANNOTATED_gene_filtered %>%
  filter(Normalized_MFE < lower_quantile) %>%
  distinct(external_gene_name, .keep_all = TRUE)

# Write the filtered data to a CSV file
write.csv(Low_MFE, file = "path_to/Low_MFE.csv", row.names = F)


# The gene list has been uploaded on https://geneontology.org/ for the Gene Set Enrichment Analysis (GSEA)
# The output has been exported as a table (txt format) and then converted into csv with excel
# Read the CSV file
GSEA_Low_MFE <- read.csv("path_to/Low_MFE_genes.csv", header = TRUE, stringsAsFactors = FALSE)
# Remove the first 5 rows
GSEA_Low_MFE <- GSEA_Low_MFE[-c(1:5), ]
# Replace the headers with the first row
colnames(GSEA_Low_MFE) <- GSEA_Low_MFE[1, ]
GSEA_Low_MFE <- GSEA_Low_MFE[-1, ]  # Remove the first row now that it's used as headers
# Reset the row numbers (optional, for clean numbering)
rownames(GSEA_Low_MFE) <- NULL

# Set a threshold to the fold enrichment values
GSEA_Low_MFE_filtered <- dplyr::filter(GSEA_Low_MFE, `upload_1 (fold Enrichment)` > 1.9)


# Use gsub to remove the pattern in the specified column
GSEA_Low_MFE_filtered$`GO biological process complete` <- gsub(
  "\\s*\\(GO:\\d+\\)",  # Pattern to match (GO:XXXXX), including the space before it
  "",                   # Replacement: an empty string
  GSEA_Low_MFE_filtered$`GO biological process complete`) # Target column


# Convert `upload_1 (fold Enrichment)` to numeric (force conversion and handle NA values)
GSEA_Low_MFE_filtered$`upload_1 (fold Enrichment)` <- as.numeric(GSEA_Low_MFE_filtered$`upload_1 (fold Enrichment)`)

# Check for NAs in `upload_1 (fold Enrichment)` and handle them (optional, e.g., remove rows with NAs)
GSEA_Low_MFE_filtered <- GSEA_Low_MFE_filtered %>%
  filter(!is.na(`upload_1 (fold Enrichment)`))

# Generate bar plot
GSEA_low_MFE <- ggplot(
  GSEA_Low_MFE_filtered, 
  aes(x = -log10(`upload_1 (fold Enrichment)`), y = reorder(`GO biological process complete`, -log10(`upload_1 (fold Enrichment)`)))) +
  geom_bar(stat = "identity", fill = "steelblue") +  # Bar color
  scale_x_reverse() +  # Reverse the x-axis
  scale_y_discrete(limits = rev) +  # Reverse the y-axis order
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),  # Rotate x-axis labels
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 8)
  ) +
  xlab("Fold Enrichment (-log10)") +
  ylab("") +
  ggtitle("GSEA of the Most Stable RNAs")

# Print the plot
print(GSEA_low_MFE)


























# Fig.1B ####
# Analysis MFE tRNAs and ribosomal RNAs vs protein coding RNAs ####


library(Biostrings)
library(tibble)
library(dplyr)
library(stringr)
library(LncFinder)
library(tidyr)
library(seqinr)
library(devtools)
library(tidyverse)
library(FastaUtils)
library(biomaRt)
library(cbioportalR)
library(R.utils)
library(rtracklayer)




Gencode_TX_norm_MFE_HC_table <- read.delim("path_to/Gencode_TX_norm_MFE_HC_table.csv",
                                           sep = ",", header = T, stringsAsFactors = F)
Table_ANNOTATED_gene_names <- Gencode_TX_norm_MFE_HC_table
# This file contains:
# Protein coding transcripts
# lncRNAs
# ribosomal RNAs



# Select the ribosomal RNAs (rRNAs)
Gencode_rRNAs <- Table_ANNOTATED_gene_names %>%
  dplyr::filter(transcript_biotype == "rRNA")

# Select the lncRNAs
Gencode_lncRNAs <- Table_ANNOTATED_gene_names %>%
  dplyr::filter(transcript_biotype == "lncRNA")

# Select the protein coding transcripts
Gencode_Protein_Coding <- Table_ANNOTATED_gene_names %>%
  dplyr::filter(transcript_biotype == "protein_coding")

# Select the tRNAs from the Gencode transcriptome
# The Gencode tRNAs are located in a different GTF file, so I need to load them separately
# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.tRNAs.gtf.gz

# Load tRNA annotations (filtering unwanted types)
tRNAs <- readGFF("path_to/gencode.v47.tRNAs.gtf.gz") %>%
  filter(!transcript_type %in% c("Pseudo_tRNA", "Undet_tRNA", "Sup_tRNA", "SeC_tRNA", "SeC(e)_tRNA"))

# Access the FASTA files
# The tRNAs GTF entries have transcript_id values like 20658, whereas Gencode_transcripts_FASTA uses canonical Ensembl transcript IDs like ENST00000xxxxx. This mismatch happens because:
# The predicted tRNAs file (gencode.v48.tRNAs.gtf) does not use ENST IDs; it uses numerical IDs (internal ENSEMBL IDs).
# The transcript FASTA file (gencode.v47.transcripts.fa) contains only standard ENST transcripts, not the predicted tRNAs.


# Use the GTF coordinates to extract tRNAs
# Download the genome FASTA file for Gencode v47
# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/GRCh38.p14.genome.fa.gz

# Load the file
genome_path <- "path_to/GRCh38.p14.genome.fa.gz"

# Read the compressed genome FASTA
genome <- Biostrings::readDNAStringSet(genome_path)


# Clean chromosome names (remove any descriptions after space)
names(genome) <- sub(" .*", "", names(genome))

# 3. Extract sequences manually (subseq + reverseComplement for negative strand)
tRNA_sequences <- DNAStringSet(mapply(function(chr, start, end, strand) {
  seq <- subseq(genome[[chr]], start = start, end = end)
  if (strand == "-") seq <- reverseComplement(seq)
  seq
}, tRNAs$seqid, tRNAs$start, tRNAs$end, tRNAs$strand))

# 4. Assign transcript IDs as names
names(tRNA_sequences) <- tRNAs$transcript_id

# 5. Save to FASTA with transcript_id as header
output_path <- "path_to/tRNAs_sequences_v47.fa"
writeXStringSet(tRNA_sequences, filepath = output_path)


library(seqinr)
library(LncFinder)
library(dplyr)
library(tibble)
library(stringr)

# 1. Read the FASTA file of extracted tRNAs (preserve IDs)
Human_t_RNAs <- seqinr::read.fasta(
  file = output_path,
  seqtype = "DNA",
  as.string = TRUE,
  set.attributes = FALSE
)

# 2. Predict secondary structure and MFE using RNAfold
Secondary_structure <- run_RNAfold(
  Human_t_RNAs,
  RNAfold.path = "RNAfold",        # Ensure RNAfold is in your PATH
  parallel.cores = 3
)

# 3. Predict coding potential using LncFinder (with SS features)
Prediction_on_Secondary_structure_t_RNAs <- lnc_finder(
  Secondary_structure,
  SS.features = TRUE,
  format = "SS",
  frequencies.file = "human",
  svm.model = "human",
  parallel.cores = 3
)

# 4. Convert prediction results to a data frame
Prediction_Gencode_tRNAs_df <- as.data.frame(Prediction_on_Secondary_structure_t_RNAs)
Prediction_Gencode_tRNAs_df <- tibble::rownames_to_column(Prediction_Gencode_tRNAs_df, var = "ID")

# Remove leading 'X' if present in ID
Prediction_Gencode_tRNAs_df$ID <- sub("^X", "", Prediction_Gencode_tRNAs_df$ID)


# 5. Calculate length of each tRNA sequence
tRNA_length_df <- data.frame(
  ID = names(tRNA_sequences),
  Length = width(tRNA_sequences)
)

# 6. Add transcript_name from GTF
tRNA_name_df <- tRNAs %>%
  dplyr::select(transcript_id, transcript_name) %>%
  dplyr::rename(ID = transcript_id, tRNA_name = transcript_name)

# 7. Merge predictions with length and name
Prediction_Gencode_tRNAs_df <- Prediction_Gencode_tRNAs_df %>%
  dplyr::left_join(tRNA_length_df, by = "ID") %>%
  dplyr::left_join(tRNA_name_df, by = "ID")

# Calculate Normalized_MFE
Gencode_tRNAs <- Prediction_Gencode_tRNAs_df %>%
  dplyr::mutate(Normalized_MFE = MFE / Length)









# tRNAs (from Prediction_Gencode_tRNAs_df)
tRNA_MFE_labeled <- Gencode_tRNAs %>%
  dplyr::select(Normalized_MFE) %>%
  dplyr::mutate(Group = "tRNA")

# rRNAs (from main annotated table)
rRNA_MFE_labeled <- Gencode_rRNAs %>%
  dplyr::select(Normalized_MFE) %>%
  dplyr::mutate(Group = "rRNA")

# Protein-coding RNAs
ProteinCoding_MFE_labeled <- Table_ANNOTATED_gene_names %>%
  dplyr::filter(transcript_biotype == "protein_coding") %>%
  dplyr::select(Normalized_MFE) %>%
  dplyr::mutate(Group = "Protein_coding")

# lncRNAs
lncRNA_MFE_labeled <- Table_ANNOTATED_gene_names %>%
  dplyr::filter(transcript_biotype == "lncRNA") %>%
  dplyr::select(Normalized_MFE) %>%
  dplyr::mutate(Group = "lncRNA")



# Generate the violin plot with the protein-coding, tRNAs and rRNAs
combined_MFE <- bind_rows(
  tRNA_MFE_labeled,
  rRNA_MFE_labeled,
  ProteinCoding_MFE_labeled
) %>%
  mutate(Group = factor(Group, levels = c("Protein_coding", "tRNA", "rRNA")))

# Violin plot
# Count transcripts per group
counts <- combined_MFE %>%
  group_by(Group) %>%
  summarise(n = n(), .groups = "drop")

# Create new labels with n=
group_labels <- paste0(levels(combined_MFE$Group), "\n(n=", counts$n[match(levels(combined_MFE$Group), counts$Group)], ")")

# Violin plot with counts
MFE_violin <- ggplot(combined_MFE, aes(x = Group, y = Normalized_MFE, fill = Group)) +
  geom_violin(trim = FALSE, colour = "black", linewidth = 0.3, alpha = 0.70) +
  geom_boxplot(width = 0.15, outlier.shape = NA,
               fill  = "white", colour = "black", linewidth = 0.3) +
  scale_fill_manual(values = rep("lightgrey", 3)) +
  scale_x_discrete(labels = group_labels) +        # <-- Add labels with counts
  theme_classic() +
  labs(x = NULL, y = "Normalized MFE (kcal·mol⁻¹·nt⁻¹)") +
  theme(
    legend.position = "none",
    axis.title.y    = element_text(size = 16),
    axis.text       = element_text(size = 14)
  )

print(MFE_violin)







# # Extract vectors without NAs
tRNA_vals    <- na.omit(tRNA_MFE_labeled$Normalized_MFE)
rRNA_vals    <- na.omit(rRNA_MFE_labeled$Normalized_MFE)
protein_vals <- na.omit(ProteinCoding_MFE_labeled$Normalized_MFE)

# Means
cat("Means:\n",
    "tRNA:", mean(tRNA_vals), "\n",
    "rRNA:", mean(rRNA_vals), "\n",
    "Protein-coding:", mean(protein_vals), "\n\n")

# KS tests
cat("Kolmogorov–Smirnov tests:\n")
print(ks.test(tRNA_vals, rRNA_vals))
print(ks.test(tRNA_vals, protein_vals))
print(ks.test(rRNA_vals, protein_vals))

# Wilcoxon rank-sum tests
cat("\nWilcoxon rank-sum tests:\n")
print(wilcox.test(tRNA_vals, rRNA_vals))
print(wilcox.test(tRNA_vals, protein_vals))
print(wilcox.test(rRNA_vals, protein_vals))






# Density plot with median lines (Protein-coding → tRNA → rRNA)
ggplot(combined_MFE, aes(x = Normalized_MFE, fill = Group, colour = Group)) +
  geom_density(aes(y = ..scaled..), alpha = 0.5) +  # scaled densities
  geom_vline(
    data = combined_MFE %>%
      group_by(Group) %>%
      summarise(median = median(Normalized_MFE, na.rm = TRUE), .groups = "drop"),
    aes(xintercept = median, colour = Group),
    linetype = "dashed"
  ) +
  scale_fill_manual(values = c(
    "Protein_coding" = "#1f78b4",  # blue
    "tRNA"           = "#e31a1c",  # red
    "rRNA"           = "#ff7f00"   # orange
  )) +
  scale_colour_manual(values = c(
    "Protein_coding" = "#1f78b4",
    "tRNA"           = "#e31a1c",
    "rRNA"           = "#ff7f00"
  )) +
  theme_classic() +
  labs(x = "Normalized MFE (kcal·mol⁻¹·nt⁻¹)", y = "Scaled density") +
  theme(
    axis.title   = element_text(size = 16),
    axis.text    = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = "top"
  )


















# Fig.1C ####
# Analysis tissue-specific genes

# Load the protein atals data which provide information about the tissue-specific gene expression.
# Load the tissue specific gene expression data
# These data have been downloaded from: https://www.proteinatlas.org/about/download
# The file name is: Data from the Human Protein Atlas in tab-separated format (proteinatlas.tsv.zip)
Proteinatlas_data <- read.delim("path_to/Proteinatlas.txt") %>%
  dplyr::rename(ensembl_gene_id = Ensembl)




library(stringr)
Skeletal_muscle_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "skeletal muscle")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id") %>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Select the genes which are specific for the retina
Retina_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "retina")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id") %>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Heart specific genes
Heart_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "heart")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id") %>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Adrenal gland specific genes
Adrenal_gland_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "adrenal gland")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id") %>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Parathyroid gland specific genes
Parathyroid_gland_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "parathyroid gland")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id") %>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Thyroid gland specific genes
Thyroid_gland_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "thyroid gland")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id") %>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Pituitary gland specific genes
Pituitary_gland_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "pituitary gland")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id") %>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Lung specific genes
Lung_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "lung")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id") %>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Bone marrow specific genes
Bone_marrow_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "bone marrow")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id") %>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Lymphoid tissue specific genes
Lymphoid_tissue_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "lymphoid tissue")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id") %>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Liver specific genes
Liver_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "liver")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id") %>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Gallbladder specific genes
Gallbladder_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "gallbladder")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id") %>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Testis specific genes
Testis_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "testis")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id") %>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Epididymis specific genes
Epididymis_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "epididymis")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id") %>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Prostate specific genes
Prostate_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "prostate")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id") %>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Seminal vesicle specific genes
Seminal_vesicle_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "seminal vesicle")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id") %>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Adipose tissue specific genes
Adipose_tissue_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "adipose tissue")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id")%>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Brain specific genes
Brain_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "brain")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id")%>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Choroid plexus specific genes
Choroid_plexus_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "choroid plexus")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id")%>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Salivary gland specific genes
Salivary_gland_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "salivary gland")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id")%>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Esophagus specific genes
Esophagus_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "esophagus")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id")%>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Tongue specific genes
Tongue_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "tongue")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id")%>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Stomach specific genes
Stomach_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "stomach")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id")%>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Intestine specific genes
Intestinal_tissue_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "intestine")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id")%>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Pancreas specific genes
Pancreas_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "pancreas")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id")%>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Kidney specific genes
Kidney_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "kidney")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id")%>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Urinary bladder specific genes
Urinary_bladder_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "urinary bladder")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id")%>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Breast specific genes
Breast_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "breast")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id")%>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Vagina specific genes
Vagina_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "vagina")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id")%>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Cervix specific genes
Cervix_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "cervix")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id")%>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Endometrium specific genes
Endometrium_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "endometrium")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id")%>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Fallopian tube specific genes
Fallopian_tube_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "fallopian tube")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id")%>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Ovary specific genes
Ovary_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "ovary")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id")%>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Placenta specific genes
Placenta_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "placenta")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id")%>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)

# Skin specific genes
Skin_specific_genes <- Proteinatlas_data %>%
  dplyr::filter(str_detect(RNA.tissue.specific.nTPM, "skin")) %>%
  dplyr::select(ensembl_gene_id, RNA.tissue.specific.nTPM) %>%
  dplyr::left_join(Table_ANNOTATED_gene_names, by = "ensembl_gene_id")%>%
  dplyr::select(ensembl_transcript_id_version, Normalized_MFE)



# Generate a box plot with the data of each tissue on the x axis and the Normalized_MFE values on the y axis

library(tidyverse)

# Define the function to add the columns with mean and confidence intervals to each dataframe
add_mean_ci_columns <- function(df, confidence_level = 0.95) {
  if ("Normalized_MFE" %in% colnames(df) && !all(is.na(df$Normalized_MFE))) {
    # Calculate mean
    df$mean_value <- mean(df$Normalized_MFE, na.rm = TRUE)
    
    # Calculate standard error
    se <- sd(df$Normalized_MFE, na.rm = TRUE) / sqrt(length(df$Normalized_MFE))
    
    # Calculate margin of error
    moe <- qnorm((1 + confidence_level) / 2) * se
    
    # Calculate confidence interval
    df$lower_ci <- df$mean_value - moe
    df$upper_ci <- df$mean_value + moe
  } else {
    df$mean_value <- NA
    df$lower_ci <- NA
    df$upper_ci <- NA
  }
  return(df)
}



# List of tissue-specific gene data frames
tissue_specific_genes_list <- list(
  Skeletal_muscle = Skeletal_muscle_specific_genes,
  Retina = Retina_specific_genes,
  Heart = Heart_specific_genes,
  Adrenal_gland = Adrenal_gland_specific_genes,
  Parathyroid_gland = Parathyroid_gland_specific_genes,
  Thyroid_gland = Thyroid_gland_specific_genes,
  Pituitary_gland = Pituitary_gland_specific_genes,
  Lung = Lung_specific_genes,
  Bone_marrow = Bone_marrow_specific_genes,
  Lymphoid_tissue = Lymphoid_tissue_specific_genes,
  Liver = Liver_specific_genes,
  Gallbladder = Gallbladder_specific_genes,
  Testis = Testis_specific_genes,
  Epididymis = Epididymis_specific_genes,
  Prostate = Prostate_specific_genes,
  Seminal_vesicle = Seminal_vesicle_specific_genes,
  Adipose_tissue = Adipose_tissue_specific_genes,
  Brain = Brain_specific_genes,
  Choroid_plexus = Choroid_plexus_specific_genes,
  Salivary_gland = Salivary_gland_specific_genes,
  Esophagus = Esophagus_specific_genes,
  Tongue = Tongue_specific_genes,
  Stomach = Stomach_specific_genes,
  Intestinal_tissue = Intestinal_tissue_specific_genes,
  Pancreas = Pancreas_specific_genes,
  Kidney = Kidney_specific_genes,
  Urinary_bladder = Urinary_bladder_specific_genes,
  Breast = Breast_specific_genes,
  Vagina = Vagina_specific_genes,
  Cervix = Cervix_specific_genes,
  Endometrium = Endometrium_specific_genes,
  Fallopian_tube = Fallopian_tube_specific_genes,
  Ovary = Ovary_specific_genes,
  Placenta = Placenta_specific_genes,
  Skin = Skin_specific_genes
)

# Function to add tissue column to each dataframe
add_tissue_column <- function(tissue_df, tissue_name) {
  tissue_df$tissue <- tissue_name
  return(tissue_df)
}

# Add tissue column and mean column to each dataframe
for (tissue_name in names(tissue_specific_genes_list)) {
  tissue_specific_genes_list[[tissue_name]] <- tissue_specific_genes_list[[tissue_name]] %>%
    add_tissue_column(tissue_name) %>%
    add_mean_ci_columns()
}

# Combine all dataframes into one dataframe
all_tissue_genes <- bind_rows(tissue_specific_genes_list)





# Plot the data using a boxplot
# Calculate the average of all mean_values
Mean_MFE <- all_tissue_genes %>% dplyr:: summarise(Mean_MFE_all_tissues = mean(mean_value)) 

# Reorder samples based on mean_value
all_tissue_genes <- all_tissue_genes[order(all_tissue_genes$mean_value), ]

# Create a vector to specify which samples to fill
fill_vector <- ifelse(all_tissue_genes$tissue %in% c("Skeletal_muscle_specific_genes", 
                                                     "Liver_specific_genes", 
                                                     "Kidney_specific_genes",
                                                     "Pituitary_gland_specific_genes",
                                                     "Heart_specific_genes",
                                                     "Adipose_tissue_specific_genes",
                                                     "Brain_specific_genes"), as.factor(all_tissue_genes$tissue), "Other")

MFE_tissues <- ggplot(all_tissue_genes, aes(x = reorder(tissue, -mean_value), y = Normalized_MFE, 
                               #fill = fill_vector # To fill the selected samples
                               )) +
geom_boxplot() +
stat_summary(fun=mean, geom="point", shape=19, size=1, color="black") + # Add mean points
theme_classic() +
geom_hline(yintercept = Mean_MFE$Mean_MFE_all_tissues, color = "black") +
labs(x = "", y = "") +
theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
          axis.text.y = element_blank(),,
          axis.title.x = element_text(size = 18)) +
scale_fill_manual(values = c("Skeletal_muscle_specific_genes" = "lightgrey", 
                                 "Liver_specific_genes" = "lightgrey", 
                                 "Kidney_specific_genes" = "lightgrey", 
                                 "Pituitary_gland_specific_genes" = "lightgrey",
                                 "Heart_specific_genes" = "lightgrey",
                                 "Adipose_tissue_specific_genes" = "lightgrey",
                                 "Brain_specific_genes" = "lightgrey",
                                 "Other" = "white"), guide = FALSE) +
ggtitle("") +
coord_flip()

print(MFE_tissues)
  
































# Fig.1D ####
# Analysis MFE curated signaling pathways
# More in general, I wonder whether the MFE had any impact on the evolution of the canonical signaling pathways.
# In other words, do specific signaling pathways have a total high or low MFE which is consistent with the rate of proliferation, during which the specific signaling pathway is active? 
# CALCULATE THE MFE OF CANONICAL SIGNALING PATHWAYS

# This is the script which contains the curated gene sets
source("path_to/All_gene_sets_AC.R")





# These are the commands to filter the 'Table_ANNOTATED_gene_names' file (which contains the Normalized MFE of each transcript) with each curated gene set

# Generate the function that filters each gene set in the target data frame
# This function will be used in the following commands
filtering <- function(x) {
  Filtering<-dplyr::filter(Table_ANNOTATED_gene_names, external_gene_name %in% x)
  return(Filtering)
}


# Generate the list containing all the gene sets to be used for the filtering.
CURATED_PATHWAYS.input <- list(MTOR_signaling, Cell_migration, Cell_motility, Gene_set_cell_adhesion, Chemokines_and_chemotaxis, Gene_set_Splicing_Factors, Homologous_recombination, Nucleotide_excision_repair, Non_hologous_end_joining, Fanconi_anemia, Base_excision_repair, Mismatch_repair, High_fidelity_polymerases, Error_prone_DNA_polymerases, ABC_SLC_transporters, calcium_channel_activity, TNF_signaling_pathway, Chaperone_mediated_protein_folding, TGFb_BMP_signaling_pathway, NF_kB_Target_genes, TGFb_signaling_targets, Notch_signaling_pathway, Hedgehog_Signaling_Pathway, WNT_signaling_pathway, WNT_signaling_targets, AMPK_signaling, cAMP_calcium_signaling_pathway, EGF_PDGF_signaling_pathway, Hypoxia_signaling_pathway, IL6_STAT3_signaling_pathway, Insulin_signaling_pathway, JAK_STAT_signaling_pathways, MAP_Kinase_signaling_pathways, p53_signaling_pathway, PI3K_AKT_signaling_pathway, TNF_signaling_pathway, Toll_Like_receptor_signaling_pathway, VEGF_signaling, Cellular_Stress_Responses, Antibacterial_response, Antifungal_response, Antiviral_response, Apoptosis, Autophagy, Cancer_drug_resistance, JNK_signaling_Vivanco, Pro_survival, Chromatin_remodelling)

# This creates the names that will be used to for the final data frames of each file resulting from the filtering of the target data frame with each gene set
CURATED_PATHWAYS.output <- c("MTOR_signaling", "Cell_migration", "Cell_motility", "Gene_set_cell_adhesion", "Chemokines_and_chemotaxis", "Gene_set_Splicing_Factors", "Homologous_recombination", "Nucleotide_excision_repair", "Non_hologous_end_joining", "Fanconi_anemia", "Base_excision_repair", "Mismatch_repair", "High_fidelity_polymerases", "Error_prone_DNA_polymerases", "ABC_SLC_transporters", "calcium_channel_activity", "TNF_signaling_pathway", "Chaperone_mediated_protein_folding", "TGFb_BMP_signaling_pathway", "NF_kB_Target_genes", "TGFb_signaling_targets", "Notch_signaling_pathway", "Hedgehog_Signaling_Pathway", "WNT_signaling_pathway", "WNT_signaling_targets", "AMPK_signaling", "cAMP_calcium_signaling_pathway", "EGF_PDGF_signaling_pathway", "Hypoxia_signaling_pathway", "IL6_STAT3_signaling_pathway", "Insulin_signaling_pathway", "JAK_STAT_signaling_pathways", "MAP_Kinase_signaling_pathways", "p53_signaling_pathway", "PI3K_AKT_signaling_pathway", "TNF_signaling_pathway", "Toll_Like_receptor_signaling_pathway", "VEGF_signaling", "Cellular_Stress_Responses", "Antibacterial_response", "Antifungal_response", "Antiviral_response", "Apoptosis", "Autophagy", "Cancer_drug_resistance", "JNK_signaling_Vivanco", "Pro_survival", "Chromatin_remodelling")

# Assign the names of the gene sets to each entry in the list
names(CURATED_PATHWAYS.input) <- CURATED_PATHWAYS.output

# Apply the filtering function to create a list of data.frames
Curated_test_df <- lapply(CURATED_PATHWAYS.input, function(x) filtering(x))

# creating a vector with desired names of filtered data.frames
Curated_dfs.output <- c("df_MTOR_signaling", "df_Cell_migration", "df_Cell_motility", "df_Gene_set_cell_adhesion", "df_Chemokines_and_chemotaxis", "df_Gene_set_Splicing_Factors", "df_Homologous_recombination", "df_Nucleotide_excision_repair", "df_Non_hologous_end_joining", "df_Fanconi_anemia", "df_Base_excision_repair", "df_Mismatch_repair", "df_High_fidelity_polymerases", "df_Error_prone_DNA_polymerases", "df_ABC_SLC_transporters", "df_calcium_channel_activity", "df_TNF_signaling_pathway", "df_Chaperone_mediated_protein_folding", "df_TGFb_BMP_signaling_pathway", "df_NF_kB_Target_genes", "df_TGFb_signaling_targets", "df_Notch_signaling_pathway", "df_Hedgehog_Signaling_Pathway", "df_WNT_signaling_pathway", "df_WNT_signaling_targets", "df_AMPK_signaling", "df_cAMP_calcium_signaling_pathway", "df_EGF_PDGF_signaling_pathway", "df_Hypoxia_signaling_pathway", "df_IL6_STAT3_signaling_pathway", "df_Insulin_signaling_pathway", "df_JAK_STAT_signaling_pathways", "df_MAP_Kinase_signaling_pathways", "df_p53_signaling_pathway", "df_PI3K_AKT_signaling_pathway", "df_TNF_signaling_pathway", "df_Toll_Like_receptor_signaling_pathway", "df_VEGF_signaling", "df_Cellular_Stress_Responses", "df_Antibacterial_response", "df_Antifungal_response", "df_Antiviral_response", "df_Apoptosis", "df_Autophagy", "df_Cancer_drug_resistance", "df_JNK_signaling_Vivanco", "df_Pro_survival", "df_Chromatin_remodelling")

# loop to assign data.frames to the desired variable name
for (i in (seq_along(Curated_dfs.output))) {
  do.call("<-",list(Curated_dfs.output[i], Curated_test_df[[i]]))
}

# Now simply type a gene set (for example: 'df_Cell_migration' to have the full table of the filtered data)
# Each gene set is a data frame that can be used to generate the plots.





# Generate new functions 'mean' and 'sd' and loop them through the 'Normalized_MFE' columns of each data frame generated by the automatic filtering here above

# Generate the function for to calculate the MEAN of the normalized_MFE of each gene set resulting from the fuiltering here above.
Mean_curated <- function(x) {
  mean_curated<-mean(x$Normalized_MFE)
  return(mean_curated)
}

# Create an input of the data frames generated here above that I will use to calculate the mean of a specific column
Curated_input_dfs <- list(df_MTOR_signaling, df_Cell_migration, df_Cell_motility, df_Gene_set_cell_adhesion, df_Chemokines_and_chemotaxis, df_Gene_set_Splicing_Factors, df_Homologous_recombination, df_Nucleotide_excision_repair, df_Non_hologous_end_joining, df_Fanconi_anemia, df_Base_excision_repair, df_Mismatch_repair, df_High_fidelity_polymerases, df_Error_prone_DNA_polymerases, df_ABC_SLC_transporters, df_calcium_channel_activity, df_TNF_signaling_pathway, df_Chaperone_mediated_protein_folding, df_TGFb_BMP_signaling_pathway, df_NF_kB_Target_genes, df_TGFb_signaling_targets, df_Notch_signaling_pathway, df_Hedgehog_Signaling_Pathway, df_WNT_signaling_pathway, df_WNT_signaling_targets, df_AMPK_signaling, df_cAMP_calcium_signaling_pathway, df_EGF_PDGF_signaling_pathway, df_Hypoxia_signaling_pathway, df_IL6_STAT3_signaling_pathway, df_Insulin_signaling_pathway, df_JAK_STAT_signaling_pathways, df_MAP_Kinase_signaling_pathways, df_p53_signaling_pathway, df_PI3K_AKT_signaling_pathway, df_TNF_signaling_pathway, df_Toll_Like_receptor_signaling_pathway, df_VEGF_signaling, df_Cellular_Stress_Responses, df_Antibacterial_response, df_Antifungal_response, df_Antiviral_response, df_Apoptosis, df_Autophagy, df_Cancer_drug_resistance, df_JNK_signaling_Vivanco, df_Pro_survival, df_Chromatin_remodelling)

# Apply the mean function to create a list of data.frames
Curated_mean_test_df <- lapply(Curated_input_dfs, function(x) Mean_curated(x))

# creating a vector with desired names of the files containing the mean of each gene set Normalized MFE
Curated_mean_dfs.output <- c("Mean_Normalized_MFE_MTOR_signaling", "Mean_Normalized_MFE_Cell_migration", "Mean_Normalized_MFE_Cell_motility", "Mean_Normalized_MFE_Gene_set_cell_adhesion", "Mean_Normalized_MFE_Chemokines_and_chemotaxis", "Mean_Normalized_MFE_Gene_set_Splicing_Factors", "Mean_Normalized_MFE_Homologous_recombination", "Mean_Normalized_MFE_Nucleotide_excision_repair", "Mean_Normalized_MFE_Non_hologous_end_joining", "Mean_Normalized_MFE_Fanconi_anemia", "Mean_Normalized_MFE_Base_excision_repair", "Mean_Normalized_MFE_Mismatch_repair", "Mean_Normalized_MFE_High_fidelity_polymerases", "Mean_Normalized_MFE_Error_prone_DNA_polymerases", "Mean_Normalized_MFE_ABC_SLC_transporters", "Mean_Normalized_MFE_calcium_channel_activity", "Mean_Normalized_MFE_TNF_signaling_pathway", "Mean_Normalized_MFE_Chaperone_mediated_protein_folding", "Mean_Normalized_MFE_TGFb_BMP_signaling_pathway", "Mean_Normalized_MFE_NF_kB_Target_genes", "Mean_Normalized_MFE_TGFb_signaling_targets", "Mean_Normalized_MFE_Notch_signaling_pathway", "Mean_Normalized_MFE_Hedgehog_Signaling_Pathway", "Mean_Normalized_MFE_WNT_signaling_pathway", "Mean_Normalized_MFE_WNT_signaling_targets", "Mean_Normalized_MFE_AMPK_signaling", "Mean_Normalized_MFE_cAMP_calcium_signaling_pathway", "Mean_Normalized_MFE_EGF_PDGF_signaling_pathway", "Mean_Normalized_MFE_Hypoxia_signaling_pathway", "Mean_Normalized_MFE_IL6_STAT3_signaling_pathway", "Mean_Normalized_MFE_Insulin_signaling_pathway", "Mean_Normalized_MFE_JAK_STAT_signaling_pathways", "Mean_Normalized_MFE_MAP_Kinase_signaling_pathways", "Mean_Normalized_MFE_p53_signaling_pathway", "Mean_Normalized_MFE_PI3K_AKT_signaling_pathway", "Mean_Normalized_MFE_TNF_signaling_pathway", "Mean_Normalized_MFE_Toll_Like_receptor_signaling_pathway", "Mean_Normalized_MFE_VEGF_signaling", "Mean_Normalized_MFE_Cellular_Stress_Responses", "Mean_Normalized_MFE_Antibacterial_response", "Mean_Normalized_MFE_Antifungal_response", "Mean_Normalized_MFE_Antiviral_response", "Mean_Normalized_MFE_Apoptosis", "Mean_Normalized_MFE_Autophagy", "Mean_Normalized_MFE_Cancer_drug_resistance", "Mean_Normalized_MFE_JNK_signaling_Vivanco", "Mean_Normalized_MFE_Pro_survival", "Mean_Normalized_MFE_Chromatin_remodelling")

# loop to assign data.frames to the desired variable name
for (i in (seq_along(Curated_mean_dfs.output))) {
  do.call("<-",list(Curated_mean_dfs.output[i], Curated_mean_test_df[[i]]))
}




# Generate the function for to calculate the STANDRAD DEVIATION of the normalized_MFE of each gene set resulting from the fuiltering here above.
STDEV_curated <- function(x) {
  stdev_curated<-sd(x$Normalized_MFE)
  return(stdev_curated)
}

# Apply the mean function to create a list of data.frames
Curated_stdev_test_df <- lapply(Curated_input_dfs, function(x) STDEV_curated(x))

# creating a vector with desired names of the files containing the mean of each gene set Normalized MFE
Curated_STDEV_dfs.output <- c("STDEV_Normalized_MFE_MTOR_signaling", "STDEV_Normalized_MFE_Cell_migration", "STDEV_Normalized_MFE_Cell_motility", "STDEV_Normalized_MFE_Gene_set_cell_adhesion", "STDEV_Normalized_MFE_Chemokines_and_chemotaxis", "STDEV_Normalized_MFE_Gene_set_Splicing_Factors", "STDEV_Normalized_MFE_Homologous_recombination", "STDEV_Normalized_MFE_Nucleotide_excision_repair", "STDEV_Normalized_MFE_Non_hologous_end_joining", "STDEV_Normalized_MFE_Fanconi_anemia", "STDEV_Normalized_MFE_Base_excision_repair", "STDEV_Normalized_MFE_Mismatch_repair", "STDEV_Normalized_MFE_High_fidelity_polymerases", "STDEV_Normalized_MFE_Error_prone_DNA_polymerases", "STDEV_Normalized_MFE_ABC_SLC_transporters", "STDEV_Normalized_MFE_calcium_channel_activity", "STDEV_Normalized_MFE_TNF_signaling_pathway", "STDEV_Normalized_MFE_Chaperone_mediated_protein_folding", "STDEV_Normalized_MFE_TGFb_BMP_signaling_pathway", "STDEV_Normalized_MFE_NF_kB_Target_genes", "STDEV_Normalized_MFE_TGFb_signaling_targets", "STDEV_Normalized_MFE_Notch_signaling_pathway", "STDEV_Normalized_MFE_Hedgehog_Signaling_Pathway", "STDEV_Normalized_MFE_WNT_signaling_pathway", "STDEV_Normalized_MFE_WNT_signaling_targets", "STDEV_Normalized_MFE_AMPK_signaling", "STDEV_Normalized_MFE_cAMP_calcium_signaling_pathway", "STDEV_Normalized_MFE_EGF_PDGF_signaling_pathway", "STDEV_Normalized_MFE_Hypoxia_signaling_pathway", "STDEV_Normalized_MFE_IL6_STAT3_signaling_pathway", "STDEV_Normalized_MFE_Insulin_signaling_pathway", "STDEV_Normalized_MFE_JAK_STAT_signaling_pathways", "STDEV_Normalized_MFE_MAP_Kinase_signaling_pathways", "STDEV_Normalized_MFE_p53_signaling_pathway", "STDEV_Normalized_MFE_PI3K_AKT_signaling_pathway", "STDEV_Normalized_MFE_TNF_signaling_pathway", "STDEV_Normalized_MFE_Toll_Like_receptor_signaling_pathway", "STDEV_Normalized_MFE_VEGF_signaling", "STDEV_Normalized_MFE_Cellular_Stress_Responses", "STDEV_Normalized_MFE_Antibacterial_response", "STDEV_Normalized_MFE_Antifungal_response", "STDEV_Normalized_MFE_Antiviral_response", "STDEV_Normalized_MFE_Apoptosis", "STDEV_Normalized_MFE_Autophagy", "STDEV_Normalized_MFE_Cancer_drug_resistance", "STDEV_Normalized_MFE_JNK_signaling_Vivanco", "STDEV_Normalized_MFE_Pro_survival", "STDEV_Normalized_MFE_Chromatin_remodelling")

# loop to assign data.frames to the desired variable name
for (i in (seq_along(Curated_STDEV_dfs.output))) {
  do.call("<-",list(Curated_STDEV_dfs.output[i], Curated_stdev_test_df[[i]]))
}



# Generate the data frame for the plot
df<-data.frame(mean=c(Mean_Normalized_MFE_MTOR_signaling, Mean_Normalized_MFE_Cell_migration, Mean_Normalized_MFE_Cell_motility, Mean_Normalized_MFE_Gene_set_cell_adhesion, Mean_Normalized_MFE_Chemokines_and_chemotaxis, Mean_Normalized_MFE_Gene_set_Splicing_Factors, Mean_Normalized_MFE_Homologous_recombination, Mean_Normalized_MFE_Nucleotide_excision_repair, Mean_Normalized_MFE_Non_hologous_end_joining, Mean_Normalized_MFE_Fanconi_anemia, Mean_Normalized_MFE_Base_excision_repair, Mean_Normalized_MFE_Mismatch_repair, Mean_Normalized_MFE_High_fidelity_polymerases, Mean_Normalized_MFE_Error_prone_DNA_polymerases, Mean_Normalized_MFE_ABC_SLC_transporters, Mean_Normalized_MFE_calcium_channel_activity, Mean_Normalized_MFE_TNF_signaling_pathway, Mean_Normalized_MFE_Chaperone_mediated_protein_folding, Mean_Normalized_MFE_TGFb_BMP_signaling_pathway, Mean_Normalized_MFE_NF_kB_Target_genes, Mean_Normalized_MFE_TGFb_signaling_targets, Mean_Normalized_MFE_Notch_signaling_pathway, Mean_Normalized_MFE_Hedgehog_Signaling_Pathway, Mean_Normalized_MFE_WNT_signaling_pathway, Mean_Normalized_MFE_WNT_signaling_targets, Mean_Normalized_MFE_AMPK_signaling, Mean_Normalized_MFE_cAMP_calcium_signaling_pathway, Mean_Normalized_MFE_EGF_PDGF_signaling_pathway, Mean_Normalized_MFE_Hypoxia_signaling_pathway, Mean_Normalized_MFE_IL6_STAT3_signaling_pathway, Mean_Normalized_MFE_Insulin_signaling_pathway, Mean_Normalized_MFE_JAK_STAT_signaling_pathways, Mean_Normalized_MFE_MAP_Kinase_signaling_pathways, Mean_Normalized_MFE_p53_signaling_pathway, Mean_Normalized_MFE_PI3K_AKT_signaling_pathway, Mean_Normalized_MFE_TNF_signaling_pathway, Mean_Normalized_MFE_Toll_Like_receptor_signaling_pathway, Mean_Normalized_MFE_VEGF_signaling, Mean_Normalized_MFE_Cellular_Stress_Responses, Mean_Normalized_MFE_Antibacterial_response, Mean_Normalized_MFE_Antifungal_response, Mean_Normalized_MFE_Antiviral_response, Mean_Normalized_MFE_Apoptosis, Mean_Normalized_MFE_Autophagy, Mean_Normalized_MFE_Cancer_drug_resistance, Mean_Normalized_MFE_JNK_signaling_Vivanco, Mean_Normalized_MFE_Pro_survival, Mean_Normalized_MFE_Chromatin_remodelling),
               sd=c(STDEV_Normalized_MFE_MTOR_signaling, STDEV_Normalized_MFE_Cell_migration, STDEV_Normalized_MFE_Cell_motility, STDEV_Normalized_MFE_Gene_set_cell_adhesion, STDEV_Normalized_MFE_Chemokines_and_chemotaxis, STDEV_Normalized_MFE_Gene_set_Splicing_Factors, STDEV_Normalized_MFE_Homologous_recombination, STDEV_Normalized_MFE_Nucleotide_excision_repair, STDEV_Normalized_MFE_Non_hologous_end_joining, STDEV_Normalized_MFE_Fanconi_anemia, STDEV_Normalized_MFE_Base_excision_repair, STDEV_Normalized_MFE_Mismatch_repair, STDEV_Normalized_MFE_High_fidelity_polymerases, STDEV_Normalized_MFE_Error_prone_DNA_polymerases, STDEV_Normalized_MFE_ABC_SLC_transporters, STDEV_Normalized_MFE_calcium_channel_activity, STDEV_Normalized_MFE_TNF_signaling_pathway, STDEV_Normalized_MFE_Chaperone_mediated_protein_folding, STDEV_Normalized_MFE_TGFb_BMP_signaling_pathway, STDEV_Normalized_MFE_NF_kB_Target_genes, STDEV_Normalized_MFE_TGFb_signaling_targets, STDEV_Normalized_MFE_Notch_signaling_pathway, STDEV_Normalized_MFE_Hedgehog_Signaling_Pathway, STDEV_Normalized_MFE_WNT_signaling_pathway, STDEV_Normalized_MFE_WNT_signaling_targets, STDEV_Normalized_MFE_AMPK_signaling, STDEV_Normalized_MFE_cAMP_calcium_signaling_pathway, STDEV_Normalized_MFE_EGF_PDGF_signaling_pathway, STDEV_Normalized_MFE_Hypoxia_signaling_pathway, STDEV_Normalized_MFE_IL6_STAT3_signaling_pathway, STDEV_Normalized_MFE_Insulin_signaling_pathway, STDEV_Normalized_MFE_JAK_STAT_signaling_pathways, STDEV_Normalized_MFE_MAP_Kinase_signaling_pathways, STDEV_Normalized_MFE_p53_signaling_pathway, STDEV_Normalized_MFE_PI3K_AKT_signaling_pathway, STDEV_Normalized_MFE_TNF_signaling_pathway, STDEV_Normalized_MFE_Toll_Like_receptor_signaling_pathway, STDEV_Normalized_MFE_VEGF_signaling, STDEV_Normalized_MFE_Cellular_Stress_Responses, STDEV_Normalized_MFE_Antibacterial_response, STDEV_Normalized_MFE_Antifungal_response, STDEV_Normalized_MFE_Antiviral_response, STDEV_Normalized_MFE_Apoptosis, STDEV_Normalized_MFE_Autophagy, STDEV_Normalized_MFE_Cancer_drug_resistance, STDEV_Normalized_MFE_JNK_signaling_Vivanco, STDEV_Normalized_MFE_Pro_survival, STDEV_Normalized_MFE_Chromatin_remodelling),
               Category=as.factor(c("MTOR_signaling", "Cell_migration", "Cell_motility", "Gene_set_cell_adhesion", "Chemokines_and_chemotaxis", "Gene_set_Splicing_Factors", "Homologous_recombination", "Nucleotide_excision_repair", "Non_hologous_end_joining", "Fanconi_anemia", "Base_excision_repair", "Mismatch_repair", "High_fidelity_polymerases", "Error_prone_DNA_polymerases", "ABC_SLC_transporters", "calcium_channel_activity", "TNF_signaling_pathway", "Chaperone_mediated_protein_folding", "TGFb_BMP_signaling_pathway", "NF_kB_Target_genes", "TGFb_signaling_targets", "Notch_signaling_pathway", "Hedgehog_Signaling_Pathway", "WNT_signaling_pathway", "WNT_signaling_targets", "AMPK_signaling", "cAMP_calcium_signaling_pathway", "EGF_PDGF_signaling_pathway", "Hypoxia_signaling_pathway", "IL6_STAT3_signaling_pathway", "Insulin_signaling_pathway", "JAK_STAT_signaling_pathways", "MAP_Kinase_signaling_pathways", "p53_signaling_pathway", "PI3K_AKT_signaling_pathway", "TNF_signaling_pathway", "Toll_Like_receptor_signaling_pathway", "VEGF_signaling", "Cellular_Stress_Responses", "Antibacterial_response", "Antifungal_response", "Antiviral_response", "Apoptosis", "Autophagy", "Cancer_drug_resistance", "JNK_signaling_Vivanco", "Pro_survival", "Chromatin_remodelling")))
df

# Calculate the mean of all genes in each category
Mean_all_pathways <- df %>%
  dplyr::summarise(mean_all_pathways = mean(mean))

library(ggplot2)
ggplot(df, aes(x=reorder(Category, -mean), y=mean)) +
  theme_bw() +
  geom_point(position=position_dodge(width = 1), stat="identity",colour='black') +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.4) +
  theme(axis.text.x = element_text(size=12, angle=0,hjust=0.95,vjust=0.98)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 10))) +
  geom_point(colour = "black", size = 2) +
  # ylim(-0.5, -0.2) +
  geom_hline(yintercept = Mean_all_pathways$mean_all_pathways, color = "black") +
  xlab("") + ylab("Normalized_MFE (kcal/mol)") +
  ggtitle("Normalized_MFE of curated signaling pathways (all transcript isoforms)") +
  theme(plot.title = element_text(color="black", size=12, face="bold")) +
  coord_flip()








# Create the box plot
# IMPORTANT: 
# Some of the gene sets have been removed because they were not relevant for the analysis of the MFE of signaling pathways and DNA repair processes.

# The gene sets that have been removed are:
# df_Cell_migration,
# df_Cell_motility,
# df_Gene_set_cell_adhesion,
# df_Chemokines_and_chemotaxis,
# df_Gene_set_Splicing_Factors,
# df_ABC_SLC_transporters,
# df_calcium_channel_activity,
# df_Chaperone_mediated_protein_folding,
# df_NF_kB_Target_genes,
# df_TGFb_signaling_targets,
# df_WNT_signaling_targets,
# df_Cellular_Stress_Responses,
# df_Antibacterial_response,
# df_Antifungal_response,
# df_Antiviral_response,
# df_Apoptosis,
# df_Autophagy,
# df_Cancer_drug_resistance,
# df_Pro_survival,
# df_Chromatin_remodelling

# Function to select specific columns from a list of data frames
select_columns_from_dfs <- function(df_list, columns) {
  # Use lapply to apply the select function to each data frame in the list
  selected_dfs <- lapply(df_list, function(df) {
    dplyr::select(df, all_of(columns))
  })
  # Return the list of modified data frames
  return(selected_dfs)
}

# List of data frames
df_list <- list(
  df_MTOR_signaling,
  df_Homologous_recombination,
  df_Nucleotide_excision_repair,
  df_Non_hologous_end_joining,
  df_Fanconi_anemia,
  df_Base_excision_repair,
  df_Mismatch_repair,
  df_High_fidelity_polymerases,
  df_Error_prone_DNA_polymerases,
  df_TNF_signaling_pathway,
  df_TGFb_BMP_signaling_pathway,
  df_Notch_signaling_pathway,
  df_Hedgehog_Signaling_Pathway,
  df_WNT_signaling_pathway,
  df_AMPK_signaling,
  df_cAMP_calcium_signaling_pathway,
  df_EGF_PDGF_signaling_pathway,
  df_Hypoxia_signaling_pathway,
  df_IL6_STAT3_signaling_pathway,
  df_Insulin_signaling_pathway,
  df_JAK_STAT_signaling_pathways,
  df_MAP_Kinase_signaling_pathways,
  df_p53_signaling_pathway,
  df_PI3K_AKT_signaling_pathway,
  df_Toll_Like_receptor_signaling_pathway,
  df_VEGF_signaling,
  df_JNK_signaling_Vivanco
)

# Remove these samples from the list: df_Genes_UPregulated_in_dormancy_Gelman_2023




# Columns to select
columns <- c("external_gene_name", "ensembl_transcript_id_version", "Normalized_MFE")

# Apply the function to select columns from all data frames
selected_dfs <- select_columns_from_dfs(df_list, columns)

# Naming the selected data frames for easy identification
names(selected_dfs) <- c(
  "MTOR_signaling",
  "Homologous_recombination",
  "Nucleotide_excision_repair",
  "Non_hologous_end_joining",
  "Fanconi_anemia",
  "Base_excision_repair",
  "Mismatch_repair",
  "High_fidelity_polymerases",
  "Error_prone_DNA_polymerases",
  "TNF_signaling",
  "TGFb_BMP_signaling",
  "Notch_signaling",
  "Hedgehog_Signaling",
  "WNT_signaling",
  "AMPK_signaling",
  "cAMP_calcium_signaling",
  "EGF_PDGF_signaling",
  "Hypoxia_signaling",
  "IL6_STAT3_signaling",
  "Insulin_signaling",
  "JAK_STAT_signaling",
  "MAP_Kinase_signaling",
  "p53_signaling",
  "PI3K_AKT_signaling",
  "Toll_Like_receptor_signaling",
  "VEGF_signaling",
  "JNK_signaling_Vivanco"
)

# Print the first few rows of each selected data frame for verification
lapply(names(selected_dfs), function(name) {
  cat("Data frame:", name, "\n")
  print(head(selected_dfs[[name]]))
  cat("\n")
})


# Combine the selected data frames into a single data frame with a source column
combine_dfs <- function(df_list, names_list) {
  combined_df <- bind_rows(lapply(1:length(df_list), function(i) {
    df_list[[i]] %>% mutate(source = names_list[i])
  }))
  return(combined_df)
}

# Combine all the selected data frames
combined_df <- combine_dfs(selected_dfs, names(selected_dfs))

# Calculate the mean of Normalized_MFE for each source
mean_df <- combined_df %>%
  group_by(source) %>%
  summarize(Mean_Norm_MFE = mean(Normalized_MFE, na.rm = TRUE))

# Merge the mean values back into the combined data frame
combined_df <- left_join(combined_df, mean_df, by = "source")


# Create a box plot using ggplot2
MFE_curated_signaling_pathways <- ggplot(combined_df, aes(x = reorder(source, -Normalized_MFE), y = Normalized_MFE)) +
  geom_boxplot() +
  theme_classic() +
  geom_point(data = mean_df, aes(x = source, y = Mean_Norm_MFE), color = "black", size = 2, shape = 18) +  # Add a dot for the mean
  geom_hline(yintercept = Mean_all_pathways$mean_all_pathways, color = "black") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)) +
  ylim(-0.6, 0) +
  labs(title = "",
       x = "",
       y = "") +

coord_flip()

print(MFE_curated_signaling_pathways)





# Inspect the mean MFE values of the homologous recombination and the base excision repair gene set
# Add a column indicating the group (pathway name)
df_HR <- df_Homologous_recombination %>%
  dplyr::select(external_gene_name, Normalized_MFE) %>%
  mutate(Pathway = "Homologous_recombination")

df_BER <- df_Base_excision_repair %>%
  dplyr::select(external_gene_name, Normalized_MFE) %>%
  mutate(Pathway = "Base_excision_repair")

# Combine the two into a single data frame
df_HR_BER <- bind_rows(df_HR, df_BER)


# Check whether they are statistically different using the Wilcoxon test
wilcox_test <- wilcox.test(Normalized_MFE ~ Pathway, data = df_HR_BER)
print(wilcox_test)


# Calculate the mean and standard deviation for each pathway
mean_HR  <- mean(df_Homologous_recombination$Normalized_MFE, na.rm = TRUE)
sd_HR    <- sd(df_Homologous_recombination$Normalized_MFE, na.rm = TRUE)
HR_str   <- sprintf("%.3f \u00B1 %.3f kcal/mol", mean_HR, sd_HR)

mean_BER <- mean(df_Base_excision_repair$Normalized_MFE, na.rm = TRUE)
sd_BER   <- sd(df_Base_excision_repair$Normalized_MFE, na.rm = TRUE)
BER_str  <- sprintf("%.3f \u00B1 %.3f kcal/mol", mean_BER, sd_BER)

HR_str
BER_str

























# Supplementary Fig.2 ####
# GSEA gene sets
# Here below, I repeat the filtering automatic function using the GSEA gene sets

# The data here above have been obtained using curated pathways generated by myself and other authors
# It's better to confirm the data using KEGG signaling pathways
# Import the GSEA pathways from; http://www.gsea-msigdb.org/gsea/downloads.jsp

library(rjson)

# Import the specific signaling pathways from the GSEA website
# It is also possible to import the GSEA canonical pathway gene set and select the KEGG pathways for the analysis
# Example:
GSEA_canonical_pathways_gene_sets <- fromJSON(file="path_to/c2.cp.v7.5.1.json")

# To convert the json file to data frame, we need to unlist the json_file
GSEA_canonical_pathways_gene_sets <- lapply(GSEA_canonical_pathways_gene_sets, function(x) {
  x[sapply(x, is.null)] <- NA
  unlist(x)
})



GSEA_canonical_pathways_gene_sets_df<-as.data.frame(do.call("cbind", GSEA_canonical_pathways_gene_sets))

# To obtain the different signaling pathways from the data frame, I convert the specific columns of the data frame to a string vector with pull Function in the dplyr Package
# Example for the MTOR signaling pathway
GSEA_KEGG_MTOR_SIGNALING_PATHWAY<- dplyr::pull(GSEA_canonical_pathways_gene_sets_df,KEGG_MTOR_SIGNALING_PATHWAY)
# Remove elements of the string vector with empty character ""
library(stringi)
GSEA_KEGG_MTOR_SIGNALING_PATHWAY <- as.character(stri_remove_empty(GSEA_KEGG_MTOR_SIGNALING_PATHWAY, na_empty = FALSE))




# Here below I have downloaded the gene sets of the specific pathways analyzed here above, by accessing the GSEA database 
# The json files have been saved in a folder of Exp 277

GSEA_cell_motility <- fromJSON(file="path_to/GOBP_ACTIN_POLYMERIZATION_DEPENDENT_CELL_MOTILITY.v7.5.1.json")
GSEA_cell_adhesion <- fromJSON(file="path_to/CELL_CELL_ADHESION.v7.5.1.json")
GSEA_chemotaxis <- fromJSON(file="path_to/GOBP_CELL_CHEMOTAXIS.v7.5.1.json")
GSEA_splicing <- fromJSON(file="path_to/GOBP_ALTERNATIVE_MRNA_SPLICING_VIA_SPLICEOSOME.v7.5.1.json")
GSEA_homologous_recombination <- fromJSON(file="path_to/GOBP_HOMOLOGOUS_RECOMBINATION.v7.5.1.json")
GSEA_NER <- fromJSON(file="path_to/GOBP_NUCLEOTIDE_EXCISION_REPAIR.v7.5.1.json")
GSEA_NHEJ <- fromJSON(file="path_to/KEGG_NON_HOMOLOGOUS_END_JOINING.v7.5.1.json")
GSEA_FA <- fromJSON(file="path_to/REACTOME_FANCONI_ANEMIA_PATHWAY.v7.5.1.json")
GSEA_BER <- fromJSON(file="path_to/KEGG_BASE_EXCISION_REPAIR.v7.5.1.json")
GSEA_MMR <- fromJSON(file="path_to/KEGG_MISMATCH_REPAIR.v7.5.1.json")
GSEA_ABC_transporters <- fromJSON(file="path_to/KEGG_ABC_TRANSPORTERS.v7.5.1.json")
GSEA_Ca_channels_activity <- fromJSON(file="path_to/GOBP_CALCIUM_ION_TRANSMEMBRANE_TRANSPORT_VIA_HIGH_VOLTAGE_GATED_CALCIUM_CHANNEL.v7.5.1.json")
GSEA_TNF_signaling <- fromJSON(file="path_to/REACTOME_TNF_SIGNALING.v7.5.1.json")
GSEA_Chaperone_protein_folding <- fromJSON(file="path_to/GOMF_PROTEIN_FOLDING_CHAPERONE.v7.5.1.json")
GSEA_TGF_BETA <- fromJSON(file="path_to/KEGG_TGF_BETA_SIGNALING_PATHWAY.v7.5.1.json")
GSEA_BMP <- fromJSON(file="path_to/REACTOME_SIGNALING_BY_BMP.v7.5.1.json")
GSEA_NFKB <- fromJSON(file="path_to/WP_CANONICAL_NFKB_PATHWAY.v7.5.1.json")
GSEA_NOTCH <- fromJSON(file="path_to/KEGG_NOTCH_SIGNALING_PATHWAY.v7.5.1.json")
GSEA_HEDGEHOG <- fromJSON(file="path_to/KEGG_HEDGEHOG_SIGNALING_PATHWAY.v7.5.1.json")
GSEA_WNT <- fromJSON(file="path_to/KEGG_WNT_SIGNALING_PATHWAY.v7.5.1.json")
GSEA_AMPK <- fromJSON(file="path_to/REACTOME_REGULATION_OF_AMPK_ACTIVITY_VIA_LKB1.v7.5.1.json")
GSEA_EGFR <- fromJSON(file="path_to/REACTOME_SIGNALING_BY_EGFR.v7.5.1.json")
GSEA_HIPPO <- fromJSON(file="path_to/REACTOME_SIGNALING_BY_HIPPO.v7.5.1.json")
GSEA_HYPOXIA <- fromJSON(file="path_to/REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA.v7.5.1.json")
GSEA_IL6_JAK_STAT3 <- fromJSON(file="path_to/HALLMARK_IL6_JAK_STAT3_SIGNALING.v7.5.1.json")
GSEA_Insulin <- fromJSON(file="path_to/KEGG_INSULIN_SIGNALING_PATHWAY.v7.5.1.json")
GSEA_JAK_STAT <- fromJSON(file="path_to/KEGG_JAK_STAT_SIGNALING_PATHWAY.v7.5.1.json")
GSEA_MAPK <- fromJSON(file="path_to/REACTOME_RAF_MAP_KINASE_CASCADE.v7.5.1.json")
GSEA_P53 <- fromJSON(file="path_to/KEGG_P53_SIGNALING_PATHWAY.v7.5.1.json")
GSEA_PI3K_AKT <- fromJSON(file="path_to/REACTOME_PI3K_AKT_SIGNALING_IN_CANCER.v7.5.1.json")
GSEA_TOLL_LIKE_RECEPTOR <- fromJSON(file="path_to/REACTOME_TOLL_LIKE_RECEPTOR_TLR1_TLR2_CASCADE.v7.5.1.json")
GSEA_VEGF <- fromJSON(file="path_to/KEGG_VEGF_SIGNALING_PATHWAY.v7.5.1.json")
GSEA_Stress_response <- fromJSON(file="path_to/REACTOME_CELLULAR_RESPONSE_TO_CHEMICAL_STRESS.v7.5.1.json")
GSEA_Bacterial_response <- fromJSON(file="path_to/GOBP_CELLULAR_RESPONSE_TO_BACTERIAL_LIPOPROTEIN.v7.5.1.json")
GSEA_Antifungal_response <- fromJSON(file="path_to/GOBP_RESPONSE_TO_MOLECULE_OF_FUNGAL_ORIGIN.v7.5.1.json")
GSEA_Antiviral_response <- fromJSON(file="path_to/REACTOME_OAS_ANTIVIRAL_RESPONSE.v7.5.1.json")
GSEA_Apoptosis <- fromJSON(file="path_to/KEGG_APOPTOSIS.v7.5.1.json")
GSEA_Autophagy <- fromJSON(file="path_to/KEGG_REGULATION_OF_AUTOPHAGY.v7.5.1.json")
GSEA_Drug_resistance <- fromJSON(file="path_to/KESHELAVA_MULTIPLE_DRUG_RESISTANCE.v7.5.1.json")
GSEA_JNK <- fromJSON(file="path_to/GOBP_POSITIVE_REGULATION_OF_JNK_CASCADE.v7.5.1.json")
GSEA_Chromatin_remodelling <- fromJSON(file="path_to/GOBP_CHROMATIN_REMODELING.v7.5.1.json")
GSEA_CELL_CYCLE_G1_S_PHASE_TRANSITION <- fromJSON(file="path_to/GOBP_CELL_CYCLE_G1_S_PHASE_TRANSITION.v7.5.1.json")
GSEA_MITOTIC_SPINDLE_ORGANIZATION <- fromJSON(file="path_to/GOBP_MITOTIC_SPINDLE_ORGANIZATION.v7.5.1.json")
GSEA_POSITIVE_REGULATION_OF_CELL_MIGRATION <- fromJSON(file="path_to/GOBP_POSITIVE_REGULATION_OF_CELL_MIGRATION_BY_VASCULAR_ENDOTHELIAL_GROWTH_FACTOR_SIGNALING_PATHWAY.v7.5.1.json")
GSEA_CELL_MIGRATION <- fromJSON(file="path_to/GOBP_CELL_MIGRATION.v7.5.1.json")
GSEA_POSITIVE_REGULATION_OF_ERK1_AND_ERK2_CASCADE <- fromJSON(file="path_to/GOBP_POSITIVE_REGULATION_OF_ERK1_AND_ERK2_CASCADE.v7.5.1.json")
GSEA_RESPONSE_TO_TYPE_I_INTERFERON <- fromJSON(file="path_to/GOBP_RESPONSE_TO_TYPE_I_INTERFERON.v7.5.1.json")
GSEA_RESPONSE_TO_INTERFERON_GAMMA <- fromJSON(file="path_to/GOBP_CELLULAR_RESPONSE_TO_INTERFERON_GAMMA.v7.5.1.json")
GSEA_GOBP_EMT <- fromJSON(file="path_to/GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION.v7.5.1.json")



# For each GSEA gene set here above, these commands generate a list of the gene symbols and use these lists to filter the file 'Table_ANNOTATED_gene_names'
# Generate the function that filters each gene set in the target data frame
filtering <- function(x) {
  Filtering<-dplyr::filter(Table_ANNOTATED_gene_names, external_gene_name %in% x)
  return(Filtering)
}

# Generate the list containing all the gene sets that have to be analyzed.
GSEA_PATHWAYS.input <- list(GSEA_cell_motility, GSEA_cell_adhesion, GSEA_chemotaxis, GSEA_splicing, GSEA_homologous_recombination, GSEA_NER, GSEA_NHEJ, GSEA_FA, GSEA_BER, GSEA_MMR, GSEA_ABC_transporters, GSEA_Ca_channels_activity, GSEA_TNF_signaling, GSEA_Chaperone_protein_folding, GSEA_TGF_BETA, GSEA_BMP, GSEA_NFKB, GSEA_NOTCH, GSEA_HEDGEHOG, GSEA_WNT, GSEA_AMPK, GSEA_EGFR, GSEA_HIPPO, GSEA_HYPOXIA, GSEA_IL6_JAK_STAT3, GSEA_Insulin, GSEA_JAK_STAT, GSEA_MAPK, GSEA_P53, GSEA_PI3K_AKT, GSEA_TOLL_LIKE_RECEPTOR, GSEA_VEGF, GSEA_Stress_response, GSEA_Bacterial_response, GSEA_Antifungal_response, GSEA_Antiviral_response, GSEA_Apoptosis, GSEA_Autophagy, GSEA_Drug_resistance, GSEA_JNK, GSEA_Chromatin_remodelling, GSEA_CELL_CYCLE_G1_S_PHASE_TRANSITION, GSEA_MITOTIC_SPINDLE_ORGANIZATION, GSEA_POSITIVE_REGULATION_OF_CELL_MIGRATION, GSEA_CELL_MIGRATION, GSEA_POSITIVE_REGULATION_OF_ERK1_AND_ERK2_CASCADE, GSEA_RESPONSE_TO_TYPE_I_INTERFERON, GSEA_RESPONSE_TO_INTERFERON_GAMMA, GSEA_GOBP_EMT)

# This creates the names that will be used to for the final data frames of each file resulting from the filtering of the target data frame with each gene set
GSEA_PATHWAYS.output <- c("GSEA_cell_motility", "GSEA_cell_adhesion", "GSEA_chemotaxis", "GSEA_splicing", "GSEA_homologous_recombination", "GSEA_NER", "GSEA_NHEJ", "GSEA_FA", "GSEA_BER", "GSEA_MMR", "GSEA_ABC_transporters", "GSEA_Ca_channels_activity", "GSEA_TNF_signaling", "GSEA_Chaperone_protein_folding", "GSEA_TGF_BETA", "GSEA_BMP", "GSEA_NFKB", "GSEA_NOTCH", "GSEA_HEDGEHOG", "GSEA_WNT", "GSEA_AMPK", "GSEA_EGFR", "GSEA_HIPPO", "GSEA_HYPOXIA", "GSEA_IL6_JAK_STAT3", "GSEA_Insulin", "GSEA_JAK_STAT", "GSEA_MAPK", "GSEA_P53", "GSEA_PI3K_AKT", "GSEA_TOLL_LIKE_RECEPTOR", "GSEA_VEGF", "GSEA_Stress_response", "GSEA_Bacterial_response", "GSEA_Antifungal_response", "GSEA_Antiviral_response", "GSEA_Apoptosis", "GSEA_Autophagy", "GSEA_Drug_resistance", "GSEA_JNK", "GSEA_Chromatin_remodelling", "GSEA_CELL_CYCLE_G1_S_PHASE_TRANSITION", "GSEA_MITOTIC_SPINDLE_ORGANIZATION", "GSEA_POSITIVE_REGULATION_OF_CELL_MIGRATION", "GSEA_CELL_MIGRATION", "GSEA_POSITIVE_REGULATION_OF_ERK1_AND_ERK2_CASCADE", "GSEA_RESPONSE_TO_TYPE_I_INTERFERON", "GSEA_RESPONSE_TO_INTERFERON_GAMMA", "GSEA_GOBP_EMT")

# extract gene symbols from the json files
library("rrapply")
test <- lapply(seq_along(GSEA_PATHWAYS.output), function(x) rrapply(GSEA_PATHWAYS.input[[x]], how = "flatten")$geneSymbols)

# assign names to each entry in the list
names(test) <- GSEA_PATHWAYS.output

# applying filtering function to do the filtering to create a list of data.frames
test_df <- lapply(test, function(x) filtering(x))

# creating a vector with desired names of filtered data.frames
dfs.output <- c("df_GSEA_cell_motility", "df_GSEA_cell_adhesion", "df_GSEA_chemotaxis", "df_GSEA_splicing", "df_GSEA_homologous_recombination", "df_GSEA_NER", "df_GSEA_NHEJ", "df_GSEA_FA", "df_GSEA_BER", "df_GSEA_MMR", "df_GSEA_ABC_transporters", "df_GSEA_Ca_channels_activity", "df_GSEA_TNF_signaling", "df_GSEA_Chaperone_protein_folding", "df_GSEA_TGF_BETA", "df_GSEA_BMP", "df_GSEA_NFKB", "df_GSEA_NOTCH", "df_GSEA_HEDGEHOG", "df_GSEA_WNT", "df_GSEA_AMPK", "df_GSEA_EGFR", "df_GSEA_HIPPO", "df_GSEA_HYPOXIA", "df_GSEA_IL6_JAK_STAT3", "df_GSEA_Insulin", "df_GSEA_JAK_STAT", "df_GSEA_MAPK", "df_GSEA_P53", "df_GSEA_PI3K_AKT", "df_GSEA_TOLL_LIKE_RECEPTOR", "df_GSEA_VEGF", "df_GSEA_Stress_response", "df_GSEA_Bacterial_response", "df_GSEA_Antifungal_response", "df_GSEA_Antiviral_response", "df_GSEA_Apoptosis", "df_GSEA_Autophagy", "df_GSEA_Drug_resistance", "df_GSEA_JNK", "df_GSEA_Chromatin_remodelling", "df_GSEA_CELL_CYCLE_G1_S_PHASE_TRANSITION", "df_GSEA_MITOTIC_SPINDLE_ORGANIZATION", "df_GSEA_POSITIVE_REGULATION_OF_CELL_MIGRATION", "df_GSEA_CELL_MIGRATION", "df_GSEA_POSITIVE_REGULATION_OF_ERK1_AND_ERK2_CASCADE", "df_GSEA_RESPONSE_TO_TYPE_I_INTERFERON", "df_GSEA_RESPONSE_TO_INTERFERON_GAMMA", "df_GSEA_GOBP_EMT")

# loop to assign data.frames to the desired variable name
for (i in (seq_along(dfs.output))) {
  do.call("<-",list(dfs.output[i], test_df[[i]]))
}






# Generate new functions 'mean' and 'sd' and loop them through the 'Normalized_MFE' columns of each data frame generated by the automatc filtering here above

# Generate the function for to calculate the MEAN of the normalized_MFE of each gene set resulting from the fuiltering here above.
Mean_GSEA <- function(x) {
  mean_GSEA<-mean(x$Normalized_MFE)
  return(mean_GSEA)
}

# Create an input of the data frames generated here above that I will use to calculate the mean of a specific column
GSEA_input_dfs <- list(df_GSEA_cell_motility, df_GSEA_cell_adhesion, df_GSEA_chemotaxis, df_GSEA_splicing, df_GSEA_homologous_recombination, df_GSEA_NER, df_GSEA_NHEJ, df_GSEA_FA, df_GSEA_BER, df_GSEA_MMR, df_GSEA_ABC_transporters, df_GSEA_Ca_channels_activity, df_GSEA_TNF_signaling, df_GSEA_Chaperone_protein_folding, df_GSEA_TGF_BETA, df_GSEA_BMP, df_GSEA_NFKB, df_GSEA_NOTCH, df_GSEA_HEDGEHOG, df_GSEA_WNT, df_GSEA_AMPK, df_GSEA_EGFR, df_GSEA_HIPPO, df_GSEA_HYPOXIA, df_GSEA_IL6_JAK_STAT3, df_GSEA_Insulin, df_GSEA_JAK_STAT, df_GSEA_MAPK, df_GSEA_P53, df_GSEA_PI3K_AKT, df_GSEA_TOLL_LIKE_RECEPTOR, df_GSEA_VEGF, df_GSEA_Stress_response, df_GSEA_Bacterial_response, df_GSEA_Antifungal_response, df_GSEA_Antiviral_response, df_GSEA_Apoptosis, df_GSEA_Autophagy, df_GSEA_Drug_resistance, df_GSEA_JNK, df_GSEA_Chromatin_remodelling, df_GSEA_CELL_CYCLE_G1_S_PHASE_TRANSITION, df_GSEA_MITOTIC_SPINDLE_ORGANIZATION, df_GSEA_POSITIVE_REGULATION_OF_CELL_MIGRATION, df_GSEA_CELL_MIGRATION, df_GSEA_POSITIVE_REGULATION_OF_ERK1_AND_ERK2_CASCADE, df_GSEA_RESPONSE_TO_TYPE_I_INTERFERON, df_GSEA_RESPONSE_TO_INTERFERON_GAMMA, df_GSEA_GOBP_EMT)

# Apply the mean function to create a list of data.frames
GSEA_mean_test_df <- lapply(GSEA_input_dfs, function(x) Mean_GSEA(x))

# creating a vector with desired names of the files containing the mean of each gene set Normalized MFE
GSEA_mean_Normalized_MFEs.output <- c("Mean_Normalized_MFE_GSEA_cell_motility", "Mean_Normalized_MFE_GSEA_cell_adhesion", "Mean_Normalized_MFE_GSEA_chemotaxis", "Mean_Normalized_MFE_GSEA_splicing", "Mean_Normalized_MFE_GSEA_homologous_recombination", "Mean_Normalized_MFE_GSEA_NER", "Mean_Normalized_MFE_GSEA_NHEJ", "Mean_Normalized_MFE_GSEA_FA", "Mean_Normalized_MFE_GSEA_BER", "Mean_Normalized_MFE_GSEA_MMR", "Mean_Normalized_MFE_GSEA_ABC_transporters", "Mean_Normalized_MFE_GSEA_Ca_channels_activity", "Mean_Normalized_MFE_GSEA_TNF_signaling", "Mean_Normalized_MFE_GSEA_Chaperone_protein_folding", "Mean_Normalized_MFE_GSEA_TGF_BETA", "Mean_Normalized_MFE_GSEA_BMP", "Mean_Normalized_MFE_GSEA_NFKB", "Mean_Normalized_MFE_GSEA_NOTCH", "Mean_Normalized_MFE_GSEA_HEDGEHOG", "Mean_Normalized_MFE_GSEA_WNT", "Mean_Normalized_MFE_GSEA_AMPK", "Mean_Normalized_MFE_GSEA_EGFR", "Mean_Normalized_MFE_GSEA_HIPPO", "Mean_Normalized_MFE_GSEA_HYPOXIA", "Mean_Normalized_MFE_GSEA_IL6_JAK_STAT3", "Mean_Normalized_MFE_GSEA_Insulin", "Mean_Normalized_MFE_GSEA_JAK_STAT", "Mean_Normalized_MFE_GSEA_MAPK", "Mean_Normalized_MFE_GSEA_P53", "Mean_Normalized_MFE_GSEA_PI3K_AKT", "Mean_Normalized_MFE_GSEA_TOLL_LIKE_RECEPTOR", "Mean_Normalized_MFE_GSEA_VEGF", "Mean_Normalized_MFE_GSEA_Stress_response", "Mean_Normalized_MFE_GSEA_Bacterial_response", "Mean_Normalized_MFE_GSEA_Antifungal_response", "Mean_Normalized_MFE_GSEA_Antiviral_response", "Mean_Normalized_MFE_GSEA_Apoptosis", "Mean_Normalized_MFE_GSEA_Autophagy", "Mean_Normalized_MFE_GSEA_Drug_resistance", "Mean_Normalized_MFE_GSEA_JNK", "Mean_Normalized_MFE_GSEA_Chromatin_remodelling", "Mean_Normalized_MFE_GSEA_CELL_CYCLE_G1_S_PHASE_TRANSITION", "Mean_Normalized_MFE_GSEA_MITOTIC_SPINDLE_ORGANIZATION", "Mean_Normalized_MFE_GSEA_POSITIVE_REGULATION_OF_CELL_MIGRATION", "Mean_Normalized_MFE_GSEA_CELL_MIGRATION", "Mean_Normalized_MFE_GSEA_POSITIVE_REGULATION_OF_ERK1_AND_ERK2_CASCADE", "Mean_Normalized_MFE_GSEA_RESPONSE_TO_TYPE_I_INTERFERON", "Mean_Normalized_MFE_GSEA_RESPONSE_TO_INTERFERON_GAMMA", "Mean_Normalized_MFE_GSEA_GOBP_EMT")

# loop to assign data.frames to the desired variable name
for (i in (seq_along(GSEA_mean_Normalized_MFEs.output))) {
  do.call("<-",list(GSEA_mean_Normalized_MFEs.output[i], GSEA_mean_test_df[[i]]))
}







# Generate the function for to calculate the MEAN of the normalized_MFE of each gene set resulting from the fuiltering here above.
STDEV_GSEA <- function(x) {
  stdev_GSEA<-sd(x$Normalized_MFE)
  return(stdev_GSEA)
}

# Apply the STDEV function to create a list of data.frames
GSEA_STDEV_test_df <- lapply(GSEA_input_dfs, function(x) STDEV_GSEA(x))

# creating a vector with desired names of the files containing the mean of each gene set Normalized MFE
GSEA_STDEV_Normalized_MFEs.output <- c("STDEV_Normalized_MFE_GSEA_cell_motility", "STDEV_Normalized_MFE_GSEA_cell_adhesion", "STDEV_Normalized_MFE_GSEA_chemotaxis", "STDEV_Normalized_MFE_GSEA_splicing", "STDEV_Normalized_MFE_GSEA_homologous_recombination", "STDEV_Normalized_MFE_GSEA_NER", "STDEV_Normalized_MFE_GSEA_NHEJ", "STDEV_Normalized_MFE_GSEA_FA", "STDEV_Normalized_MFE_GSEA_BER", "STDEV_Normalized_MFE_GSEA_MMR", "STDEV_Normalized_MFE_GSEA_ABC_transporters", "STDEV_Normalized_MFE_GSEA_Ca_channels_activity", "STDEV_Normalized_MFE_GSEA_TNF_signaling", "STDEV_Normalized_MFE_GSEA_Chaperone_protein_folding", "STDEV_Normalized_MFE_GSEA_TGF_BETA", "STDEV_Normalized_MFE_GSEA_BMP", "STDEV_Normalized_MFE_GSEA_NFKB", "STDEV_Normalized_MFE_GSEA_NOTCH", "STDEV_Normalized_MFE_GSEA_HEDGEHOG", "STDEV_Normalized_MFE_GSEA_WNT", "STDEV_Normalized_MFE_GSEA_AMPK", "STDEV_Normalized_MFE_GSEA_EGFR", "STDEV_Normalized_MFE_GSEA_HIPPO", "STDEV_Normalized_MFE_GSEA_HYPOXIA", "STDEV_Normalized_MFE_GSEA_IL6_JAK_STAT3", "STDEV_Normalized_MFE_GSEA_Insulin", "STDEV_Normalized_MFE_GSEA_JAK_STAT", "STDEV_Normalized_MFE_GSEA_MAPK", "STDEV_Normalized_MFE_GSEA_P53", "STDEV_Normalized_MFE_GSEA_PI3K_AKT", "STDEV_Normalized_MFE_GSEA_TOLL_LIKE_RECEPTOR", "STDEV_Normalized_MFE_GSEA_VEGF", "STDEV_Normalized_MFE_GSEA_Stress_response", "STDEV_Normalized_MFE_GSEA_Bacterial_response", "STDEV_Normalized_MFE_GSEA_Antifungal_response", "STDEV_Normalized_MFE_GSEA_Antiviral_response", "STDEV_Normalized_MFE_GSEA_Apoptosis", "STDEV_Normalized_MFE_GSEA_Autophagy", "STDEV_Normalized_MFE_GSEA_Drug_resistance", "STDEV_Normalized_MFE_GSEA_JNK", "STDEV_Normalized_MFE_GSEA_Chromatin_remodelling", "STDEV_Normalized_MFE_GSEA_CELL_CYCLE_G1_S_PHASE_TRANSITION", "STDEV_Normalized_MFE_GSEA_MITOTIC_SPINDLE_ORGANIZATION", "STDEV_Normalized_MFE_GSEA_POSITIVE_REGULATION_OF_CELL_MIGRATION", "STDEV_Normalized_MFE_GSEA_CELL_MIGRATION", "STDEV_Normalized_MFE_GSEA_POSITIVE_REGULATION_OF_ERK1_AND_ERK2_CASCADE", "STDEV_Normalized_MFE_GSEA_RESPONSE_TO_TYPE_I_INTERFERON", "STDEV_Normalized_MFE_GSEA_RESPONSE_TO_INTERFERON_GAMMA", "STDEV_Normalized_MFE_GSEA_GOBP_EMT")

# loop to assign data.frames to the desired variable name
for (i in (seq_along(GSEA_STDEV_Normalized_MFEs.output))) {
  do.call("<-",list(GSEA_STDEV_Normalized_MFEs.output[i], GSEA_STDEV_test_df[[i]]))
}




# Plot the data
df<-data.frame(mean=c(Mean_Normalized_MFE_GSEA_cell_motility, Mean_Normalized_MFE_GSEA_cell_adhesion, Mean_Normalized_MFE_GSEA_chemotaxis, Mean_Normalized_MFE_GSEA_splicing, Mean_Normalized_MFE_GSEA_homologous_recombination, Mean_Normalized_MFE_GSEA_NER, Mean_Normalized_MFE_GSEA_NHEJ, Mean_Normalized_MFE_GSEA_FA, Mean_Normalized_MFE_GSEA_BER, Mean_Normalized_MFE_GSEA_MMR, Mean_Normalized_MFE_GSEA_ABC_transporters, Mean_Normalized_MFE_GSEA_Ca_channels_activity, Mean_Normalized_MFE_GSEA_TNF_signaling, Mean_Normalized_MFE_GSEA_Chaperone_protein_folding, Mean_Normalized_MFE_GSEA_TGF_BETA, Mean_Normalized_MFE_GSEA_BMP, Mean_Normalized_MFE_GSEA_NFKB, Mean_Normalized_MFE_GSEA_NOTCH, Mean_Normalized_MFE_GSEA_HEDGEHOG, Mean_Normalized_MFE_GSEA_WNT, Mean_Normalized_MFE_GSEA_AMPK, Mean_Normalized_MFE_GSEA_EGFR, Mean_Normalized_MFE_GSEA_HIPPO, Mean_Normalized_MFE_GSEA_HYPOXIA, Mean_Normalized_MFE_GSEA_IL6_JAK_STAT3, Mean_Normalized_MFE_GSEA_Insulin, Mean_Normalized_MFE_GSEA_JAK_STAT, Mean_Normalized_MFE_GSEA_MAPK, Mean_Normalized_MFE_GSEA_P53, Mean_Normalized_MFE_GSEA_PI3K_AKT, Mean_Normalized_MFE_GSEA_TOLL_LIKE_RECEPTOR, Mean_Normalized_MFE_GSEA_VEGF, Mean_Normalized_MFE_GSEA_Stress_response, Mean_Normalized_MFE_GSEA_Bacterial_response, Mean_Normalized_MFE_GSEA_Antifungal_response, Mean_Normalized_MFE_GSEA_Antiviral_response, Mean_Normalized_MFE_GSEA_Apoptosis, Mean_Normalized_MFE_GSEA_Autophagy, Mean_Normalized_MFE_GSEA_Drug_resistance, Mean_Normalized_MFE_GSEA_JNK, Mean_Normalized_MFE_GSEA_Chromatin_remodelling, Mean_Normalized_MFE_GSEA_MITOTIC_SPINDLE_ORGANIZATION, Mean_Normalized_MFE_GSEA_POSITIVE_REGULATION_OF_CELL_MIGRATION, Mean_Normalized_MFE_GSEA_CELL_MIGRATION, Mean_Normalized_MFE_GSEA_POSITIVE_REGULATION_OF_ERK1_AND_ERK2_CASCADE, Mean_Normalized_MFE_GSEA_RESPONSE_TO_TYPE_I_INTERFERON, Mean_Normalized_MFE_GSEA_RESPONSE_TO_INTERFERON_GAMMA, Mean_Normalized_MFE_GSEA_GOBP_EMT),
               sd=c(STDEV_Normalized_MFE_GSEA_cell_motility, STDEV_Normalized_MFE_GSEA_cell_adhesion, STDEV_Normalized_MFE_GSEA_chemotaxis, STDEV_Normalized_MFE_GSEA_splicing, STDEV_Normalized_MFE_GSEA_homologous_recombination, STDEV_Normalized_MFE_GSEA_NER, STDEV_Normalized_MFE_GSEA_NHEJ, STDEV_Normalized_MFE_GSEA_FA, STDEV_Normalized_MFE_GSEA_BER, STDEV_Normalized_MFE_GSEA_MMR, STDEV_Normalized_MFE_GSEA_ABC_transporters, STDEV_Normalized_MFE_GSEA_Ca_channels_activity, STDEV_Normalized_MFE_GSEA_TNF_signaling, STDEV_Normalized_MFE_GSEA_Chaperone_protein_folding, STDEV_Normalized_MFE_GSEA_TGF_BETA, STDEV_Normalized_MFE_GSEA_BMP, STDEV_Normalized_MFE_GSEA_NFKB, STDEV_Normalized_MFE_GSEA_NOTCH, STDEV_Normalized_MFE_GSEA_HEDGEHOG, STDEV_Normalized_MFE_GSEA_WNT, STDEV_Normalized_MFE_GSEA_AMPK, STDEV_Normalized_MFE_GSEA_EGFR, STDEV_Normalized_MFE_GSEA_HIPPO, STDEV_Normalized_MFE_GSEA_HYPOXIA, STDEV_Normalized_MFE_GSEA_IL6_JAK_STAT3, STDEV_Normalized_MFE_GSEA_Insulin, STDEV_Normalized_MFE_GSEA_JAK_STAT, STDEV_Normalized_MFE_GSEA_MAPK, STDEV_Normalized_MFE_GSEA_P53, STDEV_Normalized_MFE_GSEA_PI3K_AKT, STDEV_Normalized_MFE_GSEA_TOLL_LIKE_RECEPTOR, STDEV_Normalized_MFE_GSEA_VEGF, STDEV_Normalized_MFE_GSEA_Stress_response, STDEV_Normalized_MFE_GSEA_Bacterial_response, STDEV_Normalized_MFE_GSEA_Antifungal_response, STDEV_Normalized_MFE_GSEA_Antiviral_response, STDEV_Normalized_MFE_GSEA_Apoptosis, STDEV_Normalized_MFE_GSEA_Autophagy, STDEV_Normalized_MFE_GSEA_Drug_resistance, STDEV_Normalized_MFE_GSEA_JNK, STDEV_Normalized_MFE_GSEA_Chromatin_remodelling, STDEV_Normalized_MFE_GSEA_MITOTIC_SPINDLE_ORGANIZATION, STDEV_Normalized_MFE_GSEA_POSITIVE_REGULATION_OF_CELL_MIGRATION, STDEV_Normalized_MFE_GSEA_CELL_MIGRATION, STDEV_Normalized_MFE_GSEA_POSITIVE_REGULATION_OF_ERK1_AND_ERK2_CASCADE, STDEV_Normalized_MFE_GSEA_RESPONSE_TO_TYPE_I_INTERFERON, STDEV_Normalized_MFE_GSEA_RESPONSE_TO_INTERFERON_GAMMA, STDEV_Normalized_MFE_GSEA_GOBP_EMT),
               Category=as.factor(c("cell_motility", "cell_adhesion", "chemotaxis", "splicing", "HR", "NER", "NHEJ", "FA", "BER", "MMR", "ABC_transporters", "Ca_channels_activity", "TNF_signaling", "Chaperone_protein_folding", "TGF_BETA", "BMP", "NFKB", "NOTCH", "HEDGEHOG", "WNT", "AMPK", "EGFR", "HIPPO", "HYPOXIA", "IL6_JAK_STAT3", "Insulin", "JAK_STAT", "MAPK", "P53", "PI3K_AKT", "TOLL_LIKE_RECEPTOR", "VEGF", "Stress_response", "Bacterial_response", "Antifungal_response", "Antiviral_response", "Apoptosis", "Autophagy", "Drug_resistance", "JNK", "Chromatin_remodelling", "MITOTIC_SPINDLE_ORGANIZATION", "POS_REG_CELL_MIGRATION", "CELL_MIGRATION", "POS_REG_ERK1_ERK2", "RESPONSE_TO_TYPE_I_IFN", "RESPONSE_TO_IFN_GAMMA", "GOBP_EMT")))
df


# Calculate the mean of all genes in each category
Mean_all_pathways <- df %>%
  dplyr::summarise(mean_all_pathways = mean(mean))


library(ggplot2)
ggplot(df, aes(x=reorder(Category, -mean), y=mean)) +
  theme_bw() +
  geom_point(position=position_dodge(width = 1), stat="identity",colour='black') +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.4) +
  theme(axis.text.x = element_text(size=12, angle=45,hjust=0.95,vjust=0.98)) +
  theme(axis.text.y = element_text(size=10)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 10))) +
  geom_point(colour = "black", size = 3) +
  geom_hline(yintercept = Mean_all_pathways$mean_all_pathways, color = "black") +
  xlab("") + ylab("Normalized_MFE (kcal/mol)") +
  ggtitle("Normalized_MFE GSEA signaling pathways (all transcript isoforms)") +
  theme(plot.title = element_text(color="black", size=14, face="bold")) +
  coord_flip()







# Generate box plot

# Function to select specific columns from a list of data frames
select_columns_from_dfs <- function(df_list, columns) {
  # Use lapply to apply the select function to each data frame in the list
  selected_dfs <- lapply(df_list, function(df) {
    dplyr::select(df, all_of(columns))
  })
  # Return the list of modified data frames
  return(selected_dfs)
}

# List of data frames
df_list <- list(
  df_GSEA_cell_motility, 
  # df_GSEA_cell_adhesion, 
  # df_GSEA_chemotaxis, 
  # df_GSEA_splicing, 
  df_GSEA_homologous_recombination, 
  df_GSEA_NER, 
  df_GSEA_NHEJ, 
  df_GSEA_FA, 
  df_GSEA_BER, 
  df_GSEA_MMR, 
  # df_GSEA_ABC_transporters, 
  # df_GSEA_Ca_channels_activity, 
  df_GSEA_TNF_signaling, 
  # df_GSEA_Chaperone_protein_folding, 
  df_GSEA_TGF_BETA, 
  df_GSEA_BMP, 
  df_GSEA_NFKB, 
  df_GSEA_NOTCH, 
  df_GSEA_HEDGEHOG, 
  df_GSEA_WNT, 
  df_GSEA_AMPK, 
  df_GSEA_EGFR, 
  df_GSEA_HIPPO, 
  # df_GSEA_HYPOXIA, 
  df_GSEA_IL6_JAK_STAT3, 
  df_GSEA_Insulin, 
  df_GSEA_JAK_STAT, 
  df_GSEA_MAPK, 
  df_GSEA_P53, 
  df_GSEA_PI3K_AKT, 
  df_GSEA_TOLL_LIKE_RECEPTOR, 
  df_GSEA_VEGF, 
  # df_GSEA_Stress_response, 
  # df_GSEA_Bacterial_response, 
  # df_GSEA_Antifungal_response, 
  # df_GSEA_Antiviral_response, 
  # df_GSEA_Apoptosis, 
  # df_GSEA_Autophagy, 
  # df_GSEA_Drug_resistance, 
  df_GSEA_JNK
  # df_GSEA_Chromatin_remodelling, 
  # df_GSEA_CELL_CYCLE_G1_S_PHASE_TRANSITION, 
  # df_GSEA_MITOTIC_SPINDLE_ORGANIZATION, 
  # df_GSEA_POSITIVE_REGULATION_OF_CELL_MIGRATION, 
  # df_GSEA_CELL_MIGRATION, 
  # df_GSEA_POSITIVE_REGULATION_OF_ERK1_AND_ERK2_CASCADE, 
  # df_GSEA_RESPONSE_TO_TYPE_I_INTERFERON, 
  # df_GSEA_RESPONSE_TO_INTERFERON_GAMMA, 
  # df_GSEA_GOBP_EMT
  )

# Columns to select
columns <- c("external_gene_name", "ensembl_transcript_id_version", "Normalized_MFE")

# Apply the function to select columns from all data frames
selected_dfs <- select_columns_from_dfs(df_list, columns)

# Naming the selected data frames for easy identification
names(selected_dfs) <- c(
  "Cell motility", 
  # "cell_adhesion", 
  # "chemotaxis", 
  # "splicing", 
  "Homologous_recombination", 
  "Nucleotide excision repair", 
  "Non-homologous end joining", 
  "Fanconi anemia", 
  "Base excision repair", 
  "Mismatch repair", 
  # "ABC_transporters", 
  # "Ca_channels_activity", 
  "TNF signaling", 
  # "Chaperone_protein_folding", 
  "TGF Bsta signaling", 
  "BMP signaling", 
  "NFKB signaling", 
  "NOTCH signaling", 
  "HEDGEHOG signaling", 
  "WNT signaling", 
  "AMPK signaling", 
  "EGFR signaling", 
  "HIPPO signaling", 
  # "HYPOXIA", 
  "IL6 JAK STAT3 signaling", 
  "Insulin signaling", 
  "JAK STAT signaling", 
  "MAPK signaling", 
  "P53 signaling", 
  "PI3K AKT signaling", 
  "TOLL like recetor signaling", 
  "VEGF signaling", 
  # "Stress_response", 
  # "Bacterial_response", 
  # "Antifungal_response", 
  # "Antiviral_response", 
  # "Apoptosis", 
  # "Autophagy", 
  # "Drug_resistance", 
  "JNK signaling"
  # "Chromatin_remodelling", 
  # "CELL_CYCLE_G1_S_PHASE_TRANSITION", 
  # "MITOTIC_SPINDLE_ORGANIZATION", 
  # "POS_REG_OF_CELL_MIGRATION", 
  # "CELL_MIGRATION", 
  # "POS_REG_OF_ERK1_AND_ERK2_CASCADE", 
  # "RESPONSE_TO_TYPE_I_INTERFERON", 
  # "RESPONSE_TO_INTERFERON_GAMMA", 
  # "GOBP_EMT"
)

# Print the first few rows of each selected data frame for verification
lapply(names(selected_dfs), function(name) {
  cat("Data frame:", name, "\n")
  print(head(selected_dfs[[name]]))
  cat("\n")
})


# Combine the selected data frames into a single data frame with a source column
combine_dfs <- function(df_list, names_list) {
  combined_df <- bind_rows(lapply(1:length(df_list), function(i) {
    df_list[[i]] %>% mutate(source = names_list[i])
  }))
  return(combined_df)
}

# Combine all the selected data frames
combined_df <- combine_dfs(selected_dfs, names(selected_dfs))

# Calculate the mean of Normalized_MFE for each source
mean_df <- combined_df %>%
  group_by(source) %>%
  summarize(Mean_Norm_MFE = mean(Normalized_MFE, na.rm = TRUE))

# Merge the mean values back into the combined data frame
combined_df <- left_join(combined_df, mean_df, by = "source")

# Create a box plot using ggplot2
ggplot(combined_df, aes(x = reorder(source, -Normalized_MFE), y = Normalized_MFE)) +
  geom_boxplot() +
  theme_classic() +
  geom_point(data = mean_df, aes(x = source, y = Mean_Norm_MFE), color = "blue", size = 2, shape = 18) +  # Add a dot for the mean
  geom_hline(yintercept = Mean_all_pathways$mean_all_pathways, color = "black") +
  theme(axis.text.y = element_text(size=10)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  theme(axis.title.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(size = 10)) +
  labs(title = "Normalized_MFE GSEA signaling pathways (all transcript isoforms)",
       x = "",
       y = "Normalized_MFE (kcal/mol)") +
  coord_flip()














# Fig.1E ####
# Open the curated gene sets

source("path_to/All_gene_sets_AC.R")


# Box plot
# Filter the dataset for Error-prone and High-fidelity DNA polymerases
Error_prone_genes_MFE <- dplyr::filter(Table_ANNOTATED_gene_filtered, external_gene_name %in% Error_prone_DNA_polymerases)
High_fidelity_genes_MFE <- dplyr::filter(Table_ANNOTATED_gene_filtered, external_gene_name %in% High_fidelity_polymerases)

# Add a new column to differentiate between the two groups
Error_prone_genes_MFE$Group <- "Error_prone DNA pol."
High_fidelity_genes_MFE$Group <- "High_fidelity DNA pol."

# Combine the two dataframes
combined_MFE <- rbind(Error_prone_genes_MFE, High_fidelity_genes_MFE)

# Calculate the mean of Normalized_MFE for the combined data
mean_normalized_MFE <- mean(combined_MFE$Normalized_MFE)

# Create the box plot
ggplot(combined_MFE, aes(x = Group, y = Normalized_MFE, fill = Group)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 19, size = 5, color = "blue") +  # Add mean points
  geom_hline(yintercept = mean_normalized_MFE, linetype = "dashed", color = "black", size = 1) +  # Add horizontal line for mean
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 20),
        axis.title.x = element_text(size = 20),  # Increase x title font size
        axis.title.y = element_text(size = 20),  # Increase y title font size
        legend.position = "none") +  # Remove the legend
  labs(x = "", y = "Normalized MFE (kcal/mol)", title = "Normalized MFE Error-prone and High-fidelity DNA polymerases (all TX)") +
  scale_fill_manual(values = c("Error_prone DNA pol." = "lightgrey", "High_fidelity DNA pol." = "lightgrey")) +
  theme(legend.position = "none")  # Remove the legend








# The NHEJ and error-prone DNA polymerases have less negative DG values, opposite to high-fidelity DNA repair mechanisms like BER and high-fidelity DNA polymerases.


# Box plot
# Filter the dataset for Error-prone and High-fidelity DNA polymerases
Error_prone_genes_MFE <- dplyr::filter(Table_ANNOTATED_gene_filtered, external_gene_name %in% Error_prone_DNA_polymerases)
High_fidelity_genes_MFE <- dplyr::filter(Table_ANNOTATED_gene_filtered, external_gene_name %in% High_fidelity_polymerases)
Non_holologous_end_joining_genes_MFE <- dplyr::filter(Table_ANNOTATED_gene_filtered, external_gene_name %in% Non_hologous_end_joining)
Base_excision_repair_genes_MFE <- dplyr::filter(Table_ANNOTATED_gene_filtered, external_gene_name %in% Base_excision_repair)





# Add a new column to differentiate between the two groups
Error_prone_genes_MFE$Group <- "EP DNA pol"
High_fidelity_genes_MFE$Group <- "HF DNA pol"
Non_holologous_end_joining_genes_MFE$Group <- "NHEJ"
Base_excision_repair_genes_MFE$Group <- "BER"

# Combine the two dataframes
combined_MFE <- rbind(Error_prone_genes_MFE, High_fidelity_genes_MFE, Non_holologous_end_joining_genes_MFE, Base_excision_repair_genes_MFE)

# Calculate the mean of Normalized_MFE for the combined data
mean_normalized_MFE <- mean(combined_MFE$Normalized_MFE)

# Create the box plot
EP_HF_DNArepair_replication <- ggplot(combined_MFE, aes(x = Group, y = Normalized_MFE, fill = Group)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 19, size = 3, color = "blue") +  # Add mean points
  geom_hline(yintercept = mean_normalized_MFE, linetype = "dashed", color = "black", size = 1) +  # Add horizontal line for mean
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title.x = element_text(size = 20),  # Increase x title font size
        axis.title.y = element_text(size = 10),  # Increase y title font size
        legend.position = "none") +  # Remove the legend
  labs(x = "", y = "Normalized MFE (kcal/mol)", title = "") +
  scale_fill_manual(values = c("EP DNA pol" = "lightgrey", "HF DNA pol" = "lightgrey", "BER" = "lightgrey", "NHEJ" = "lightgrey")) +
  scale_x_discrete(limits = c("NHEJ", "BER", "EP DNA pol", "HF DNA pol")) +
  theme(legend.position = "none")  # Remove the legend

print(EP_HF_DNArepair_replication)



# KS test
ks_test_NHEJ_BER <- ks.test(Non_holologous_end_joining_genes_MFE$Normalized_MFE, Base_excision_repair_genes_MFE$Normalized_MFE)
print(ks_test_NHEJ_BER)

ks_test_EP_HF_pol <- ks.test(Error_prone_genes_MFE$Normalized_MFE, High_fidelity_genes_MFE$Normalized_MFE)
print(ks_test_EP_HF_pol)

# Wilcoxon rank-sum tests (Mann–Whitney U)
wilcox_test_NHEJ_BER <- wilcox.test(na.omit(Non_holologous_end_joining_genes_MFE$Normalized_MFE),
                                    na.omit(Base_excision_repair_genes_MFE$Normalized_MFE))
print(wilcox_test_NHEJ_BER)

wilcox_test_EP_HF_pol <- wilcox.test(na.omit(Error_prone_genes_MFE$Normalized_MFE),
                                     na.omit(High_fidelity_genes_MFE$Normalized_MFE))
print(wilcox_test_EP_HF_pol)
















# Fig.1F ####
# Phylostratigraphic analysis
# The code here below calculates the Ages Index of the entire Human transcriptome.

# The original script used to calculate the Transcriptome Age Index (TAI) for a given set of gene expression data and corresponding phylostratum values is in the following path:
# /Volumes/AC_HD_8TB/Exp.\ 271\ Profiling\ of\ alternative\ splicing\ landscape_Nanopore\ sequencing\ drug\ treatment\ 72h/R_scripts/TAI_72h_calculated_on_gene_counts_AC_GOOD_TO_USE.R
# The script contains the formula:
# # The TAI formula has been published in Tigros et al. PNAS 2017
# tai <- function(ps,gene_exp) {
#   TAI<-(sum(ps*gene_exp)) / (sum(gene_exp))
#   return(TAI)
# }

# The TAI is a measure used in evolutionary biology to assess the "evolutionary age" of a transcriptome based on the gene expression levels and their associated phylostratum, which represents the evolutionary origin of each gene.


# Open the file of the Trigos et al. PNAS 2017 paper containing the phylustratum of each gene
Genes_phylostratigraphy_PNAS<-read.delim(file="path_to/PNAS_Genes_phylostratigraphy.csv", sep = ",",
                                         header = T, stringsAsFactors = F) %>%
  dplyr::rename(external_gene_name = GeneID)





# Join the data with the predicted MFE values of the Human transcriptome
Table_plylostratigraphy_MFE <- dplyr::left_join(Genes_phylostratigraphy_PNAS, Table_ANNOTATED_gene_names, by = "external_gene_name")


# Group by phylostrata and phylostrata_name calculate the mean Normalized_MFE of the genes in each phylostratum
Evolution_MFE <- Table_plylostratigraphy_MFE %>%
  dplyr::group_by(Phylostrata_N, Phylostrata_name) %>%
  dplyr::summarize(AVG_X = mean(Normalized_MFE, na.rm = T), STDEV_X = sd(Normalized_MFE, na.rm = T)) %>%
  ungroup()

# Ensure Phylostrata_name is a factor ordered by Phylostrata_N
Evolution_MFE <- Evolution_MFE %>%
  arrange(Phylostrata_N) %>%
  mutate(Phylostrata_name = factor(Phylostrata_name, levels = unique(Phylostrata_name)))

# Plot the Phylostrata_N vs the mean Normalized_MFE of the genes in each phylostratum
ggplot(Evolution_MFE, aes(x = Phylostrata_name, y = AVG_X)) +
  geom_point() +
  geom_errorbar(aes(ymin = AVG_X - STDEV_X, ymax = AVG_X + STDEV_X), width = 0.1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "", y = "Mean_Normalized_MFE (kcal/mol)", title = "Phylostrata vs Normalized_MFE") +
  coord_flip()



# Generate box plot
# Ensure Phylostrata_name is a factor ordered by Phylostrata_N
Table_plylostratigraphy_MFE_box_plot <- Table_plylostratigraphy_MFE %>%
  arrange(Phylostrata_N) %>%
  mutate(Phylostrata_name = factor(Phylostrata_name, levels = unique(Phylostrata_name)))

# Plot the Phylostrata_N vs the Normalized_MFE of the genes in each phylostratum
ggplot(Table_plylostratigraphy_MFE_box_plot, aes(x = Phylostrata_name, y = Normalized_MFE)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) +
  labs(x = "", y = "Normalized_MFE (kcal/mol)", title = "Phylostrata vs Normalized_MFE") +
  coord_flip()









# Plot only the phylostrata up to Metazoa
# Calculate the mean of the Cellular_organisms phylostrata
mean_cellular_organisms <- mean(Table_plylostratigraphy_MFE_box_plot$Normalized_MFE[Table_plylostratigraphy_MFE_box_plot$Phylostrata_name == "Cellular_organisms"], na.rm = TRUE)


MFE_phylostratigraphy <- ggplot(Table_plylostratigraphy_MFE_box_plot, aes(x = Phylostrata_name, y = Normalized_MFE)) +
  geom_boxplot(fill = "lightgrey") +
  stat_summary(fun = mean, geom = "point", shape = 19, size = 2, color = "blue") +  # Add mean points
  geom_hline(yintercept = mean_cellular_organisms, linetype = "dashed", color = "black", size = 1) +  # Add horizontal line
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "", y = "Normalized_MFE (kcal/mol)", title = "Phylostrata vs Normalized_MFE") +
  scale_x_discrete(limits = c("Cellular_organisms", "Eukaryota", "Opisthokonta", "Metazoa"))

print(MFE_phylostratigraphy)




# KS test
MFE_Cellular_organisms <- dplyr::filter(Table_plylostratigraphy_MFE, Phylostrata_name == "Cellular_organisms")
MFE_Eukaryota <- dplyr::filter(Table_plylostratigraphy_MFE, Phylostrata_name == "Eukaryota")
MFE_Opisthokonta <- dplyr::filter(Table_plylostratigraphy_MFE, Phylostrata_name == "Opisthokonta")
MFE_Metazoa <- dplyr::filter(Table_plylostratigraphy_MFE, Phylostrata_name == "Metazoa")



ks_test_Eu_vs_CO <- ks.test(MFE_Eukaryota$Normalized_MFE, MFE_Cellular_organisms$Normalized_MFE)
print(ks_test_Eu_vs_CO)

ks_test_Opi_vs_CO <- ks.test(MFE_Opisthokonta$Normalized_MFE, MFE_Cellular_organisms$Normalized_MFE)
print(ks_test_Opi_vs_CO)

ks_test_Meta_vs_CO <- ks.test(MFE_Metazoa$Normalized_MFE, MFE_Cellular_organisms$Normalized_MFE)
print(ks_test_Meta_vs_CO)



# Perform ANOVA on the selected phylostrata
# Filter the data frame for the specified phylostrata names
Filtered_phylostrata_MFE <- Table_plylostratigraphy_MFE %>%
  filter(Phylostrata_name %in% c("Cellular_organisms", "Eukaryota", "Opisthokonta", "Metazoa"))



anova_result <- aov(Normalized_MFE ~ Phylostrata_name, data = Filtered_phylostrata_MFE)

# Get the summary of the ANOVA
anova_summary <- summary(anova_result)
print(anova_summary)

# Extract the p-value
p_value <- anova_summary[[1]]$`Pr(>F)`[1]
print(p_value)


# Perform Linear Trend Analysis
linear_model <- lm(Normalized_MFE ~ Phylostrata_name, data = Filtered_phylostrata_MFE)
summary(linear_model)



