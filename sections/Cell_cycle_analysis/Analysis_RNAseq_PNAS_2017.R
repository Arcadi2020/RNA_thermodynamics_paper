
# This script has been created to identify the genes expressed in each phase of the cell cycle, identified in the paper: Yin Liu et al. PNAS 2017 'Transcriptional landscape of the human cell cycle'
# This script opens the .bw (BigWig) RNAseq files in the NCBI Gene Expression Omnibus depository (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94479)

# IMPORTANT
# These data have been also analysed in the GREIN database GSE94479, but there are only three samples (G0/G1, S and M phases) and no replicates, therefore it is not possible to perform differential expression analysis using the DESeq2 package.



library(rtracklayer)
library(tidyverse)
library(GenomicRanges)
library(derfinder)
library(dplyr)
library(biomaRt)




# Identify the genes that are expressed in each phase of the cell cycle in this cell line.
# Open the RNAseq data of the G0/G1 cell cycle phase.
# These data have been downloaded from Gene Expression Omnibus depository GSE94479.
RNAseq_G0_G1 <- "path_to/GSM2476255_RNA_G0G1.bw"
gr_G0_G1 <- import(RNAseq_G0_G1)

# Convert GRanges object to BED file
# gr_G0_G1_BED <- export.bed(gr_G0_G1, con='granges.bed')
# bed <- as.data.frame(read.table("/Volumes/AC_HD_8TB/Exp.\ 277\ Analysis\ of\ Gibbs\ free\ energy\ in\ response\ to\ drug\ selection/R_scripts/granges.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))



# Overlap the BigWig file with the annotations of the human genome
# Download the gene annotations for hg38 (this is the most recent version).
# download.file(url = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.0.ensembl_genomic.gtf.gz",
#               destfile = "~/Downloads/MANE.gtf.gz")
# gtf_38 <- readGFFAsGRanges("~/Downloads/MANE.gtf.gz")

# If I want to restrict the annotations to the genes symbols:
# gtf <- gtf[gtf$type == "gene"]



# IMPORTANT: the authors have mapped the RNAseq reads to to human genome (hg19) with tophat (v2.1.1).
# It is important to use the same reference genome.
# Download the gene annotations for GRCh37 (hg19) from: https://genome.ucsc.edu/cgi-bin/hgTables
# Assembly: Feb. 2009 (GRCh37/hg19)
# Table: refFlat
# Output format: GTF
# Choose an output name and select gzip compressed
gtf_37 <- readGFFAsGRanges("path_to/refFlat_hg19.gz")

# Find overlap between intervals and genes
my.overlaps_G0_G1 <- findOverlaps(query = gr_G0_G1, subject = gtf_37)

my.results_G0_G1 <- data.frame(to = gtf_37$gene_id[to(my.overlaps_G0_G1)],
                         from = gr_G0_G1[from(my.overlaps_G0_G1)])

# Find the gene lists
Gene.list_G0G1 <- unique(my.results_G0_G1$to)
Gene.list_G0G1_df <- data.frame(Gene.list_G0G1)






# Identify genes S phase
RNAseq_S <- "path_to/GSM2476256_RNA_S.bw"
gr_S <- import(RNAseq_S)

# Find overlap between intervals and genes
my.overlaps_S <- findOverlaps(query = gr_S, subject = gtf_37)

my.results_S <- data.frame(to = gtf_37$gene_id[to(my.overlaps_S)],
                         from = gr_S[from(my.overlaps_S)])

# Find the gene lists
Gene.list_S <- unique(my.results_S$to)
Gene.list_S_df <- data.frame(Gene.list_S)






# Identify genes M phase
RNAseq_M <- "path_to/GSM2476257_RNA_M.bw"
gr_M <- import(RNAseq_M)

# Find overlap between intervals and genes
my.overlaps_M <- findOverlaps(query = gr_M, subject = gtf_37)

my.results_M <- data.frame(to = gtf_37$gene_id[to(my.overlaps_M)],
                         from = gr_M[from(my.overlaps_M)])

# Find the gene lists
Gene.list_M <- unique(my.results_M$to)
Gene.list_M_df <- data.frame(Gene.list_M)












# Load the differential gene expression data
# These files have been downloaded from:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94479

GSE94479_Deseq_G0G1vsG2M <- read.delim("path_to/GSE94479_Deseq_G2MvsG0G1.txt", header=T, stringsAsFactors=F) %>%
  dplyr::rename(refseq_mrna = id) %>%
  dplyr::mutate(log2FoldChange = log2FoldChange * (-1))
  
GSE94479_Deseq_G2MvsS <- read.delim("path_to/GSE94479_Deseq_G2MvsS.txt", header=T, stringsAsFactors=F) %>%
  dplyr::rename(refseq_mrna = id)

GSE94479_Deseq_SvsG0G1 <- read.delim("path_to/GSE94479_Deseq_SvsG0G1.txt", header=T, stringsAsFactors=F) %>%
  dplyr::rename(refseq_mrna = id)

# Load the normalised GROseq data
GSE94479_GROCountNorm <- read.table("path_to/GSE94479_GROCountNorm.txt", header=T, sep = "") %>%
  dplyr::rename(external_gene_name = id) %>%
  dplyr::rename(refseq_mrna = refseq) %>%
  tibble::rownames_to_column(var = "Genomic_coordinates") 

# Use BiomaRt to convert refseq IDs to Ensembl ID

ensembl = useMart(biomart = "ensembl", dataset="hsapiens_gene_ensembl", host="https://www.ensembl.org")

Ensemble_IDs <- getBM(attributes=c('refseq_mrna',
                                   'ensembl_gene_id',
                                   'ensembl_gene_id_version',
                                   'ensembl_transcript_id',
                                   'ensembl_transcript_id_version',
                                   'transcript_is_canonical'),
               filters = 'refseq_mrna',
               values = GSE94479_GROCountNorm$refseq_mrna, # Here goes the list of genes
               mart = ensembl)

GSE94479_GROCountNorm_IDs <- dplyr::left_join(GSE94479_GROCountNorm, Ensemble_IDs, by = "refseq_mrna")




# Add a column with the Ensembl IDs and the gene symbols to the differential expression data
GSE94479_Deseq_G0G1vsG2M_IDs <- dplyr::left_join(GSE94479_Deseq_G0G1vsG2M, Ensemble_IDs, by = "refseq_mrna", multiple = "all")
GSE94479_Deseq_G2MvsS_IDs <- dplyr::left_join(GSE94479_Deseq_G2MvsS, Ensemble_IDs, by = "refseq_mrna", multiple = "all")
GSE94479_Deseq_SvsG0G1_IDs <- dplyr::left_join(GSE94479_Deseq_SvsG0G1, Ensemble_IDs, by = "refseq_mrna", multiple = "all")





