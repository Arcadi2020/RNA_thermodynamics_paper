library(rtracklayer)
library(tidyverse)
library(GenomicRanges)
library(derfinder)
library(dplyr)
library(biomaRt)

RNAseq_G0_G1 <- "../../data/upstream/GSE94479/GSM2476255_RNA_G0G1.bw"
gr_G0_G1 <- import(RNAseq_G0_G1)

gtf_37 <- readGFFAsGRanges("../../data/refFlat_hg19_gtf.gz")

my.overlaps_G0_G1 <- findOverlaps(query = gr_G0_G1, subject = gtf_37)

my.results_G0_G1 <- data.frame(to = gtf_37$gene_id[to(my.overlaps_G0_G1)],
                         from = gr_G0_G1[from(my.overlaps_G0_G1)])

Gene.list_G0G1 <- unique(my.results_G0_G1$to)
Gene.list_G0G1_df <- data.frame(Gene.list_G0G1)

RNAseq_S <- "../../data/upstream/GSE94479/GSM2476256_RNA_S.bw"
gr_S <- import(RNAseq_S)

my.overlaps_S <- findOverlaps(query = gr_S, subject = gtf_37)

my.results_S <- data.frame(to = gtf_37$gene_id[to(my.overlaps_S)],
                         from = gr_S[from(my.overlaps_S)])

Gene.list_S <- unique(my.results_S$to)
Gene.list_S_df <- data.frame(Gene.list_S)

RNAseq_M <- "../../data/upstream/GSE94479/GSM2476257_RNA_M.bw"
gr_M <- import(RNAseq_M)

my.overlaps_M <- findOverlaps(query = gr_M, subject = gtf_37)

my.results_M <- data.frame(to = gtf_37$gene_id[to(my.overlaps_M)],
                         from = gr_M[from(my.overlaps_M)])

Gene.list_M <- unique(my.results_M$to)
Gene.list_M_df <- data.frame(Gene.list_M)

GSE94479_Deseq_G0G1vsG2M <- read.delim("../../data/upstream/GSE94479/GSE94479_Deseq_G2MvsG0G1.txt", header=T, stringsAsFactors=F) %>%
  dplyr::rename(refseq_mrna = id) %>%
  dplyr::mutate(log2FoldChange = log2FoldChange * (-1))

GSE94479_Deseq_G2MvsS <- read.delim("../../data/upstream/GSE94479/GSE94479_Deseq_G2MvsS.txt", header=T, stringsAsFactors=F) %>%
  dplyr::rename(refseq_mrna = id)

GSE94479_Deseq_SvsG0G1 <- read.delim("../../data/upstream/GSE94479/GSE94479_Deseq_SvsG0G1.txt", header=T, stringsAsFactors=F) %>%
  dplyr::rename(refseq_mrna = id)

GSE94479_GROCountNorm <- read.table("../../data/upstream/GSE94479/GSE94479_GROCountNorm.txt", header=T, sep = "") %>%
  dplyr::rename(external_gene_name = id) %>%
  dplyr::rename(refseq_mrna = refseq) %>%
  tibble::rownames_to_column(var = "Genomic_coordinates")

ensembl = useMart(biomart = "ensembl", dataset="hsapiens_gene_ensembl", host="https://www.ensembl.org")

Ensemble_IDs <- getBM(attributes=c('refseq_mrna',
                                   'ensembl_gene_id',
                                   'ensembl_gene_id_version',
                                   'ensembl_transcript_id',
                                   'ensembl_transcript_id_version',
                                   'transcript_is_canonical'),
               filters = 'refseq_mrna',
               values = GSE94479_GROCountNorm$refseq_mrna,
               mart = ensembl)

GSE94479_GROCountNorm_IDs <- dplyr::left_join(GSE94479_GROCountNorm, Ensemble_IDs, by = "refseq_mrna")

GSE94479_Deseq_G0G1vsG2M_IDs <- dplyr::left_join(GSE94479_Deseq_G0G1vsG2M, Ensemble_IDs, by = "refseq_mrna", multiple = "all")
GSE94479_Deseq_G2MvsS_IDs <- dplyr::left_join(GSE94479_Deseq_G2MvsS, Ensemble_IDs, by = "refseq_mrna", multiple = "all")
GSE94479_Deseq_SvsG0G1_IDs <- dplyr::left_join(GSE94479_Deseq_SvsG0G1, Ensemble_IDs, by = "refseq_mrna", multiple = "all")
