suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(biomaRt)
  library(ggplot2)
  library(clusterProfiler)
  library(pheatmap)
})

annotation_file <- "path_to/Table_ANNOTATED_iLINCs.csv"
signature_directory <- "path_to/ILINCS_signatures_3.33uM"
signature_metadata_file <- "path_to/siginfo_beta_LINCS_metadata.csv"
cell_metadata_file <- "path_to/cellinfo_beta.csv"
gene_sets_file <- "path_to/All_gene_sets_AC.R"
analysis_output_dir <- "path_to/iLINCS_3.33uM_analysis_outputs"
intermediate_output_dir <- file.path(analysis_output_dir, "intermediate")
gsea_output_dir <- file.path(analysis_output_dir, "GSEA")
analysis_plot_dir <- file.path(analysis_output_dir, "plots")
final_results_file <- file.path(
  analysis_output_dir,
  "Results_all_cancer_cell_lines_3.33uM.csv"
)
zweighted_dgfold_results_file <- file.path(
  analysis_output_dir,
  "iLINCS_3.33uM_z_score_weighted_DGfold_results.csv"
)
figure4_results_file <- file.path(
  analysis_output_dir,
  "iLINCS_3.33uM_results_for_Fig_4.rds"
)

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
}

for (input_file in c(
  annotation_file,
  signature_metadata_file,
  cell_metadata_file,
  gene_sets_file
)) {
  stop_if_missing(input_file, "Input file")
}

if (!dir.exists(signature_directory)) {
  stop("Signature directory not found: ", signature_directory, call. = FALSE)
}

for (output_directory in c(
  analysis_output_dir,
  intermediate_output_dir,
  gsea_output_dir,
  analysis_plot_dir
)) {
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
}

Table_ANNOTATED_iLINCs <- read.delim(annotation_file,
                                     sep = ",", header = T, stringsAsFactors = F)
Table_ANNOTATED_iLINCs <- dplyr::rename(Table_ANNOTATED_iLINCs, ensembl_id = ensembl_gene_id)

process_gene_iLINCS <- function(iLINCS) {
  iLINCS %>%
    dplyr::filter(gene_id != "117153") %>%
    dplyr::filter(gene_id != "222161") %>%
    dplyr::filter(gene_symbol != "LOC51145") %>%
    dplyr::filter(gene_symbol != "IPW") %>%
    dplyr::filter(gene_symbol != "ZNF205-AS1") %>%
    dplyr::filter(gene_symbol != "RPS14P3") %>%
    dplyr::filter(gene_symbol != "LOC100507424") %>%
    dplyr::filter(gene_symbol != "LILRA3") %>%
    dplyr::filter(gene_symbol != "LOC202181") %>%
    dplyr::filter(gene_symbol != "OCLM") %>%
    dplyr::filter(gene_symbol != "TMEM257") %>%
    dplyr::filter(gene_symbol != "USP9Y") %>%
    dplyr::filter(gene_symbol != "C3orf36") %>%
    dplyr::filter(gene_symbol != "LINC00341") %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "80111" & ensembl_id == "ENSG00000174640", "ENSG00000288547 ", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "304" & gene_symbol == "ANXA2P2", "ENSG00000231991", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "1720" & gene_symbol == "DHFRP3", "ENSG00000168129", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "23117" & gene_symbol == "NPIPB3", "ENSG00000169246", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "56969" & gene_symbol == "RPL23AP32", "ENSG00000237887", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "118433" & gene_symbol == "RPL23AP7", "ENSG00000240356", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "653238" & gene_symbol == "GTF2H2B", "ENSG00000226259", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "386627" & gene_symbol == "SAP30L-AS1", "ENSG00000245275", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "445347" & gene_symbol == "TARP", "ENSG00000289746", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "27004" & gene_symbol == "TCL6", "ENSG00000187621", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "339044" & gene_symbol == "PKD1P1", "ENSG00000244257", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "29798" & gene_symbol == "C2orf27A", "ENSG00000287151", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "643224" & gene_symbol == "TUBBP5", "ENSG00000159247", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "9103" & gene_symbol == "FCGR2C", "ENSG00000244682", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "150244" & gene_symbol == "ZDHHC8P1", "ENSG00000133519", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "6315" & gene_symbol == "ATXN8OS", "ENSG00000230223", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "100131825" & gene_symbol == "CADM3-AS1", "ENSG00000225670", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "5383" & gene_symbol == "PMS2P5", "ENSG00000123965", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "23629" & gene_symbol == "BRD7P3", "ENSG00000169075", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "79946" & gene_symbol == "C10orf95", "ENSG00000120055", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "10151" & gene_symbol == "HNRNPA3P1", "ENSG00000226790", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "26583" & gene_symbol == "DUX4L8", "ENSG00000281720", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "100132923" & gene_symbol == "FAM66D", "ENSG00000255052", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "54581" & gene_symbol == "SCAND2P", "ENSG00000176700", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "8857" & gene_symbol == "FCGBP", "ENSG00000275395", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "389715" & gene_symbol == "FAM205BP", "ENSG00000257198", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "440352" & gene_symbol == "SNX29P2", "NSG00000271699", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "80161" & gene_symbol == "ASMTL-AS1", "ENSG00000236017", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "29075" & gene_symbol == "LINC00652", "ENSG00000179935", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "25787" & gene_symbol == "DGCR9", "ENSG00000273032", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "4583" & gene_symbol == "MUC2", "ENSG00000198788", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "5382" & gene_symbol == "PMS2P4", "ENSG00000067601", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "100129460" & gene_symbol == "DPY19L1P1", "ENSG00000229358", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "80099" & gene_symbol == "C7orf69", "ENSG00000136275", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "440915" & gene_symbol == "POTEKP", "ENSG00000204434", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "55073" & gene_symbol == "LRRC37A4P", "ENSG00000214425", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "79940" & gene_symbol == "LINC00472", "ENSG00000233237", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "10740" & gene_symbol == "RFPL1S", "ENSG00000225465", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "3044" & gene_symbol == "HBBP1", "ENSG00000229988", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "728855" & gene_symbol == "LINC00623", "ENSG00000226067", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "7754" & gene_symbol == "ZNF204P", "ENSG00000204789", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "2952" & gene_symbol == "GSTT1", "ENSG00000277656", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "219595" & gene_symbol == "FOLH1B", "ENSG00000134612", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "5430" & gene_symbol == "POLR2A", "ENSG00000181222", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "6648" & gene_symbol == "SOD2", "ENSG00000291237", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "64595" & gene_symbol == "TTTY15", "ENSG00000114374", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "84726" & gene_symbol == "PRRC2B", "ENSG00000288701", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "8337" & gene_symbol == "HIST2H2AA3", "ENSG00000288825", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "26148" & gene_symbol == "C10orf12", "ENSG00000196233", ensembl_id)) %>%
    dplyr::mutate(gene_id = ifelse(gene_symbol == "C10orf12", "84458", gene_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "56260" & gene_symbol == "C8orf44", "ENSG00000288596", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "440352" & gene_symbol == "SNX29P2", "ENSG00000271699", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "1038" & gene_symbol == "CDR1", "ENSG00000288642", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "79686" & gene_symbol == "LINC00341", "ENSG00000176438 ", ensembl_id)) %>%
    dplyr::mutate(ensembl_id = ifelse(gene_id == "9997" & gene_symbol == "SCO2", "ENSG00000284194 ", ensembl_id)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "AARS", "AARS1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "ACPP", "ACP3", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "C19orf54", "ACTMAP", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "C20orf27", "ADISSP", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "ADSS", "ADSS2", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "SPATA5L1", "AFG2B", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "C1orf109", "AIRIM", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "HHLA3", "ANKRD13C-DT", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "FAM173A", "ANTKMT", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "C18orf25", "ARK2N", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "ARSE", "ARSL", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "FAM214B", "ATOSB", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "ATP5MPL", "ATP5MJ", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "C9orf16", "BBLN", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "KIAA1109", "BLTP1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "KIAA0100", "BLTP2", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "UHRF1BP1L", "BLTP3B", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "ARNTL", "BMAL1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "ARNTL2", "BMAL2", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "C16orf45", "BMERB1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "IMPAD1", "BPNT2", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "C19orf57", "BRME1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "KIAA1107", "BTBD8", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "LINC01136", "BTG2-DT", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "C22orf46", "C22orf46P", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "CARS", "CARS1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "CNTD2", "CCNP", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "CENPBD1P1", "CENPBD2P", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "C3orf14", "CEP15", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "FGFR1OP", "CEP43", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "COL4A3BP", "CERT1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "C11orf1", "CFAP68", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "ICK", "CILK1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "EFCAB1", "CLXN", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "GPR1", "CMKLR2", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "KIAA1211L", "CRACDL", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "C11orf49", "CSTPP1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "FAM49A", "CYRIA", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "FAM49B", "CYRIB", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "DARS", "DARS1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "DEC1", "DELEC1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "ST5", "DENND2B", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "DGCR9", "DGCR5", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "LRRC6","DNAAF11", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "C16orf71","DNAAF8", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "WDR78","DNAI4", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "CASC1","DNAI7", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "DUSP13","DUSP13B", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "WDR60","DYNC2I1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "FAM102A","EEIG1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "KIAA1324","ELAPOR1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "FAM189A2","ENTREP1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "FAM189A1","ENTREP2", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "FAM189B","ENTREP3", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "CXorf40A","EOLA1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "CXorf40B","EOLA2", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "LINC00894","EOLA2-DT", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "EPRS","EPRS1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "C19orf24","FAM174C", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "ZBED8","FAM200C", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "FAM86C1","FAM86C1P", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "CDRT1","FBXW10B", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "FAM160B2","FHIP2B", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "C1orf112","FIRRM", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "G6PC","G6PC1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "KIAA0355","GARRE1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "GARS","GARS1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "GATD3A","GATD3", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "WRB","GET1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "ASNA1","GET3", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "TSTA3","GFUS", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "GK3P","GK3", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "BHLHB9","GPRASP3", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "H1F0","H1-0", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "HIST1H1A","H1-1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "H1FX","H1-10", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "HIST1H1C","H1-2", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "HIST1H1D","H1-3", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "HIST1H1E","H1-4", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "HIST1H1B","H1-5", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "HIST1H1T","H1-6", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "HIST1H2AI","H2AC13", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "HIST1H2AJ","H2AC14", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "HIST2H2AA3","H2AC18", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "HIST3H2A","H2AC25", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "HIST1H2AC","H2AC6", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "HIST1H2AE","H2AC8", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "H2AFJ","H2AJ", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "HYPM","H2AP", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "H2AFX","H2AX", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "H2AFZ","H2AZ1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "H2AFV","H2AZ2", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "HIST1H2BJ","H2BC11", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "HIST1H2BK","H2BC12", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "HIST1H2BL","H2BC13", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "HIST1H2BM","H2BC14", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "HIST1H2BN","H2BC15", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "HIST1H2BO","H2BC17", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "HIST2H2BE","H2BC21", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "HIST1H2BB","H2BC3", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "HIST1H2BD","H2BC5", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "HIST1H2BG","H2BC8", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "HIST1H2BH","H2BC9", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "H3F3B","H3-3B", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "HIST3H3","H3-4", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "HIST1H3A","H3C1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "HIST1H4A","H4C1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "HIST1H4G","H4C7", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "HARS","HARS1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "C17orf53","HROB", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "IARS","IARS1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "HSPB11","IFT25", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "TTC26","IFT56", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "TTC30A","IFT70A", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "C7orf26","INTS15", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "LRMP","IRAG2", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "KARS","KARS1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "KIAA0556","KATNIP", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "KIF1BP","KIFBP", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "C1orf68","KPLCE", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "LARS","LARS1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "C10orf12","LCOR", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "TMEM159","LDAF1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "LOR","LORICRIN", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "C5orf30","MACIR", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "H2AFY","MACROH2A1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "H2AFY2","MACROH2A2", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "MARCH1","MARCHF1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "MARCH2","MARCHF2", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "MARCH3","MARCHF3", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "MARCH5","MARCHF5", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "MARCH6","MARCHF6", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "MARCH7","MARCHF7", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "MARCH8","MARCHF8", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "MARS","MARS1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "KIAA0895","MATCAP2", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "EEF1AKNMT","METTL13", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "RRNAD1","METTL25B", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "LINC00599","MIR124-1HG", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "C1orf61","MIR9-1HG", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "MARC1","MTARC1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "MARC2","MTARC2", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "C12orf10","MYG1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "MYLPF","MYL11", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "FAM155A","NALF1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "FAM155B","NALF2", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "NARS","NARS1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "SLC9A3R1","NHERF1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "SLC9A3R2","NHERF2", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "PDZD3","NHERF4", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "FAM129A","NIBAN1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "WDYHV1","NTAQ1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "NUPL2","NUP42", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "RNF219","OBI1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "ARMC4","ODAD2", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "TMEM206","PACC1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "MPP5","PALS1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "MPP6","PALS2", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "TMEM246","PGAP4", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "TMEM8A","PGAP6", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "C16orf70","PHAF1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "PHB","PHB1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "C9orf116","PIERCE1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "C7orf69","PKD1L1-AS1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "TWISTNB","POLR1F", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "ZNRD1ASP","POLR1HASP", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "KIAA0391","PRORP", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "FAM192A","PSME3IP1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "QARS","QARS1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "RARS","RARS1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "DDX58","RIGI", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "C12orf29","RLIG1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "C16orf58","RUSF1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "SARS","SARS1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "SEPT10","SEPTIN10", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "SEPT11","SEPTIN11", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "SEPT2","SEPTIN2", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "SEPT4","SEPTIN4", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "SEPT6","SEPTIN6", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "SEPT7","SEPTIN7", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "LOC441601","SEPTIN7P11", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "SEPT8","SEPTIN8", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "SEPT9","SEPTIN9", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "C19orf66","SHFL", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "SKIV2L","SKIC2", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "TTC37","SKIC3", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "WDR61","SKIC8", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "SLC26A10","SLC26A10P", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "PQLC2","SLC66A1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "PQLC1","SLC66A2", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "PQLC3","SLC66A3", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "C6orf48","SNHG32", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "FAM205BP","SPATA31F2P", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "THEG","SPMAP2", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "C12orf49","SPRING1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "CXorf56","STEEP1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "TAZ","TAFAZZIN", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "TARS","TARS1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "CXorf21","TASL", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "TCP10","TCP10L3", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "FAM57A","TLCD3A", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "METTL7A","TMT1A", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "MAP11","TRAPPC14", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "TTTY15","USP9Y", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "VARS","VARS1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "WARS","WARS1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "YARS","YARS1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "CCDC130","YJU2B", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "C22orf24","YWHAH-AS1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "ZDHHC8P1","ZDHHC8BP", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "CYHR1","ZFTRAF1", gene_symbol)) %>%
    dplyr::mutate(gene_symbol = ifelse(gene_symbol == "SSSCA1","ZNRD2", gene_symbol)) %>%
    dplyr::mutate(gene_id = ifelse(gene_symbol == "LINC00341", "161176", gene_id))

}

process_tsv_files <- function(directory_path) {

  tsv_files <- list.files(directory_path, pattern = "\\.tsv$", full.names = TRUE, recursive = TRUE)

  data_frames <- list()

  for (tsv_file in tsv_files) {

    data <- read.delim(tsv_file, header = TRUE, sep = "\t")
    data_frames[[tsv_file]] <- data
  }

  return(data_frames)
}

directory_path <- signature_directory

iLINCs_data_frames <- process_tsv_files(directory_path)

names(iLINCs_data_frames)

data_frames <- iLINCs_data_frames
names(data_frames)

Gents_targeting_non_mitogenic_signaling <- c("obatoclax", "NU-7441", "MLN-4924", "CM03", "bortezomib",
                                             "MG-132", "MLN-2238", "MLN-9708", "XMD-1150", "vorinostat", "YM-201636",
                                             "XMD-1150", "entinostat", "TG-101348", "tamoxifen", "MLN-9708",
                                             "MLN-4924", "MLN-8054")

pattern <- paste(Gents_targeting_non_mitogenic_signaling, collapse = "|")

data_frames <- data_frames[!grepl(pattern, names(data_frames))]

data_frames_names <- unlist(lapply(names(data_frames), function(x) {

  cleaned_name <- gsub("_characteristic_direction.tsv", "", x)

  cleaned_name <- paste0(tail(unlist(strsplit(cleaned_name, "/")), n = 2), collapse = "_")
  return(cleaned_name)
}))

data_frames_names <- gsub("_[^_]*$", "", data_frames_names)
data_frames_names <- gsub("_[^_]*$", "", data_frames_names)
data_frames_names <- gsub("_[^_]*$", "", data_frames_names)

data_frames_names <- gsub(":", "_", data_frames_names)

data_frames_names <- gsub("24h_10uM_", "", data_frames_names)

names(data_frames) <- data_frames_names

process_data_frame <- function(df) {
  iLINCS <- process_gene_iLINCS(df)
  return(iLINCS)
}

result_list <- lapply(data_frames, function(x) process_data_frame(x))

directory <- intermediate_output_dir
saveRDS(result_list, file = file.path(directory, "iLINCs_data_unfiltered.rds"))

result_list <- readRDS(file.path(directory, "iLINCs_data_unfiltered.rds"))

Complet_iLINCs_datasets_10uM <- bind_rows(result_list, .id = "name")

Complet_iLINCs_datasets_10uM$extracted_name <- sub("^[^_]*_([^_]*)_.*$", "\\1", Complet_iLINCs_datasets_10uM$name)

Drugs_iLINCs_10uM <- unique(Complet_iLINCs_datasets_10uM$extracted_name)
print(Drugs_iLINCs_10uM)
length(Drugs_iLINCs_10uM)

Complet_iLINCs_datasets_10uM$cell_line <- sub("^(.*?)_.*$", "\\1", Complet_iLINCs_datasets_10uM$name)

Cell_lines_iLINCs_10uM_complete <- unique(Complet_iLINCs_datasets_10uM$cell_line)
print(Cell_lines_iLINCs_10uM_complete)
length(Cell_lines_iLINCs_10uM_complete)

result_list <- result_list[!grepl("XC", names(result_list))]

Lines_unknown_cancer_types <- c("C42", "HIMG002", "NKDBA", "SKB", "OAW28")
result_list <- result_list[!grepl(paste(Lines_unknown_cancer_types, collapse = "|"), names(result_list))]

Normal_cell_lines <- c("MCF10A", "1HAE", "AALE", "ASC", "CD34", "HA1E", "HEK293", "HEK293T", "HEKTE", "HFL1", "HME", "HPTEC", "HUES3", "HUVEC", "IMR90", "LHSAR", "NAMEC8", "NEU", "NL20", "NPC", "PHH", "SKL", "WA09", "WI38", "HIMG001")

result_list <- result_list[!grepl(paste(Normal_cell_lines, collapse = "|"), names(result_list))]

result_list_2 <- result_list

names(result_list_2) <- gsub("_24h_3.33uM", "", names(result_list_2))

Sig_info_iLINCs_metadata <- read.csv(signature_metadata_file, header = TRUE, sep = ",", stringsAsFactors = FALSE)

Sig_info_iLINCs_metadata <- dplyr::filter(Sig_info_iLINCs_metadata, pert_idose == "3.33 uM") %>%
  dplyr::filter(pert_itime == "24 h")

Sig_info_iLINCs_metadata$sig_id <- gsub("/", "_", Sig_info_iLINCs_metadata$sig_id)

Sig_info_iLINCs_metadata$sig_id <- paste0(Sig_info_iLINCs_metadata$cmap_name, "_", Sig_info_iLINCs_metadata$sig_id)

Sig_info_iLINCs_metadata$sig_id <- paste0(Sig_info_iLINCs_metadata$cell_iname, "_", Sig_info_iLINCs_metadata$sig_id)

not_in_metadata <- setdiff(names(result_list_2), Sig_info_iLINCs_metadata$sig_id)

not_in_result_list <- setdiff(Sig_info_iLINCs_metadata$sig_id, names(result_list_2))
Samples_not_in_metadata <- as.data.frame(not_in_result_list)

result_list_2 <- result_list_2[names(result_list_2) %in% Sig_info_iLINCs_metadata$sig_id]

directory <- intermediate_output_dir
saveRDS(result_list_2, file = file.path(directory, "iLINCs_filtered_samples.rds"))

directory <- intermediate_output_dir
result_list_2 <- readRDS(file.path(directory, "iLINCs_filtered_samples.rds"))

get_entrez_ids_once <- function(result_list_2) {

  all_gene_symbols <- unique(unlist(lapply(result_list_2, function(df) df$gene_symbol)))

  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org")
  res_gene_entrez_id <- getBM(attributes = c('entrezgene_id', 'external_gene_name'),
                              filters = 'external_gene_name',
                              values = all_gene_symbols,
                              mart = mart)

  res_gene_entrez_id <- dplyr::filter(res_gene_entrez_id, !is.na(entrezgene_id)) %>%
    dplyr::rename(gene_symbol = external_gene_name)

  return(res_gene_entrez_id)
}

process_and_save_plots_with_entrez <- function(data_frame, z_score_col, plot_title, output_dir, res_gene_entrez_id) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  file_name <- file.path(output_dir, paste0(plot_title, "_GSEA_KEGG.pdf"))

  data_frame <- dplyr::left_join(data_frame, res_gene_entrez_id, by = "gene_symbol")

  all_values_entrez <- data_frame %>%
    dplyr::select(entrezgene_id, !!rlang::sym(z_score_col)) %>%
    dplyr::filter(!is.na(entrezgene_id), !is.na(!!rlang::sym(z_score_col))) %>%
    dplyr::mutate(entrezgene_id = as.character(entrezgene_id))

  collapsed <- all_values_entrez %>%
    dplyr::group_by(entrezgene_id) %>%
    dplyr::summarise(.z = mean(!!rlang::sym(z_score_col), na.rm = TRUE), .groups = "drop")

  collapsed <- collapsed %>%
    dplyr::filter(is.finite(.z)) %>%
    dplyr::arrange(dplyr::desc(.z))

  stats <- collapsed$.z
  names(stats) <- collapsed$entrezgene_id

  cat("Sample:", plot_title, "\n")
  cat("Genes before merge/filter:", nrow(all_values_entrez), "\n")
  cat("Unique Entrez after collapse:", length(stats), "\n")
  if (length(stats) == 0) {
    cat("No genes available for", plot_title, "\n")
    return(NULL)
  }
  if (any(duplicated(names(stats)))) {
    stop("Still have duplicated Entrez IDs after collapsing—please inspect your input.")
  }

  kk2 <- clusterProfiler::gseKEGG(
    geneList = stats,
    organism = "hsa",
    keyType = "ncbi-geneid",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    verbose = TRUE
  )

  if (is.null(kk2) || nrow(kk2) == 0) {
    cat("No significant pathways found for", plot_title, "\n")
    return(NULL)
  }

  p <- enrichplot::dotplot(
    kk2,
    x = "GeneRatio",
    color = "p.adjust",
    showCategory = 50,
    size = NULL,
    split = ".sign",
    font.size = 8,
    title = plot_title,
    orderBy = "x",
    label_format = 100
  ) +
    ggplot2::facet_grid(~.sign) +
    ggplot2::theme(panel.spacing = grid::unit(0.5, "cm"))

  ggplot2::ggsave(file_name, plot = p, width = 8, height = 12)
  return(kk2)
}

process_result_list_with_entrez <- function(result_list_2, output_dir, res_gene_entrez_id) {
  plot_titles <- names(result_list_2)
  kk2_list <- list()

  for (plot_title in plot_titles) {
    data_frame <- result_list_2[[plot_title]]
    z_score_col <- "z_score"

    kk2 <- process_and_save_plots_with_entrez(data_frame, z_score_col, plot_title, output_dir, res_gene_entrez_id)

    if (!is.null(kk2)) {
      kk2_list[[plot_title]] <- kk2
    }
  }

  combined_kk2 <- do.call(rbind, lapply(names(kk2_list), function(name) {
    df <- as.data.frame(kk2_list[[name]])
    df$Sample <- name
    return(df)
  }))

  return(combined_kk2)
}

res_gene_entrez_id <- get_entrez_ids_once(result_list_2)

output_directory <- gsea_output_dir
combined_results <- process_result_list_with_entrez(result_list_2, output_directory, res_gene_entrez_id)

write.csv(combined_results, file = paste0(output_directory, "/combined_kk2_results.csv"), row.names = FALSE)

proliferation_science_genes <- c(
  "TPX2", "PRC1", "FOXM1", "CDK1", "C20orf24", "MCM2", "H2AFZ", "TOP2A", "PCNA", "UBE2C",
  "MELK", "TRIP13", "NCAPD2", "MCM7", "RNASEH2A", "RAD51AP1", "KIF20A", "CDC45", "MAD2L1", "ESPL1",
  "CCNB2", "FEN1", "TTK", "CCT5", "RFC4", "ATAD2", "CKAP5", "NUP205", "CDC20", "CKS2",
  "RRM2", "ELAVL1", "CCNB1", "RRM1", "AURKB", "MSH6", "EZH2", "CTPS", "DKC1", "OIP5",
  "CDCA8", "PTTG1", "CEP55", "H2AFX", "CMAS", "NCAPH", "MCM10", "LSM4", "NCAPG2", "ASF1B",
  "ZWINT", "PBK", "ZWILCH", "CDCA3", "ECT2", "CDC6", "UNG", "MTCH2", "RAD21", "ACTL6A",
  "GPI", "SRSF2", "HDGF", "NXT1", "NEK2", "DHCR7", "AURKA", "NDUFAB1", "TMEM194A", "KIF4A"
)

meta_pcna_genes <- c(
  "PCNA", "CKS2", "NUSAP1", "RRM2", "ZWINT", "PRC1", "TFDP1", "CCNA2", "CCNB1", "MELK",
  "TPX2", "BIRC5", "NCAPG2", "RFWD3", "TACC3", "CDK1", "PCLAF", "MCM2", "MCM7", "MCM5",
  "ASF1B", "CCNB2", "GTPBP2", "KIF20A", "PTTG1", "AURKA", "CDC20", "FOXM1", "GINS2", "MAD2L1",
  "UBE2C", "GINS1", "MCM6", "NCAPD2", "NUP37", "CKS1B", "MCM7", "ARID3A", "AURKB", "CDCA8",
  "SNRPD1", "SNRPD2", "TROAP", "MIS18A", "DDX39A", "AHSP", "ESPL1", "HMBS", "LSM6", "MCM4",
  "NCAPD3", "PPIH", "CDC45", "DTL", "EPB42", "HMGN2", "MCM3", "PF4", "PPBP", "RFC4",
  "TOP2A", "BUB1B", "CDT1", "FEN1", "GATA1", "GYPB", "HMGB2", "KIF22", "KLF1", "CENPU",
  "RFC3", "RHCE", "TAL1", "TCF3", "TCF7L1", "ALAS2", "TSPO2", "CDCA3", "CENPA", "CKLF",
  "GTSE1", "NFE2", "OIP5", "SHCBP1", "SNF8", "SNRPB", "SPTA1", "KEL", "KIF2C", "LYL1",
  "PSMD9", "RACGAP1", "RPIA", "TIMELESS", "TRMT5", "TYMS", "VRK1", "FBXO7", "H3-3A", "NUDT1",
  "PLEK", "POLE2", "RHAG", "NSD2", "APOBEC3B", "BPGM", "CHAF1A", "DNAJC9", "FBXO5", "KIF4A",
  "KLF15", "MKI67", "PGD", "RPA3", "SRSF2", "TRIM58", "ADAMTS13", "CDCA4", "CDKN3", "FECH",
  "LBR", "LIG1", "LMNB1", "MICB", "NUP210", "ORC6", "RAD51AP1", "RHD", "SMC4", "GYPA",
  "RPP30", "TRIM10"
)

essential_cell_cycle_genes <- c(
  "ANAPC2", "ATP6V1B2", "ATR", "AURKA", "BIRC5", "BUB3", "CCDC101", "CDC16", "CDC23", "CDC45",
  "CDC5L", "CDC6", "CDK1", "CDT1", "CHAF1A", "CHEK1", "DAD1", "DCTN2", "DDB1", "DDX54",
  "DGCR14", "DHX30", "DONSON", "DTL", "EEF2", "EIF3C", "EP300", "ESPL1", "FAU", "GINS2",
  "HINFP", "INTS3", "KIF11", "KIF23", "NDC80", "NELFB", "NOL9", "NOP9", "NSF", "ORAOV1",
  "ORC3", "PCID2", "PLK1", "POLA1", "POLE2", "PPP1CB", "PPP2CA", "PRMT1", "PRPF19", "PSMA3",
  "PSMB7", "RACGAP1", "RFC2", "RPA1", "RPS13", "RPS2", "RPS6", "RRM1", "RRM2", "SEC61A1",
  "SLC7A6OS", "SMARCA4", "SRSF3", "TARS", "TICRR", "TIMELESS", "TONSL", "TOP2A", "TPX2", "U2AF2",
  "UPF1", "WDR74", "WEE1", "YME1L1", "ZNF830"
)

TERM2GENE <- data.frame(
  term = c(rep("Proliferation_Science_paper", length(proliferation_science_genes)),
           rep("Meta_PCNA", length(meta_pcna_genes)),
           rep("Essential_cell_cycle_Viner_Breuer_Cell_Division_2019", length(essential_cell_cycle_genes))),
  Gene_symbol = c(proliferation_science_genes, meta_pcna_genes, essential_cell_cycle_genes)
)

TERM2GENE <- TERM2GENE %>%
  dplyr::mutate(term = dplyr::recode(term,
                       "Proliferation_Science_paper" = "Signt_1",
                       "Meta_PCNA" = "Signt_2",
                       "Essential_cell_cycle_Viner_Breuer_Cell_Division_2019" = "Signt_3"))

print(TERM2GENE)

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="https://www.ensembl.org")
gene_ids <- (TERM2GENE$Gene_symbol)

res_gene_entrez_id <- getBM(attributes = c('entrezgene_id',
                                           'external_gene_name'),
                            filters = 'external_gene_name',
                            values = gene_ids,
                            mart = mart)
res_gene_entrez_id <- dplyr::filter(res_gene_entrez_id, !is.na(entrezgene_id)) %>%
  dplyr::rename(Gene_symbol = external_gene_name) %>%
  dplyr::distinct(Gene_symbol, .keep_all = TRUE)

TERM2GENE_entrez <- dplyr::left_join(TERM2GENE, res_gene_entrez_id, by = "Gene_symbol") %>%
  dplyr::filter(!is.na(entrezgene_id)) %>%
  dplyr::select(-Gene_symbol)
TERM2GENE <- TERM2GENE_entrez
print(TERM2GENE)

get_entrez_ids_all <- function(result_list_2) {
  all_gene_symbols <- unique(unlist(lapply(result_list_2, function(df) df$gene_symbol)))

  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org")

  res_gene_entrez_id <- getBM(attributes = c('entrezgene_id', 'external_gene_name'),
                              filters = 'external_gene_name',
                              values = all_gene_symbols,
                              mart = mart)

  res_gene_entrez_id <- dplyr::filter(res_gene_entrez_id, !is.na(entrezgene_id)) %>%
    dplyr::rename(gene_symbol = external_gene_name)

  return(res_gene_entrez_id)
}

get_entrez_ids_all <- function(result_list_2) {
  all_gene_symbols <- unique(unlist(lapply(result_list_2, function(df) df$gene_symbol)))

  mart <- useMart(biomart = "ensembl",
                  dataset = "hsapiens_gene_ensembl",
                  host = "https://www.ensembl.org")

  res_gene_entrez_id <- getBM(
    attributes = c("entrezgene_id", "external_gene_name"),
    filters    = "external_gene_name",
    values     = all_gene_symbols,
    mart       = mart
  )

  res_gene_entrez_id %>%
    dplyr::filter(!is.na(entrezgene_id)) %>%
    dplyr::rename(gene_symbol = external_gene_name) %>%
    dplyr::mutate(entrezgene_id = as.character(entrezgene_id)) %>%
    dplyr::distinct(gene_symbol, entrezgene_id)
}

process_and_save_results <- function(data_frame, z_score_col, plot_title, output_dir,
                                     TERM2GENE, res_gene_entrez_id) {

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  df <- dplyr::left_join(data_frame, res_gene_entrez_id, by = "gene_symbol") %>%
    dplyr::select(entrezgene_id, !!rlang::sym(z_score_col)) %>%
    dplyr::filter(!is.na(entrezgene_id), !is.na(!!rlang::sym(z_score_col))) %>%
    dplyr::mutate(entrezgene_id = as.character(entrezgene_id))

  collapsed <- df %>%
    dplyr::group_by(entrezgene_id) %>%
    dplyr::summarise(stat = mean(!!rlang::sym(z_score_col), na.rm = TRUE), .groups = "drop")

  collapsed <- collapsed %>% dplyr::filter(is.finite(stat)) %>% dplyr::arrange(dplyr::desc(stat))

  stats <- collapsed$stat
  names(stats) <- collapsed$entrezgene_id

  cat("Sample:", plot_title, "\n")
  cat("Genes before collapse:", nrow(df), " | Unique Entrez after collapse:", length(stats), "\n")

  if (length(stats) == 0) {
    cat("No genes available for", plot_title, "\n")
    return(NULL)
  }
  if (any(duplicated(names(stats)))) stop("Duplicated Entrez IDs after collapsing.")

  gsea_res <- clusterProfiler::GSEA(
    geneList      = stats,
    TERM2GENE     = TERM2GENE,
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    minGSSize     = 10,
    maxGSSize     = 500,
    verbose       = TRUE
  )

  if (is.null(gsea_res) || nrow(as.data.frame(gsea_res)) == 0) {
    cat("No significant enrichment found for GSEA (TERM2GENE) in", plot_title, "\n")
  }

  list(gsea_res = gsea_res)
}

process_result_list <- function(result_list_2, output_dir, TERM2GENE, res_gene_entrez_id) {
  plot_titles <- names(result_list_2)
  results_list <- list()

  for (plot_title in plot_titles) {
    df <- result_list_2[[plot_title]]
    res <- process_and_save_results(df, "z_score", plot_title, output_dir, TERM2GENE, res_gene_entrez_id)

    if (!is.null(res) && !is.null(res$gsea_res) && nrow(as.data.frame(res$gsea_res)) > 0) {
      results_list[[plot_title]] <- res
    } else {
      cat("No significant GSEA results for:", plot_title, "\n")
    }
  }

  combined_gsea <- do.call(rbind, lapply(names(results_list), function(name) {
    df <- as.data.frame(results_list[[name]]$gsea_res)
    if (!is.null(df) && nrow(df) > 0) df$Sample <- name
    df
  }))

  list(combined_gsea = combined_gsea)
}

res_gene_entrez_id <- get_entrez_ids_all(result_list_2)

output_directory <- gsea_output_dir
results <- process_result_list(result_list_2, output_directory, TERM2GENE, res_gene_entrez_id)

write.csv(results$combined_gsea, file = paste0(output_directory, "/combined_gsea_results.csv"), row.names = FALSE)

GSEA_combined_kk2 <- read.delim(file.path(gsea_output_dir, "combined_kk2_results.csv"),
                                sep = ",", header = T, stringsAsFactors = F)

GSEA_combined_iLINCs <- read.delim(file.path(gsea_output_dir, "combined_gsea_results.csv"),
                                sep = ",", header = T, stringsAsFactors = F)

iLINCs_complete_GSEA <- dplyr::bind_rows(GSEA_combined_iLINCs, GSEA_combined_kk2)

required_ids <- c("hsa04110", "Signt_1", "Signt_2", "Signt_3")

valid_samples <- iLINCs_complete_GSEA %>%
  dplyr::filter(ID %in% required_ids) %>%
  dplyr::group_by(Sample) %>%
  dplyr::filter(all(required_ids %in% ID)) %>%
  dplyr::pull(Sample) %>%
  unique()

iLINCs_complete_GSEA <- iLINCs_complete_GSEA %>%
  filter(Sample %in% valid_samples)

Selected_KEGG_pathways <- c("hsa04150", "hsa04062", "hsa04668", "hsa04657", "hsa04064", "hsa04110", "hsa03030", "hsa03410", "hsa03430", "hsa03420",
                            "hsa03320", "hsa03460", "hsa03440", "hsa04014", "hsa04066", "hsa04218", "hsa04512", "hsa04024", "hsa04071", "hsa04070",
                            "hsa04068", "hsa04310", "hsa03082", "hsa04015", "hsa04625", "hsa04371", "hsa04390", "hsa04010", "hsa04622", "hsa04350",
                            "hsa04152", "hsa04926", "hsa04630", "hsa04072", "hsa04115", "hsa04151", "hsa04910", "hsa04330", "hsa04012", "hsa04022",
                            "hsa04020", "hsa04370", "hsa04210", "hsa03450", "Signt_1", "Signt_2", "Signt_3")

GSEA_combined_kk2_filtered <- dplyr::filter(iLINCs_complete_GSEA, ID %in% Selected_KEGG_pathways)

generate_heatmap <- function(kk2_df) {

  kk2_df <- kk2_df %>%
    dplyr::select(Description, Sample, NES) %>%
    tidyr::spread(Sample, NES)

  kk2_df[is.na(kk2_df)] <- 0

  rownames(kk2_df) <- kk2_df$Description
  kk2_df <- kk2_df[,-1]

  my_palette <- colorRampPalette(c("blue", "white", "red"))(100)

  max_val <- max(kk2_df)
  min_val <- min(kk2_df)

  breaks <- c(seq(min_val, -0.001, length.out = 50), 0, seq(0.001, max_val, length.out = 50))

  pheatmap(as.matrix(kk2_df),
           color = my_palette,
           breaks = breaks,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           display_numbers = TRUE,
           number_format = "%.2g",
           main = "Heatmap of KEGG pathways all samples",
           fontsize_row = 10,
           fontsize_col = 2,
           fontsize_number = 0)
}

generate_heatmap(GSEA_combined_kk2_filtered)

Samples_N1 <- unique(GSEA_combined_kk2_filtered$Sample)

iLINCs_proliferating <- GSEA_combined_kk2_filtered[(GSEA_combined_kk2_filtered$ID %in% c("hsa04110", "Signt_1", "Signt_2", "Signt_3") & GSEA_combined_kk2_filtered$NES >= 0), ]
iLINCs_proliferating <- dplyr::distinct(iLINCs_proliferating, Sample, .keep_all = TRUE)

GSEA_combined_kk2_filtered_cell_cycle_arrest <- GSEA_combined_kk2_filtered[!(GSEA_combined_kk2_filtered$Sample %in% iLINCs_proliferating$Sample), ]

generate_heatmap <- function(kk2_df) {

  kk2_df <- kk2_df %>%
    dplyr::select(Description, Sample, NES) %>%
    tidyr::spread(Sample, NES)

  kk2_df[is.na(kk2_df)] <- 0

  rownames(kk2_df) <- kk2_df$Description
  kk2_df <- kk2_df[,-1]

  my_palette <- colorRampPalette(c("blue", "white", "red"))(100)

  max_val <- max(kk2_df, na.rm = TRUE)
  min_val <- min(kk2_df, na.rm = TRUE)

  breaks <- c(seq(min_val, -0.001, length.out = 50),
              0,
              seq(0.001, max_val, length.out = 50))

  pheatmap(as.matrix(kk2_df),
           color = my_palette,
           breaks = breaks,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           display_numbers = TRUE,
           number_format = "%.2g",
           main = "Heatmap of KEGG pathways samples cell cycle arrest",
           fontsize_row = 10,
           fontsize_col = 2,
           fontsize_number = 0)
}

generate_heatmap(GSEA_combined_kk2_filtered_cell_cycle_arrest)

Samples_N2 <- unique(GSEA_combined_kk2_filtered_cell_cycle_arrest$Sample)

ids_to_check <- c(
  "hsa04150", "hsa04062", "hsa04668", "hsa04657", "hsa04064", "hsa03320",
  "hsa04014", "hsa04066", "hsa04068", "hsa04310", "hsa04015", "hsa04070",
  "hsa04625", "hsa04371", "hsa04390", "hsa04010", "hsa04622", "hsa04350",
  "hsa04152", "hsa04926", "hsa04630", "hsa04072", "hsa04151", "hsa04512",
  "hsa04910", "hsa04330", "hsa04012", "hsa04022", "hsa04020", "hsa03030",
  "hsa03420", "hsa04071"
)

samples_to_exclude <- GSEA_combined_kk2_filtered_cell_cycle_arrest %>%
  dplyr::filter(ID %in% ids_to_check, NES >= 0) %>%
  dplyr::select(Sample) %>%
  dplyr::distinct()

iLINCS_cell_cycle_arrest_mitogenic_signaling <- dplyr::filter(GSEA_combined_kk2_filtered_cell_cycle_arrest,
                                                              !Sample %in% samples_to_exclude$Sample)

generate_heatmap <- function(kk2_df) {

  kk2_df <- kk2_df %>%
    dplyr::select(Description, Sample, NES) %>%
    tidyr::spread(Sample, NES)

  kk2_df[is.na(kk2_df)] <- 0

  rownames(kk2_df) <- kk2_df$Description
  kk2_df <- kk2_df[,-1]

  pheatmap(as.matrix(kk2_df),
           color = colorRampPalette(c("blue", "white"))(50),
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           display_numbers = TRUE,
           number_format = "%.2g",
           main = "Heatmap of KEGG pathways cell cycle arrest and no mitogenic signaling",
           fontsize_row = 11,
           fontsize_col = 2,
           fontsize_number = 0)
}

Heatmap_cell_cycle <- generate_heatmap(iLINCS_cell_cycle_arrest_mitogenic_signaling)
print(Heatmap_cell_cycle)

Samples_N3 <- unique(iLINCS_cell_cycle_arrest_mitogenic_signaling$Sample)

ggsave("Heatmap_cell_cycle.pdf", plot = Heatmap_cell_cycle, width = 10, height = 8, path = analysis_plot_dir)

sample_names_to_retain <- unique(iLINCS_cell_cycle_arrest_mitogenic_signaling$Sample)
result_list_2_filtered <- result_list_2[names(result_list_2) %in% sample_names_to_retain]

source(gene_sets_file)

Table_ANNOTATED_iLINCs <- read.delim(annotation_file,
                                     sep = ",", header = T, stringsAsFactors = F)
Table_ANNOTATED_iLINCs <- dplyr::rename(Table_ANNOTATED_iLINCs, ensembl_id = ensembl_gene_id)

safe_ks <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) < 2 || length(y) < 2) return(list(p.value = NA_real_))
  ks.test(x, y)
}
safe_wilcox <- function(x, y, ...) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) < 1 || length(y) < 1) return(list(p.value = NA_real_))
  suppressWarnings(wilcox.test(x, y, ...))
}

all_results <- data.frame(
  SampleName = character(0),
  Net_MFE_change = numeric(0),
  MeanMFE_UP = numeric(0),
  MeanMFE_DOWN = numeric(0),
  TTEST_pValue = numeric(0),
  wilcox_test_all_UP_vs_DOWN_pValue = numeric(0),
  ks_test_all_UP_vs_DOWN_pValue = numeric(0),
  ks_test_MFE_075_above_pValue = numeric(0),
  ks_test_MFE_025_below_pValue = numeric(0),
  Proliferation_Science = numeric(0),
  PCNA = numeric(0),
  Essential_cell_cycle_genes = numeric(0),
  Homologous_recombination = numeric(0),
  Fanconi_anemia = numeric(0),
  High_fid_pol = numeric(0)
)

plot_dir <- analysis_plot_dir
results_csv <- final_results_file

try(dir.create(plot_dir, recursive = TRUE), silent = TRUE)

process_result_list <- function(Sample_name) {
  message("Processing ", Sample_name)
  df <- result_list_2_filtered[[Sample_name]]

  title_plot_1 <- paste0(Sample_name, "_Analysis MFE 24h p_value < 0.05")
  title_plot_2 <- paste0("Analysis cell proliferation iLINCS 24h_", Sample_name)
  message("Saving plots under ", plot_dir)

  Sample_name_TX <- dplyr::right_join(df, Table_ANNOTATED_iLINCs, by = "ensembl_id")
  Sample_name_TX_filtered <- dplyr::filter(Sample_name_TX, p_value < 0.05)

  MFE_Sample_name_UP   <- dplyr::filter(Sample_name_TX_filtered, z_score > 0)
  MFE_Sample_name_DOWN <- dplyr::filter(Sample_name_TX_filtered, z_score < 0)

  Mean_MFE_Sample_name_UP   <- mean(MFE_Sample_name_UP$Normalized_MFE, na.rm = TRUE)
  Mean_MFE_Sample_name_DOWN <- mean(MFE_Sample_name_DOWN$Normalized_MFE, na.rm = TRUE)

  TTEST <- tryCatch(
    t.test(MFE_Sample_name_DOWN$Normalized_MFE, MFE_Sample_name_UP$Normalized_MFE, paired = FALSE),
    error = function(e) list(p.value = NA_real_)
  )

  wilcox_test_all_UP_vs_DOWN <- safe_wilcox(
    MFE_Sample_name_DOWN$Normalized_MFE,
    MFE_Sample_name_UP$Normalized_MFE,
    paired = FALSE,
    alternative = "two.sided",
    exact = FALSE
  )

  ks_test_all_UP_vs_DOWN <- safe_ks(MFE_Sample_name_DOWN$Normalized_MFE, MFE_Sample_name_UP$Normalized_MFE)

  q75_up   <- tryCatch(quantile(MFE_Sample_name_UP$Normalized_MFE,   probs = 0.75, na.rm = TRUE), error = function(e) NA_real_)
  q75_down <- tryCatch(quantile(MFE_Sample_name_DOWN$Normalized_MFE, probs = 0.75, na.rm = TRUE), error = function(e) NA_real_)
  ks_test_MFE_075_above <- if (is.finite(q75_up) && is.finite(q75_down)) {
    safe_ks(
      MFE_Sample_name_UP$Normalized_MFE[  MFE_Sample_name_UP$Normalized_MFE   >= q75_up],
      MFE_Sample_name_DOWN$Normalized_MFE[MFE_Sample_name_DOWN$Normalized_MFE >= q75_down]
    )
  } else list(p.value = NA_real_)

  q25_up   <- tryCatch(quantile(MFE_Sample_name_UP$Normalized_MFE,   probs = 0.25, na.rm = TRUE), error = function(e) NA_real_)
  q25_down <- tryCatch(quantile(MFE_Sample_name_DOWN$Normalized_MFE, probs = 0.25, na.rm = TRUE), error = function(e) NA_real_)
  ks_test_MFE_025_below <- if (is.finite(q25_up) && is.finite(q25_down)) {
    safe_ks(
      MFE_Sample_name_UP$Normalized_MFE[  MFE_Sample_name_UP$Normalized_MFE   <= q25_up],
      MFE_Sample_name_DOWN$Normalized_MFE[MFE_Sample_name_DOWN$Normalized_MFE <= q25_down]
    )
  } else list(p.value = NA_real_)

  Sample_name_TX_net_MFE_change <- Sample_name_TX_filtered %>%
    dplyr::mutate(Net_MFE_change = (Normalized_MFE * z_score) * (-1))
  Net_MFE_change <- sum(Sample_name_TX_net_MFE_change$Net_MFE_change, na.rm = TRUE)

  MFE_Sample_name_plot <- ggplot(Sample_name_TX_filtered, aes(x = factor(z_score > 0), y = Normalized_MFE)) +
    geom_violin(fill = "lightgrey", alpha = 0.7) +
    geom_jitter(width = 0.1, size = 1, alpha = 0.3) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 20),
      plot.title = element_text(size = 10),
      axis.text.x = element_text(size = 16, color = "black")
    ) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    labs(y = "Normalized_MFE", x = NULL, title = title_plot_1) +
    scale_x_discrete(labels = c("z_score < 0", "z_score > 0"))

  plot_fname_1 <- file.path(plot_dir, paste0("Plot_", Sample_name, "_MFE.pdf"))
  ggsave(filename = plot_fname_1, plot = MFE_Sample_name_plot, width = 20, height = 15, units = "cm")

  filtering <- function(x) {
    dplyr::filter(Sample_name_TX, gene_symbol %in% x)
  }

  CURATED_PATHWAYS.input <- list(
    MTOR_positive_regulators, MTOR_negative_regulators, Proliferation_Science_paper, Meta_PCNA,
    Proliferation_GO_genes_Magill_Sack_2018, Proliferation_STOP_genes_Magill_Sack_2018,
    Reactome_G1_phase, Proliferation_Carey_2021, Essential_cell_cycle_Viner_Breuer_Cell_Division_2019,
    Homologous_recombination, Fanconi_anemia, High_fidelity_polymerases
  )
  CURATED_PATHWAYS.output <- c(
    "MTOR_POS","MTOR_NEG","Proliferation_Science","Meta_PCNA",
    "Proliferation_GO","Proliferation_STOP","G1_phase","Proliferation_Carey",
    "Essential_cell_cycle_Viner_Breuer_Cell_Division_2019","HR","Fanconi_anemia","HF_POL"
  )
  names(CURATED_PATHWAYS.input) <- CURATED_PATHWAYS.output

  Curated_test_df <- lapply(CURATED_PATHWAYS.input, filtering)
  Curated_dfs.output <- c(
    "df_MTOR_POS","df_MTOR_NEG","df_Proliferation_Science","df_Meta_PCNA",
    "df_Proliferation_GO","df_Proliferation_STOP","df_G1_phase","df_Proliferation_Carey",
    "df_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019","df_HR","df_FA","df_HF_POL"
  )
  for (i in seq_along(Curated_dfs.output)) {
    do.call("<-", list(Curated_dfs.output[i], Curated_test_df[[i]]))
  }

  Mean_curated <- function(x) mean(x$z_score, na.rm = TRUE)
  STDEV_curated <- function(x) sd(x$z_score, na.rm = TRUE)

  Curated_input_dfs <- list(
    df_MTOR_POS, df_MTOR_NEG, df_Proliferation_Science, df_Meta_PCNA,
    df_Proliferation_GO, df_Proliferation_STOP, df_G1_phase, df_Proliferation_Carey,
    df_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019, df_HR, df_FA, df_HF_POL
  )

  Curated_mean_test_df  <- lapply(Curated_input_dfs, Mean_curated)
  Curated_stdev_test_df <- lapply(Curated_input_dfs, STDEV_curated)

  names(Curated_mean_test_df)  <- CURATED_PATHWAYS.output
  names(Curated_stdev_test_df) <- CURATED_PATHWAYS.output

  Mean_MTOR_POS  <- Curated_mean_test_df[["MTOR_POS"]]
  Mean_MTOR_NEG  <- Curated_mean_test_df[["MTOR_NEG"]]
  Mean_Proliferation_Science <- Curated_mean_test_df[["Proliferation_Science"]]
  Mean_Meta_PCNA <- Curated_mean_test_df[["Meta_PCNA"]]
  Mean_Proliferation_GO <- Curated_mean_test_df[["Proliferation_GO"]]
  Mean_Proliferation_STOP <- Curated_mean_test_df[["Proliferation_STOP"]]
  Mean_G1_phase <- Curated_mean_test_df[["G1_phase"]]
  Mean_Proliferation_Carey <- Curated_mean_test_df[["Proliferation_Carey"]]
  Mean_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019 <- Curated_mean_test_df[["Essential_cell_cycle_Viner_Breuer_Cell_Division_2019"]]
  Mean_HR <- Curated_mean_test_df[["HR"]]
  Mean_FA <- Curated_mean_test_df[["Fanconi_anemia"]]
  Mean_HF_POL <- Curated_mean_test_df[["HF_POL"]]

  STDEV_MTOR_POS  <- Curated_stdev_test_df[["MTOR_POS"]]
  STDEV_MTOR_NEG  <- Curated_stdev_test_df[["MTOR_NEG"]]
  STDEV_Proliferation_Science <- Curated_stdev_test_df[["Proliferation_Science"]]
  STDEV_Meta_PCNA <- Curated_stdev_test_df[["Meta_PCNA"]]
  STDEV_Proliferation_GO <- Curated_stdev_test_df[["Proliferation_GO"]]
  STDEV_Proliferation_STOP <- Curated_stdev_test_df[["Proliferation_STOP"]]
  STDEV_G1_phase <- Curated_stdev_test_df[["G1_phase"]]
  STDEV_Proliferation_Carey <- Curated_stdev_test_df[["Proliferation_Carey"]]
  STDEV_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019 <- Curated_stdev_test_df[["Essential_cell_cycle_Viner_Breuer_Cell_Division_2019"]]
  STDEV_HR <- Curated_stdev_test_df[["HR"]]
  STDEV_FA <- Curated_stdev_test_df[["Fanconi_anemia"]]
  STDEV_HF_POL <- Curated_stdev_test_df[["HF_POL"]]

  df_names <- c("df_MTOR_POS","df_MTOR_NEG","df_Proliferation_Science","df_Meta_PCNA","df_Proliferation_GO",
                "df_Proliferation_STOP","df_G1_phase","df_Proliferation_Carey",
                "df_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019","df_HR","df_FA","df_HF_POL")
  for (nm in df_names) assign(nm, transform(get(nm), Gene_set = nm))

  iLINCS_gene_sets <- bind_rows(
    df_MTOR_POS, df_MTOR_NEG, df_Proliferation_Science, df_Meta_PCNA, df_Proliferation_GO,
    df_Proliferation_STOP, df_G1_phase, df_Proliferation_Carey,
    df_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019, df_HR, df_FA, df_HF_POL
  )

  Gene_sets_plot <- ggplot(iLINCS_gene_sets, aes(Gene_set, z_score)) +
    geom_boxplot(fill = "lightgrey") +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 0.95, vjust = 0.98),
      plot.title = element_text(size = 9)
    ) +
    labs(x = "", y = "z-score", title = title_plot_2) +
    scale_x_discrete(limits = c("df_Proliferation_Science","df_Meta_PCNA","df_Proliferation_GO","df_Proliferation_STOP",
                                "df_Proliferation_Carey","df_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019",
                                "df_G1_phase","df_MTOR_POS","df_MTOR_NEG","df_HR","df_FA","df_HF_POL"))

  plot_fname_2 <- file.path(plot_dir, paste0("Plot_", Sample_name, "_Gene_sets.pdf"))
  ggsave(filename = plot_fname_2, plot = Gene_sets_plot, width = 20, height = 15, units = "cm")

  sample_results <- data.frame(
    SampleName = Sample_name,
    Net_MFE_change = Net_MFE_change,
    MeanMFE_UP = Mean_MFE_Sample_name_UP,
    MeanMFE_DOWN = Mean_MFE_Sample_name_DOWN,
    TTEST_pValue = TTEST$p.value,
    wilcox_test_all_UP_vs_DOWN_pValue = wilcox_test_all_UP_vs_DOWN$p.value,
    ks_test_all_UP_vs_DOWN_pValue = ks_test_all_UP_vs_DOWN$p.value,
    ks_test_MFE_075_above_pValue = ks_test_MFE_075_above$p.value,
    ks_test_MFE_025_below_pValue = ks_test_MFE_025_below$p.value,
    Proliferation_Science = Mean_Proliferation_Science,
    PCNA = Mean_Meta_PCNA,
    Essential_cell_cycle_genes = Mean_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019,
    Homologous_recombination = Mean_HR,
    Fanconi_anemia = Mean_FA,
    High_fid_pol = Mean_HF_POL
  )

  all_results <<- rbind(all_results, sample_results)
}

Data_result_list <- lapply(names(result_list_2_filtered), function(x) process_result_list(x))

readr::write_csv(all_results, results_csv)

iLINCS_DGfold_annotation <- Table_ANNOTATED_iLINCs %>%
  dplyr::mutate(
    ensembl_id = trimws(ensembl_id)
  ) %>%
  dplyr::select(
    ensembl_id,
    Normalized_MFE
  ) %>%
  dplyr::filter(
    is.finite(Normalized_MFE)
  ) %>%
  dplyr::distinct(
    ensembl_id,
    .keep_all = TRUE
  )

safe_weighted_mean <- function(x, w) {

  x <- as.numeric(x)
  w <- as.numeric(w)

  keep <- is.finite(x) & is.finite(w) & w > 0

  if (sum(keep) == 0) {
    return(NA_real_)
  }

  sum(x[keep] * w[keep], na.rm = TRUE) /
    sum(w[keep], na.rm = TRUE)
}

safe_spearman_p_value <- function(x, y) {

  x <- as.numeric(x)
  y <- as.numeric(y)

  keep <- is.finite(x) & is.finite(y)

  if (sum(keep) < 3) {
    return(NA_real_)
  }

  suppressWarnings(
    cor.test(
      x[keep],
      y[keep],
      method = "spearman",
      exact = FALSE
    )$p.value
  )
}

calculate_z_weighted_DGfold_one_universe <- function(
    sample_tbl,
    Sample_name,
    universe_label,
    p_cutoff = 0.05,
    B_perm = 1000,
    seed = 123
) {

  if (universe_label == "All finite z-score genes") {

    analysis_tbl <- sample_tbl

  } else if (universe_label == "Significant genes p_value < 0.05") {

    analysis_tbl <- sample_tbl %>%
      dplyr::filter(
        is.finite(p_value),
        p_value < p_cutoff
      )

  } else {

    stop("Unknown universe_label.")
  }

  analysis_tbl <- analysis_tbl %>%
    dplyr::filter(
      is.finite(z_score),
      is.finite(Normalized_MFE),
      z_score != 0
    )

  n_genes_used <- dplyr::n_distinct(analysis_tbl$ensembl_id)
  n_positive_z <- dplyr::n_distinct(analysis_tbl$ensembl_id[analysis_tbl$z_score > 0])
  n_negative_z <- dplyr::n_distinct(analysis_tbl$ensembl_id[analysis_tbl$z_score < 0])

  if (n_genes_used < 3 || n_positive_z < 1 || n_negative_z < 1) {

    return(
      tibble::tibble(
        SampleName = Sample_name,
        universe = universe_label,
        n_genes_used = n_genes_used,
        n_positive_z = n_positive_z,
        n_negative_z = n_negative_z,
        z_weighted_MFE_UP = NA_real_,
        z_weighted_MFE_DOWN = NA_real_,
        delta_z_weighted_Normalized_MFE_UP_minus_DOWN = NA_real_,
        z_weighted_centered_DGfold_perturbation_index = NA_real_,
        z_weighted_standardized_DGfold_perturbation_index = NA_real_,
        spearman_z_score_Normalized_MFE = NA_real_,
        spearman_z_score_Normalized_MFE_p_value = NA_real_,
        permutation_p_delta_UP_minus_DOWN = NA_real_,
        permutation_p_centered_index = NA_real_,
        B_perm = B_perm
      )
    )
  }

  MFE_mean <- mean(analysis_tbl$Normalized_MFE, na.rm = TRUE)
  MFE_sd <- sd(analysis_tbl$Normalized_MFE, na.rm = TRUE)

  analysis_tbl <- analysis_tbl %>%
    dplyr::mutate(
      abs_z_score = abs(z_score),
      Normalized_MFE_centered =
        Normalized_MFE - MFE_mean
    )

  if (is.finite(MFE_sd) && MFE_sd > 0) {

    analysis_tbl <- analysis_tbl %>%
      dplyr::mutate(
        Normalized_MFE_standardized =
          (Normalized_MFE - MFE_mean) / MFE_sd
      )

  } else {

    analysis_tbl <- analysis_tbl %>%
      dplyr::mutate(
        Normalized_MFE_standardized = NA_real_
      )
  }

  z_weighted_MFE_UP <- analysis_tbl %>%
    dplyr::filter(
      z_score > 0
    ) %>%
    dplyr::summarise(
      value = safe_weighted_mean(
        x = Normalized_MFE,
        w = abs_z_score
      )
    ) %>%
    dplyr::pull(value)

  z_weighted_MFE_DOWN <- analysis_tbl %>%
    dplyr::filter(
      z_score < 0
    ) %>%
    dplyr::summarise(
      value = safe_weighted_mean(
        x = Normalized_MFE,
        w = abs_z_score
      )
    ) %>%
    dplyr::pull(value)

  observed_delta_UP_minus_DOWN <-
    z_weighted_MFE_UP - z_weighted_MFE_DOWN

  observed_centered_index <- sum(
    analysis_tbl$z_score *
      analysis_tbl$Normalized_MFE_centered,
    na.rm = TRUE
  ) / sum(
    analysis_tbl$abs_z_score,
    na.rm = TRUE
  )

  observed_standardized_index <- sum(
    analysis_tbl$z_score *
      analysis_tbl$Normalized_MFE_standardized,
    na.rm = TRUE
  ) / sum(
    analysis_tbl$abs_z_score,
    na.rm = TRUE
  )

  spearman_cor <- suppressWarnings(
    cor(
      analysis_tbl$z_score,
      analysis_tbl$Normalized_MFE,
      method = "spearman",
      use = "complete.obs"
    )
  )

  spearman_p <- safe_spearman_p_value(
    analysis_tbl$z_score,
    analysis_tbl$Normalized_MFE
  )

  set.seed(seed)

  permutation_delta <- replicate(
    B_perm,
    {
      MFE_perm <- sample(
        analysis_tbl$Normalized_MFE,
        size = nrow(analysis_tbl),
        replace = FALSE
      )

      perm_tbl <- analysis_tbl %>%
        dplyr::mutate(
          Normalized_MFE_perm = MFE_perm
        )

      perm_UP <- perm_tbl %>%
        dplyr::filter(
          z_score > 0
        ) %>%
        dplyr::summarise(
          value = safe_weighted_mean(
            x = Normalized_MFE_perm,
            w = abs_z_score
          )
        ) %>%
        dplyr::pull(value)

      perm_DOWN <- perm_tbl %>%
        dplyr::filter(
          z_score < 0
        ) %>%
        dplyr::summarise(
          value = safe_weighted_mean(
            x = Normalized_MFE_perm,
            w = abs_z_score
          )
        ) %>%
        dplyr::pull(value)

      perm_UP - perm_DOWN
    }
  )

  set.seed(seed)

  permutation_centered_index <- replicate(
    B_perm,
    {
      MFE_perm <- sample(
        analysis_tbl$Normalized_MFE,
        size = nrow(analysis_tbl),
        replace = FALSE
      )

      MFE_perm_centered <- MFE_perm - mean(MFE_perm, na.rm = TRUE)

      sum(
        analysis_tbl$z_score * MFE_perm_centered,
        na.rm = TRUE
      ) / sum(
        analysis_tbl$abs_z_score,
        na.rm = TRUE
      )
    }
  )

  permutation_p_delta <- (
    sum(
      abs(permutation_delta) >= abs(observed_delta_UP_minus_DOWN),
      na.rm = TRUE
    ) + 1
  ) / (B_perm + 1)

  permutation_p_centered_index <- (
    sum(
      abs(permutation_centered_index) >= abs(observed_centered_index),
      na.rm = TRUE
    ) + 1
  ) / (B_perm + 1)

  tibble::tibble(
    SampleName = Sample_name,
    universe = universe_label,
    n_genes_used = n_genes_used,
    n_positive_z = n_positive_z,
    n_negative_z = n_negative_z,
    z_weighted_MFE_UP = z_weighted_MFE_UP,
    z_weighted_MFE_DOWN = z_weighted_MFE_DOWN,
    delta_z_weighted_Normalized_MFE_UP_minus_DOWN =
      observed_delta_UP_minus_DOWN,
    z_weighted_centered_DGfold_perturbation_index =
      observed_centered_index,
    z_weighted_standardized_DGfold_perturbation_index =
      observed_standardized_index,
    spearman_z_score_Normalized_MFE =
      spearman_cor,
    spearman_z_score_Normalized_MFE_p_value =
      spearman_p,
    permutation_p_delta_UP_minus_DOWN =
      permutation_p_delta,
    permutation_p_centered_index =
      permutation_p_centered_index,
    B_perm = B_perm
  )
}

calculate_z_weighted_DGfold_one_sample <- function(
    Sample_name,
    result_list_object,
    annotation_tbl,
    p_cutoff = 0.05,
    B_perm = 1000,
    seed = 123
) {

  message("Processing z-score-weighted DGfold index: ", Sample_name)

  df <- result_list_object[[Sample_name]] %>%
    dplyr::mutate(
      ensembl_id = trimws(ensembl_id),
      z_score = as.numeric(z_score),
      p_value = as.numeric(p_value)
    )

  sample_tbl <- df %>%
    dplyr::inner_join(
      annotation_tbl,
      by = "ensembl_id"
    ) %>%
    dplyr::filter(
      is.finite(z_score),
      is.finite(Normalized_MFE)
    ) %>%
    dplyr::group_by(
      ensembl_id
    ) %>%
    dplyr::summarise(
      z_score = mean(z_score, na.rm = TRUE),
      p_value = min(p_value, na.rm = TRUE),
      Normalized_MFE = dplyr::first(Normalized_MFE),
      .groups = "drop"
    ) %>%
    dplyr::filter(
      is.finite(z_score),
      is.finite(Normalized_MFE)
    )

  output_tbl <- dplyr::bind_rows(
    calculate_z_weighted_DGfold_one_universe(
      sample_tbl = sample_tbl,
      Sample_name = Sample_name,
      universe_label = "All finite z-score genes",
      p_cutoff = p_cutoff,
      B_perm = B_perm,
      seed = seed
    ),
    calculate_z_weighted_DGfold_one_universe(
      sample_tbl = sample_tbl,
      Sample_name = Sample_name,
      universe_label = "Significant genes p_value < 0.05",
      p_cutoff = p_cutoff,
      B_perm = B_perm,
      seed = seed
    )
  )

  return(output_tbl)
}

iLINCS_3.33uM_z_score_weighted_DGfold_results <- dplyr::bind_rows(
  lapply(
    names(result_list_2_filtered),
    function(x) {
      calculate_z_weighted_DGfold_one_sample(
        Sample_name = x,
        result_list_object = result_list_2_filtered,
        annotation_tbl = iLINCS_DGfold_annotation,
        p_cutoff = 0.05,
        B_perm = 1000,
        seed = 123
      )
    }
  )
) %>%
  dplyr::mutate(
    dose = "3.33uM",
    cell_line = sub("_.*", "", SampleName),
    Drug = sub("^[^_]*_([^_]*).*", "\\1", SampleName),
    iLINCs_cell_line_drug = paste(
      cell_line,
      Drug,
      sep = "_"
    ),
    direction_delta_UP_minus_DOWN = dplyr::case_when(
      delta_z_weighted_Normalized_MFE_UP_minus_DOWN > 0 ~
        "Positive z-score genes have higher / less negative DGfold",
      delta_z_weighted_Normalized_MFE_UP_minus_DOWN < 0 ~
        "Positive z-score genes have lower / more negative DGfold",
      TRUE ~
        "No directional difference"
    ),
    direction_centered_index = dplyr::case_when(
      z_weighted_centered_DGfold_perturbation_index > 0 ~
        "Positive z-scores enriched for higher-than-average DGfold",
      z_weighted_centered_DGfold_perturbation_index < 0 ~
        "Positive z-scores enriched for lower-than-average DGfold",
      TRUE ~
        "No directional enrichment"
    )
  ) %>%
  dplyr::relocate(
    dose,
    cell_line,
    Drug,
    iLINCs_cell_line_drug,
    .after = SampleName
  )

iLINCS_3.33uM_z_score_weighted_DGfold_results_csv <-
  zweighted_dgfold_results_file

readr::write_csv(
  iLINCS_3.33uM_z_score_weighted_DGfold_results,
  iLINCS_3.33uM_z_score_weighted_DGfold_results_csv
)

iLINCS_3.33uM_z_score_weighted_DGfold_results_for_STEP14 <-
  iLINCS_3.33uM_z_score_weighted_DGfold_results %>%
  dplyr::filter(
    universe == "Significant genes p_value < 0.05"
  ) %>%
  dplyr::select(
    SampleName,
    z_weighted_MFE_UP,
    z_weighted_MFE_DOWN,
    delta_z_weighted_Normalized_MFE_UP_minus_DOWN,
    z_weighted_centered_DGfold_perturbation_index,
    z_weighted_standardized_DGfold_perturbation_index,
    spearman_z_score_Normalized_MFE,
    spearman_z_score_Normalized_MFE_p_value,
    permutation_p_delta_UP_minus_DOWN,
    permutation_p_centered_index,
    n_genes_used,
    n_positive_z,
    n_negative_z
  ) %>%
  dplyr::rename(
    z_weighted_DGfold_n_genes_used = n_genes_used,
    z_weighted_DGfold_n_positive_z = n_positive_z,
    z_weighted_DGfold_n_negative_z = n_negative_z
  )

print(iLINCS_3.33uM_z_score_weighted_DGfold_results)

print(iLINCS_3.33uM_z_score_weighted_DGfold_results_for_STEP14)

all_results <- read.csv(final_results_file,
                                     header = TRUE, sep = ",", stringsAsFactors = FALSE)

all_results_2 <- all_results

all_results_2$SampleName <- gsub("_MFE.pdf", "", all_results_2$SampleName)

Sig_info_iLINCs_metadata <- read.csv(signature_metadata_file,
                                     header = TRUE, sep = ",", stringsAsFactors = FALSE)

Sig_info_iLINCs_metadata <- dplyr::filter(Sig_info_iLINCs_metadata, pert_idose == "3.33 uM") %>%
  dplyr::filter(pert_itime == "24 h") %>%
  dplyr::rename(SampleName = sig_id) %>%
  dplyr::rename(Drug = cmap_name, cell_line = cell_iname)

Sig_info_iLINCs_metadata$SampleName <- gsub("/", "_", Sig_info_iLINCs_metadata$SampleName)

Sig_info_iLINCs_metadata$SampleName <- paste0(Sig_info_iLINCs_metadata$Drug, "_", Sig_info_iLINCs_metadata$SampleName)

Sig_info_iLINCs_metadata$SampleName <- paste0(Sig_info_iLINCs_metadata$cell_line, "_", Sig_info_iLINCs_metadata$SampleName)

all_results_2 <- dplyr::left_join(all_results_2, Sig_info_iLINCs_metadata, by ="SampleName") %>%
  dplyr::select(SampleName, Net_MFE_change, MeanMFE_UP, MeanMFE_DOWN, TTEST_pValue, wilcox_test_all_UP_vs_DOWN_pValue, ks_test_all_UP_vs_DOWN_pValue, ks_test_MFE_075_above_pValue, ks_test_MFE_025_below_pValue,
                Proliferation_Science, PCNA, Essential_cell_cycle_genes, Homologous_recombination, Fanconi_anemia, High_fid_pol, Drug, cell_line)

Cell_info_beta <- read.csv(cell_metadata_file,
                           header = TRUE, sep = ",", stringsAsFactors = FALSE) %>%
  dplyr::select(cell_iname, primary_disease, subtype) %>%
  dplyr::rename(cell_line = cell_iname)

all_results_2 <- dplyr::left_join(all_results_2, Cell_info_beta, by = "cell_line")

all_results_2_complete <- dplyr::mutate(all_results_2, Log10_MFE_change = log10(MeanMFE_DOWN/MeanMFE_UP)) %>%
  dplyr::relocate(Log10_MFE_change, .before = TTEST_pValue)

print(nrow(all_results_2_complete))

write.csv(all_results_2_complete, file = final_results_file, row.names = FALSE)

figure4_results <- list(
  iLINCS_3.33uM_z_score_weighted_DGfold_results =
    iLINCS_3.33uM_z_score_weighted_DGfold_results,
  all_results_2_complete = all_results_2_complete
)

saveRDS(figure4_results, figure4_results_file)

cat("\n3.33 uM iLINCS analysis complete.\n")
cat("Figure 4 input objects written to:\n", figure4_results_file, "\n")
