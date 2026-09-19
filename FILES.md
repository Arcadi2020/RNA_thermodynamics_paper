# Zenodo file manifest

Record [10.5281/zenodo.17363927](https://doi.org/10.5281/zenodo.17363927) -- 36 files, 1,545 MB.

Zenodo presents a record as a flat file list, so this table is how you tell which file belongs to which analysis. It is generated from the scripts by `make_file_manifest.py`, not maintained by hand.

Fetch everything with `./download_zenodo_data.sh`; files land in `data/` and every one is verified against its checksum.

## Figure 1

| File | Size | Read by |
|---|---:|---|
| `02_analyse_full_shuffle_vs_native_master_tables.RData` | 6.56 MB | Fig_1_R_script.R |
| `02_analyse_full_shuffle_vs_native_selected_objects.rds` | 2.09 MB | Fig_1_R_script.R |
| `Compare_ViennaRNA_vs_RNAstructure_MFE.RData` | 13.82 MB | Fig_1_R_script.R |
| `GENCODE_filtered_transcripts.txt.gz` | 7.63 MB | 03_prepare_ViennaRNA_position_table_for_icSHAPE.R |
| `GENCODE_filtered_transcripts_the_last_one.txt.gz` | 7.88 MB | 03_prepare_ViennaRNA_position_table_for_icSHAPE.R |
| `HEK293T_Amo_mRNA_enriched_CRSSANT_combined_128shards_dg.bedpe` | 2.68 MB | HEK293T_Amo_PARIS2_CRSSANT_BEDPE_enriched_mRNA_to_transcriptome_ViennaRNA_pairset_preprocessor.R |
| `PARIS2_enriched_mRNA_ViennaRNA_pair_overlap_analysis_final_environment.rds` | 0.85 MB | Fig_1_R_script.R |
| `plot2_transcript_level_shift_summary.csv` | 0.00 MB | Analysis_icSHAPE_ViennaRNA_across_five_cell_lines.R; Fig_1_R_script.R |
| `Prediction_MFE_Gencode_37C.txt.gz` | 103.92 MB | 03_prepare_ViennaRNA_position_table_for_icSHAPE.R |
| `transcript_summary_filtered_all_cell_lines.tsv.gz` | 15.60 MB | Analysis_icSHAPE_ViennaRNA_across_five_cell_lines.R; Fig_1_R_script.R |
| `ViennaRNA_structures_recovered.tsv.gz` | 172.78 MB | HEK293T_Amo_PARIS2_CRSSANT_ViennaRNA_enriched_mRNA_transcript_pair_overlap_quantification.R |
| `ViennaRNA_transcript_level_pairs.tsv.gz` | 883.40 MB | HEK293T_Amo_PARIS2_CRSSANT_BEDPE_enriched_mRNA_to_transcriptome_ViennaRNA_pairset_preprocessor.R; HEK293T_Amo_PARIS2_CRSSANT_ViennaRNA_enriched_mRNA_transcript_pair_overlap_quantification.R |

## Shared: Figure 1; Figure 2; Figure 3

| File | Size | Read by |
|---|---:|---|
| `Gencode_TX_norm_MFE_HC_table.csv` | 26.37 MB | Analysis_cell_cycle_genes_MFE.R; Cross_validation_MFE.R; Fig_1_R_script.R |

## Figure 2

| File | Size | Read by |
|---|---:|---|
| `cell_cycle_analysis_results_for_Fig_2.rds` | 31.96 MB | Analysis_cell_cycle_genes_MFE.R; R_script_Fig_2.R |
| `GSE104736_GeneLevel_Raw_data.csv` | 3.56 MB | Analysis_cell_cycle_genes_MFE.R |
| `GSE104736_TranscriptLevel_Raw_data.csv` | 35.93 MB | Analysis_cell_cycle_genes_MFE.R |
| `GSE104736_TranscriptLevel_TPM_data.csv` | 21.29 MB | Analysis_cell_cycle_genes_MFE.R |
| `GSE94479_TranscriptLevel_Raw_data.csv` | 6.55 MB | Analysis_cell_cycle_genes_MFE.R |
| `GSE94479_TranscriptLevel_TPM_data.csv` | 7.09 MB | Analysis_cell_cycle_genes_MFE.R |
| `refFlat_hg19_gtf.gz` | 9.20 MB | Analysis_RNAseq_PNAS_2017.R |

## Figure 3

| File | Size | Read by |
|---|---:|---|
| `pharmacologic_inhibition_analysis_results_for_Fig_3.rds` | 57.34 MB | R_script_Fig_3.R |
| `Pharmacologic_inhibition_of_cell_proliferation_part1.zip` | 38.56 MB | upstream GREIN inputs, read from the unpacked archive |
| `Pharmacologic_inhibition_of_cell_proliferation_part2.zip` | 5.67 MB | upstream GREIN inputs, read from the unpacked archive |

## Figure 4

| File | Size | Read by |
|---|---:|---|
| `iLINCS_analysis_results_for_Fig_4.rds` | 0.07 MB | R_script_Fig_4.R |
| `Table_ANNOTATED_iLINCs.csv` | 2.75 MB | Analysis_MFE_iLINCS_all_datasets_10uM_cancer_lines.R; Analysis_MFE_iLINCS_all_datasets_3.33uM_cancer_lines.R |

## Figure 5

| File | Size | Read by |
|---|---:|---|
| `GSE224463_analysis_results_for_Fig_5.rds` | 2.24 MB | Analysis_GSE224463_stationary_phase.R; R_script_Fig_5.R |
| `GSE226643_analysis_results_for_Fig_5.rds` | 0.42 MB | Analysis_GSE226643_stationary_phase.R; R_script_Fig_5.R |
| `MG1655_gene_info.csv` | 0.51 MB | Analysis_GSE224463_stationary_phase.R; Analysis_GSE226643_stationary_phase.R |
| `Prediction_Ecoli_MFE_HC_df.csv` | 2.07 MB | Analysis_GSE224463_stationary_phase.R; Analysis_GSE226643_stationary_phase.R |

## Not figure-specific

| File | Size | Read by |
|---|---:|---|
| `Supplementary_table_1.csv` | 58.64 MB | manuscript supplementary table, not read by any script |
| `Supplementary_table_2.csv` | 15.84 MB | manuscript supplementary table, not read by any script |
| `Supplementary_Table_3.xlsx` | 0.03 MB | manuscript supplementary table, not read by any script |
| `Supplementary_Table_4.csv` | 0.01 MB | manuscript supplementary table, not read by any script |
| `Supplementary_table_5.csv` | 0.00 MB | manuscript supplementary table, not read by any script |
| `Supplementary_table_6.csv` | 0.11 MB | manuscript supplementary table, not read by any script |
| `Supplementary_Table_7.csv` | 1.23 MB | manuscript supplementary table, not read by any script |
