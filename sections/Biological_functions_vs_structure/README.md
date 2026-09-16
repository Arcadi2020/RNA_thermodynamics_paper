# Figure 1: biological functions and predicted RNA structure

These scripts examine transcript-level predicted RNA folding stability and compare the predictions with sequence-shuffled controls, RNAstructure, PARIS2 and icSHAPE data.

## Scripts

- `Fig_1_R_script.R`: generates the Figure 1 analyses and plots.
- `PARIS2_data_analysis/`: contains the PARIS2 preprocessing and overlap analyses.
- `icSHAPE_data_analysis/`: contains the five-cell-line icSHAPE preprocessing and comparison analyses.

`Normalized_MFE` is the predicted minimum folding free energy divided by transcript length. More-negative values indicate greater predicted folding stability.

## Input data

Data repository: **Zenodo**

Place these files in `../Zenodo_Biological_functions_vs_structure/`:

- `Gencode_TX_norm_MFE_HC_table.csv`
- `Compare_ViennaRNA_vs_RNAstructure_MFE.RData`
- `PARIS2_enriched_mRNA_ViennaRNA_pair_overlap_analysis_final_environment.rds`
- `transcript_summary_filtered_all_cell_lines.tsv.gz`
- `plot2_transcript_level_shift_summary.csv`
- `02_analyse_full_shuffle_vs_native_master_tables.RData`
- `02_analyse_full_shuffle_vs_native_selected_objects.rds`

These processed files are sufficient for `Fig_1_R_script.R`. The PARIS2 and icSHAPE folders document the upstream analyses and require the original public inputs named in their configuration sections.

## Run

Confirm that the input variables in `Fig_1_R_script.R` point to the Zenodo directory above, set `output_dir`, and run from `GitHub_files`:

```bash
Rscript Fig_1_R_script.R
```
