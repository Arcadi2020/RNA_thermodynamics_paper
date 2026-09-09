# Figure 2: cell-cycle analysis

These scripts examine changes in predicted transcript-level RNA folding stability across the cell cycle in MCF7 and HeLa cells.

## Scripts

- `Analysis_RNAseq_PNAS_2017.R`: prepares the GSE94479 MCF7 data.
- `Analysis_cell_cycle_genes_MFE.R`: performs the cell-cycle and abundance-weighted MFE analyses.
- `R_script_Fig_2.R`: generates Figure 2 from the processed Zenodo data.

## Input data

Data repository: **Zenodo**

Place the following file in `../Zenodo_Cell_cycle_analysis/`:

- `cell_cycle_analysis_results_for_Fig_2.rds`

The two analysis scripts document the upstream processing of GSE94479 and GSE104736. The RDS file is sufficient to rerun the final figure script.

## Run

In `R_script_Fig_2.R`, set:

```r
figure2_results_file <- "../Zenodo_Cell_cycle_analysis/cell_cycle_analysis_results_for_Fig_2.rds"
output_dir <- "Figure_2_outputs"
```

Then run from `GitHub_files`:

```bash
Rscript R_script_Fig_2.R
```
