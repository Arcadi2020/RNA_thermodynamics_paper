# Figure 2: cell-cycle analysis

These scripts examine changes in predicted transcript-level RNA folding stability across the cell cycle in MCF7 and HeLa cells.

## Scripts

- `Analysis_RNAseq_PNAS_2017.R`: prepares the GSE94479 MCF7 data.
- `Analysis_cell_cycle_genes_MFE.R`: performs the cell-cycle and abundance-weighted MFE analyses.
- `R_script_Fig_2.R`: generates Figure 2 from the processed Zenodo data.

## Input data

Data repository: **Zenodo**

Download it with `./download_zenodo_data.sh` from the repository root; it lands in `data/` and the scripts read it from there:

- `cell_cycle_analysis_results_for_Fig_2.rds`

The two analysis scripts document the upstream processing of GSE94479 and GSE104736. The RDS file is sufficient to rerun the final figure script.

## Run

The script reads its input from `data/` by default, so no configuration is needed. Run from this section's directory:

```bash
Rscript R_script_Fig_2.R
```

Outputs are written to `Figure_2_outputs/`.
