# Figure 5: *Escherichia coli* stationary-phase analysis

These scripts analyze predicted transcriptome-wide RNA folding stability during entry into stationary phase in GSE224463 and during nitrogen or sulfur limitation in GSE226643.

## Scripts

- `Analysis_GSE224463_stationary_phase.R`: performs the GSE224463 time-course analysis.
- `Analysis_GSE226643_stationary_phase.R`: performs the GSE226643 nutrient-limitation analysis.
- `R_script_Fig_5.R`: generates Figure 5 from the processed Zenodo data.

## Input data

Data repository: **Zenodo**

Download them with `./download_zenodo_data.sh` from the repository root; they land in `data/` and the scripts read them from there:

- `GSE224463_analysis_results_for_Fig_5.rds`
- `GSE226643_analysis_results_for_Fig_5.rds`

These two RDS files are sufficient to rerun the final figure script. The dataset-specific scripts document the complete upstream analyses and identify their original public inputs.

## Run

The script reads its inputs from `data/` by default. Run from this section's directory:

```bash
Rscript R_script_Fig_5.R
```

Outputs are written to `Figure_5_outputs/`. Explicit paths can still be passed instead, and all three are required together:

```bash
Rscript R_script_Fig_5.R <GSE224463.rds> <GSE226643.rds> <output_dir>
```
