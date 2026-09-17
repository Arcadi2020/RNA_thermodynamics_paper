# Figure 4: iLINCS data analysis

These scripts analyze drug-induced transcriptional responses and predicted RNA folding stability in cancer cell lines at 3.33 µM and 10 µM.

## Scripts

- `Analysis_MFE_iLINCS_all_datasets_3.33uM_cancer_lines.R`: performs the 3.33 µM analysis.
- `Analysis_MFE_iLINCS_all_datasets_10uM_cancer_lines.R`: performs the 10 µM analysis.
- `R_script_Fig_4.R`: generates the combined Figure 4 outputs.
- `Shell_scripts_to_download_iLINCs_data/` and `lincs_maayanlab.py`: optional utilities for downloading the original iLINCS inputs.

## Input data

Data repository: **Zenodo**

Download it with `./download_zenodo_data.sh` from the repository root; it lands in `data/` and the scripts read it from there:

- `iLINCS_analysis_results_for_Fig_4.rds`

This RDS file is sufficient to rerun the final figure script. The dose-specific analysis scripts require the original iLINCS inputs specified in their configuration sections.

## Run

The script reads its input from `data/` by default. Run from this section's directory:

```bash
Rscript R_script_Fig_4.R
```

Outputs are written to `Figure_4_outputs/`. Explicit paths can still be passed instead:

```bash
Rscript R_script_Fig_4.R <results.rds> <output_dir>
```
