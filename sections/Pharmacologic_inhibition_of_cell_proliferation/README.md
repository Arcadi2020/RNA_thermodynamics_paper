# Figure 3: pharmacologic inhibition of cell proliferation

These scripts analyze transcriptome-wide predicted RNA folding stability after pharmacologic inhibition of cancer-cell proliferation.

## Scripts

- `All_gene_sets_AC.R`: defines the gene sets used in the analysis.
- `Cross_validation_MFE.R`: performs the complete upstream analysis of the public RNA-seq datasets.
- `R_script_Fig_3.R`: generates Figure 3 from the processed Zenodo data.

## Input data

Data repository: **Zenodo**

Download it with `./download_zenodo_data.sh` from the repository root; it lands in `data/` and the scripts read it from there:

- `pharmacologic_inhibition_analysis_results_for_Fig_3.rds`

This RDS file is sufficient to rerun the final figure script. The complete upstream analysis additionally requires the public inputs specified in `Cross_validation_MFE.R`.

## Run

The script reads its input from `data/` by default, so no configuration is needed. Run from this section's directory:

```bash
Rscript R_script_Fig_3.R
```

Outputs are written to `Figure_3_outputs/`.

`R_script_Fig_3.R` uses the precomputed workspace from Zenodo and does not execute `Cross_validation_MFE.R`, but it does check that the file is present in the same directory, so run it from within this section rather than from a copy of the data alone.
