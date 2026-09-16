# Figure 3: pharmacologic inhibition of cell proliferation

These scripts analyze transcriptome-wide predicted RNA folding stability after pharmacologic inhibition of cancer-cell proliferation.

## Scripts

- `All_gene_sets_AC.R`: defines the gene sets used in the analysis.
- `Cross_validation_MFE.R`: performs the complete upstream analysis of the public RNA-seq datasets.
- `R_script_Fig_3.R`: generates Figure 3 from the processed Zenodo data.

## Input data

Data repository: **Zenodo**

Place the following file in `../Zenodo_Pharmacologic_inhibition_of_cell_proliferation/`:

- `pharmacologic_inhibition_analysis_results_for_Fig_3.rds`

This RDS file is sufficient to rerun the final figure script. The complete upstream analysis additionally requires the public inputs specified in `Cross_validation_MFE.R`.

## Run

In `R_script_Fig_3.R`, use:

```r
input_dir <- "../Zenodo_Pharmacologic_inhibition_of_cell_proliferation"
scripts_dir <- "."
figure_output_dir <- "Figure_3_outputs"
use_precomputed_workspace <- TRUE
```

Then run from `GitHub_files`:

```bash
Rscript R_script_Fig_3.R
```
