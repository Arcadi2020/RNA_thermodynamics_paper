# Figure 1: biological functions and predicted RNA structure

These scripts examine transcript-level predicted RNA folding stability and compare the predictions with sequence-shuffled controls, RNAstructure, PARIS2 and icSHAPE data.

## Scripts

- `Fig_1_R_script.R`: generates the Figure 1 analyses and plots.
- `PARIS2_data_analysis/`: contains the PARIS2 preprocessing and overlap analyses.
- `icSHAPE_data_analysis/`: contains the five-cell-line icSHAPE preprocessing and comparison analyses.

`Normalized_MFE` is the predicted minimum folding free energy divided by transcript length. More-negative values indicate greater predicted folding stability.

## Input data

Data repository: **Zenodo**

Download them with `./download_zenodo_data.sh` from the repository root; they land in `data/` and the scripts read them from there:

- `Gencode_TX_norm_MFE_HC_table.csv`
- `Compare_ViennaRNA_vs_RNAstructure_MFE.RData`
- `PARIS2_enriched_mRNA_ViennaRNA_pair_overlap_analysis_final_environment.rds`
- `transcript_summary_filtered_all_cell_lines.tsv.gz`
- `plot2_transcript_level_shift_summary.csv`
- `02_analyse_full_shuffle_vs_native_master_tables.RData`
- `02_analyse_full_shuffle_vs_native_selected_objects.rds`

These processed files are sufficient for `Fig_1_R_script.R`. The PARIS2 and icSHAPE folders document the upstream analyses and need additional inputs, described below.

### Upstream inputs for `icSHAPE_data_analysis/`

The three RNAfold outputs (`Prediction_MFE_Gencode_37C.txt.gz`, `GENCODE_filtered_transcripts.txt.gz`, `GENCODE_filtered_transcripts_the_last_one.txt.gz`) come from Zenodo with the rest, and `gencode.v47.annotation.gtf` is fetched by `./download_upstream_data.sh --annotation`.

The icSHAPE measurements themselves must be **downloaded manually** from the RASP database at <http://rasp.zhanglab.net/>. RASP assembles archives per request through its web interface, so there is no stable URL and `download_upstream_data.sh` cannot fetch them. Retrieve one archive per cell line and place them where the stage-01 scripts expect:

| File | Read by |
|---|---|
| `RASP_files_H9.tar` | `01_prepare_H9_CellResearch_icSHAPE_bed.R:7` |
| `RASP_files_HEK293T.tar` | `01_prepare_HEK293T_icSHAPE_bed.R:7` |
| `RASP_files_HeLa_icSHAPE.tar` | `01_prepare_HeLa_CellResearch_icSHAPE_bed.R:7` |
| `RASP_files_HepG2.tar` | `01_prepare_HepG2_CellResearch_icSHAPE_bed.R:7` |
| `RASP_files_K562_Cell_Research.tar` | `01_prepare_K562_CellResearch_icSHAPE_bed.R:7` |

The stages must then be run in order, because each reads what the previous one wrote: `01_prepare_*` → `02_map_*` → `03_join_*` (or `04_join_HEK293T`) → `05_HEK293T_..._followup` / `Analysis_icSHAPE_ViennaRNA_across_five_cell_lines.R`. The intermediate files named in each script's configuration block are produced by the preceding stage, not downloaded.

### Upstream inputs for `PARIS2_data_analysis/`

All three come from Zenodo via `./download_zenodo_data.sh`, alongside `gencode.v47.annotation.gtf` from `./download_upstream_data.sh --annotation`:

| File | Read by |
|---|---|
| `HEK293T_Amo_mRNA_enriched_CRSSANT_combined_128shards_dg.bedpe` | preprocessor `:8` |
| `ViennaRNA_transcript_level_pairs.tsv.gz` | preprocessor `:10`, quantification `:7` |
| `ViennaRNA_structures_recovered.tsv.gz` | quantification `:8` |

Run the two stages in order — the quantification reads `mapped_intra_transcript_duplexes.tsv.gz` from the preprocessor's output directory:

```bash
Rscript HEK293T_Amo_PARIS2_CRSSANT_BEDPE_enriched_mRNA_to_transcriptome_ViennaRNA_pairset_preprocessor.R
Rscript HEK293T_Amo_PARIS2_CRSSANT_ViennaRNA_enriched_mRNA_transcript_pair_overlap_quantification.R
```

Note the preprocessor writes a ~900 MB copy of the pairs table into its own output directory, so allow roughly 2 GB of scratch beyond the inputs.

## Run

The script reads its inputs from `data/` by default, so no configuration is needed. Run from this section's directory:

```bash
Rscript Fig_1_R_script.R
```

Outputs are written to `Figure_1_outputs/`.
