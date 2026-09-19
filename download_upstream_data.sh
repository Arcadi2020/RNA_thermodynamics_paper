#!/usr/bin/env bash
#
# Download the third-party upstream inputs used by the analysis scripts into
# ./data/upstream/. These are the public datasets the paper builds on -- GEO
# series, GENCODE/UCSC annotation, and the LINCS 2020 beta release. They are
# deliberately NOT part of the Zenodo deposit, which carries only processed
# results; use download_zenodo_data.sh for those.
#
#   ./download_upstream_data.sh --list        what each group covers, download nothing
#   ./download_upstream_data.sh               annotation + GEO + LINCS metadata (~1.2 GB)
#   ./download_upstream_data.sh --annotation  GENCODE and UCSC only (~63 MB)
#   ./download_upstream_data.sh --geo         GEO series only (~1 GB)
#   ./download_upstream_data.sh --lincs       LINCS metadata only (~1.1 GB)
#   ./download_upstream_data.sh --lincs-gctx  ALSO the 111 GB level 3 matrix
#   ./download_upstream_data.sh --check       report what is present
#
# IMPORTANT -- this script does not fetch everything the upstream scripts want,
# but nothing is unobtainable any more. Run --list for the full picture:
#
#   fetchable here      the GEO, GENCODE and LINCS files below
#   on Zenodo           the authors' own ViennaRNA outputs and the GREIN count
#                       tables -- use download_zenodo_data.sh
#   in the repo already All_gene_sets_AC.R
#   regenerable         ~18 intermediates that earlier scripts in the same
#                       pipeline produce (run the icSHAPE/PARIS2 stages in order)
#   manual download     the five RASP icSHAPE archives, which have no stable URL
#
# Nothing here is verified against a checksum, because these sources do not
# publish one. Sizes are checked against the server's Content-Length instead.
#
# The sources themselves live in upstream_sources.csv next to this script,
# one row per file: group,subdir,filename,url,post. Add or change a source
# by editing that file -- no shell changes needed. 'post' is no, yes
# (gunzip after download) or untar.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="$ROOT/data/upstream"

MANIFEST_CSV="$ROOT/upstream_sources.csv"

if [ ! -f "$MANIFEST_CSV" ]; then
  echo "Source list not found: $MANIFEST_CSV"
  echo "It ships alongside this script; fetch both from the repository."
  exit 1
fi

# group,subdir,filename,url,post  -- header row skipped.
# Edit upstream_sources.csv to add or change a source; no shell changes needed.
MANIFEST=$(tail -n +2 "$MANIFEST_CSV" | tr ',' '|')

show_list() {
  cat <<'EOF'
GROUPS THIS SCRIPT CAN DOWNLOAD
===============================

  --annotation   GENCODE v47 GTF + tRNA GTF                              ~56 MB
                 needed by: icSHAPE 02_map_*, PARIS2 preprocessor

  --geo          GSE94479  three GSM bigwigs + four DESeq/GRO tables
                 GSE104736 ReadCounts + RPKM
                 GSE224463 RAW.tar   GSE226643 X.csv                    ~1 GB
                 needed by: Analysis_RNAseq_PNAS_2017.R,
                            Analysis_GSE224463/GSE226643_stationary_phase.R

  --lincs        LINCS 2020 beta metadata: cellinfo, geneinfo,
                 siginfo (444 MB), instinfo (644 MB)                    ~1.1 GB
                 needed by: lincs_maayanlab.py, and by
                 Analysis_MFE_iLINCS_*.R which reads cellinfo_beta and
                 siginfo_beta_LINCS_metadata. These have stable CLUE URLs,
                 so they are fetched here rather than deposited.

  --lincs-gctx   LINCS 2020 level 3 expression matrix                   ~111 GB
                 needed by: lincs_maayanlab.py. Opt-in: it is enormous.


WHAT THIS SCRIPT CANNOT GET
===========================

  Already in the repository
    All_gene_sets_AC.R
      -> sections/Pharmacologic_inhibition_of_cell_proliferation/

  NOT UCSC's refFlat, despite the name. Analysis_RNAseq_PNAS_2017.R reads
  refFlat_hg19.gz with readGFFAsGRanges(), so it needs the UCSC refFlat table
  converted to GTF (9 columns, source "hg19_refFlat", gene_id/transcript_id
  attributes) -- not the raw 11-column genePred that UCSC serves as
  refFlat.txt.gz. This script previously downloaded the raw table under that
  name, which fails with "line 1 has more than 9 tab-separated columns". The
  GTF-format file is deposited on Zenodo as refFlat_hg19_gtf.gz instead.

  Intermediates produced by earlier scripts in the same pipeline.
  Run the stages in order rather than downloading these.
    *_icSHAPE_bed_clean_sorted.tsv.gz              <- icSHAPE 01_prepare_*
    *_icSHAPE_positions_mapped_to_transcripts.tsv.gz <- icSHAPE 02_map_*
    *_ViennaRNA_joined_per_position_table.tsv.gz   <- icSHAPE 03/04_join_*
    ViennaRNA_per_position_state_table_for_icSHAPE.tsv.gz <- icSHAPE 03_prepare_*
    mapped_intra_transcript_duplexes.tsv.gz        <- PARIS2 preprocessor
    Analysis_RNAseq_PNAS_2017_environment.RData    <- Analysis_RNAseq_PNAS_2017.R

  GREIN reprocessings of public GEO series. The underlying GEO data is public,
  but these particular files are not: GEO does not distribute them in this form
  (it ships GSE94479_Deseq_*/GROCountNorm and GSE104736_RPKM/ReadCounts), and
  GREIN is an interactive Shiny app with no REST API or stable download URL, so
  they cannot be fetched by script. A fresh GREIN export is also not guaranteed
  to match byte-for-byte, since its pipeline and annotation versions move.

  For that reason they belong in the Zenodo deposit rather than here, and most
  of them are already there:
    Pharmacologic_..._part1.zip  carries the GREIN_raw_counts and
                                 GREIN_series_matrices_data trees (46 entries)
    deposited individually       GSE94479_TranscriptLevel_{Raw,TPM}_data.csv
                                 GSE104736_TranscriptLevel_{Raw,TPM}_data.csv
                                 GSE104736_GeneLevel_Raw_data.csv

  If you do need to regenerate one, export it per accession from
  https://www.ilincs.org/apps/grein/ and expect minor differences.

  The two Escherichia coli scripts expect a project directory per accession:

      GSE224463/  metadata/  GEO_download/  Zenodo_8083605/  SRA_selected_samples/
      GSE226643/  metadata/  GEO_download/  author_GitHub_repository/  SRA_selected_samples/

  This script fills GEO_download/ only, and that is all it needs to do. The
  scripts are otherwise self-sufficient: they fetch Zenodo_8083605/ and
  author_GitHub_repository/ at run time (Zenodo API, and the GitHub zipball
  endpoint for hyungyulim/Ecoli_transition_stationary), write their own
  metadata/*.tsv manifests, and resolve SRA_selected_samples/FASTQ/ from ENA
  themselves. The only external inputs they still lack are
  Prediction_Ecoli_MFE.csv and MG1655_gene_info.csv, listed below.

  Assembled through the RASP web interface; no stable direct URL.
  Fetch from http://rasp.zhanglab.net/ per cell line.
    RASP_files_H9.tar  RASP_files_HEK293T.tar  RASP_files_HeLa_icSHAPE.tar
    RASP_files_HepG2.tar  RASP_files_K562_Cell_Research.tar

  DEPOSITED ON ZENODO -- no longer a gap. These are the authors' own upstream
  outputs. They were absent from both the deposit and every public archive
  until the authors supplied them, and they are now in the deposit, so
  download_zenodo_data.sh fetches them into data/ with the rest. Several are
  deposited under a different name from the one the scripts originally
  declared; the scripts have been updated to match.

    declared originally                         deposited as
    combined_CRSSANT_DG.bedpe                   HEK293T_Amo_mRNA_enriched_CRSSANT_
                                                combined_128shards_dg.bedpe
                                                (the original was a placeholder name)
    Prediction_Ecoli_MFE.csv                    Prediction_Ecoli_MFE_HC_df.csv
    ViennaRNA_structures_recovered.tsv          ViennaRNA_structures_recovered.tsv.gz
    Prediction_MFE_Gencode_37C.txt              Prediction_MFE_Gencode_37C.txt.gz
    GENCODE_filtered_transcripts.txt            GENCODE_filtered_transcripts.txt.gz
    GENCODE_filtered_transcripts_the_last_one   ...the_last_one.txt.gz
    ViennaRNA_transcript_level_pairs.tsv.gz     (same)
    MG1655_gene_info.csv                        (same)
    Table_ANNOTATED_iLINCs.csv                  (same)

  The .gz forms are read directly: read.delim() and data.table::fread() both
  decompress transparently, and the parsed result is identical.

  Verified end to end after these were deposited: the PARIS2 preprocessor and
  the pair-overlap quantification both run to completion from data/ alone.

EOF
}

WANT_ANNOTATION=0; WANT_GEO=0; WANT_LINCS=0; WANT_GCTX=0; CHECK_ONLY=0

if [ $# -eq 0 ]; then
  WANT_ANNOTATION=1; WANT_GEO=1; WANT_LINCS=1
else
  for a in "$@"; do
    case "$a" in
      --list)       show_list; exit 0 ;;
      --annotation) WANT_ANNOTATION=1 ;;
      --geo)        WANT_GEO=1 ;;
      --lincs)      WANT_LINCS=1 ;;
      --lincs-gctx) WANT_LINCS=1; WANT_GCTX=1 ;;
      --check)      CHECK_ONLY=1; WANT_ANNOTATION=1; WANT_GEO=1; WANT_LINCS=1 ;;
      *) echo "Unknown option: $a"; echo "Try --list"; exit 2 ;;
    esac
  done
fi

wanted() {
  case "$1" in
    annotation) [ "$WANT_ANNOTATION" -eq 1 ] ;;
    geo)        [ "$WANT_GEO" -eq 1 ] ;;
    lincs)      [ "$WANT_LINCS" -eq 1 ] ;;
    lincs-gctx) [ "$WANT_GCTX" -eq 1 ] ;;
    *) false ;;
  esac
}

echo "Destination: $DEST"
echo

present=0; fetched=0; failed=0; skipped=0

while IFS='|' read -r group sub name url gz; do
  [ -z "${group:-}" ] && continue
  if ! wanted "$group"; then skipped=$((skipped + 1)); continue; fi

  full_url="$url"
  outdir="$DEST/$sub"
  target="$outdir/$name"
  # a file that gets gunzipped is "present" under its decompressed name
  final="$target"; [ "$gz" = "yes" ] && final="${target%.gz}"
  if [ "$gz" = "untar" ] && [ -d "$outdir" ] && [ -n "$(ls -A "$outdir" 2>/dev/null)" ]; then
    printf "  present   %-52s %s\n" "$name (extracted)" "$(du -sh "$outdir" | cut -f1)"
    present=$((present + 1)); continue
  fi

  if [ -s "$final" ]; then
    printf "  present   %-52s %s\n" "$name" "$(du -h "$final" | cut -f1)"
    present=$((present + 1)); continue
  fi

  if [ "$CHECK_ONLY" -eq 1 ]; then
    printf "  ABSENT    %-52s %s\n" "$name" "$group"
    failed=$((failed + 1)); continue
  fi

  mkdir -p "$outdir"
  printf "  fetching  %-52s ... " "$name"
  # No -C -: a resumed transfer onto a stale file cannot be validated here,
  # since these sources publish no checksum.
  rm -f "$target"
  if curl -sSL --retry 3 --retry-delay 5 --fail -o "$target" "$full_url"; then
    if [ "$gz" = "yes" ]; then gunzip -f "$target"; fi
    if [ "$gz" = "untar" ]; then tar -xf "$target" -C "$outdir" && rm -f "$target"; fi
    echo "ok"
    fetched=$((fetched + 1))
  else
    echo "FAILED"
    rm -f "$target"
    failed=$((failed + 1))
  fi
done <<< "$MANIFEST"

echo
echo "  present: $present   downloaded: $fetched   failed: $failed   not requested: $skipped"
if [ "$WANT_GCTX" -eq 0 ] && [ "$CHECK_ONLY" -eq 0 ]; then
  echo
  echo "  The 111 GB LINCS level 3 matrix was not downloaded."
  echo "  Add --lincs-gctx if you need to rerun lincs_maayanlab.py."
fi
echo "  Run --list for what this script cannot fetch and where to get it."
[ "$failed" -gt 0 ] && exit 1 || exit 0
