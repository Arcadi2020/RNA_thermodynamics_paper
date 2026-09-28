#!/usr/bin/env bash
#
# Download this paper's processed data from Zenodo into ./data/.
#
#   ./download_zenodo_data.sh              # download everything, skip what is current
#   ./download_zenodo_data.sh --check      # verify what is already there, download nothing
#   ./download_zenodo_data.sh --figure 3   # only what Figure 3 needs
#   ./download_zenodo_data.sh --jobs 8     # more parallel transfers (default 4)
#
# The largest single file is ViennaRNA_transcript_level_pairs.tsv.gz (883 MB,
# Figure 1). Figure 3's GREIN counts are one zip per GEO accession (GSE*.zip,
# ~713 MB together) plus GREIN_supporting_inputs.zip (the per-series metadata).
# If you only want one figure, --figure reads file_manifest.csv and fetches
# just that figure's files plus the ones no figure claims.
#
# Figure 3's GREIN archives land here as zips; expand them into the directory
# Cross_validation_MFE.R reads with:
#   Rscript sections/Pharmacologic_inhibition_of_cell_proliferation/unpack_GREIN_counts.R \
#           ./data <input_dir>
#
# On a healthy link --jobs makes no measurable difference: 6 files / 90 MB took
# 5s serial and 6s at --jobs 6, because a single connection already saturates
# the link at ~18 MB/s. It is worth raising only when per-connection throughput
# is the bottleneck -- Zenodo's download path was measured at 462 KB/s during a
# degraded spell, and that is the case parallelism helps with.
#
# Every file is checked against the MD5 recorded by Zenodo. Files that are
# already present and match are skipped, so re-running is cheap and a partial
# download is repaired rather than left truncated.
#
# If an upload appears stuck as "pending (checksum not yet calculated)":
# Zenodo runs two APIs. This script reads the legacy deposit API, but the web
# UI uses the InvenioRDM one, and a file whose transfer never completed is left
# there as a zero-byte entry with status "pending". That entry holds the key, so
# every later attempt to upload the same filename conflicts with it and fails --
# including attempts that look like server throttling. The legacy API does not
# show pending entries at all, so the deposit looks clean while uploads keep
# failing. Check for and clear them before retrying:
#
#   curl -H "Authorization: Bearer $ZENODO_TOKEN" \
#     https://zenodo.org/api/records/17363927/draft/files \
#     | python3 -c 'import json,sys; [print(e["key"], e["status"], e["size"]) \
#         for e in json.load(sys.stdin)["entries"] if e["status"] != "completed"]'
#
#   curl -X DELETE -H "Authorization: Bearer $ZENODO_TOKEN" \
#     https://zenodo.org/api/records/17363927/draft/files/<filename>
#
# Clearing one such entry took an 883 MB upload from ten consecutive failures
# (best 66%, degrading to 196 B/s) to succeeding first time at 1.7 MB/s.
#
# The record is currently an UNPUBLISHED DRAFT, which the public API cannot
# see. Until it is published you need a Zenodo token with deposit:write from
# an account that owns the record:
#
#   export ZENODO_TOKEN=...
#
# Once the record is published the token is no longer needed and the public
# endpoint is used automatically.

set -euo pipefail

RECORD_ID=17363927
API=https://zenodo.org/api
DEST="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/data"

CHECK_ONLY=0
JOBS=4
FIGURE=""
while [ $# -gt 0 ]; do
  case "$1" in
    --check)  CHECK_ONLY=1 ;;
    --jobs)   JOBS="${2:?--jobs needs a number}"; shift ;;
    --figure) FIGURE="${2:?--figure needs a number, e.g. 3}"; shift ;;
    *) echo "Unknown option: $1"; exit 2 ;;
  esac
  shift
done

MANIFEST="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/file_manifest.csv"
if [ -n "$FIGURE" ] && [ ! -f "$MANIFEST" ]; then
  echo "--figure needs file_manifest.csv beside this script."
  echo "Generate it with ./make_file_manifest.py, or omit --figure."
  exit 2
fi

mkdir -p "$DEST"

# ---- locate the record -------------------------------------------------------
# Published records answer on /records/<id>; an unpublished draft only answers on
# the deposit endpoint, and only with a token.
meta=$(mktemp)
trap 'rm -f "$meta"' EXIT

code=$(curl -sS -o "$meta" -w '%{http_code}' "$API/records/$RECORD_ID")

if [ "$code" != "200" ]; then
  if [ -z "${ZENODO_TOKEN:-}" ]; then
    echo "Record $RECORD_ID is not publicly available (HTTP $code)."
    echo
    echo "It is still an unpublished draft. Set a Zenodo token from an account"
    echo "that owns the record and re-run:"
    echo
    echo "    export ZENODO_TOKEN=..."
    echo "    $0"
    exit 1
  fi
  code=$(curl -sS -H "Authorization: Bearer $ZENODO_TOKEN" \
              -o "$meta" -w '%{http_code}' \
              "$API/deposit/depositions/$RECORD_ID")
  if [ "$code" != "200" ]; then
    echo "Could not read record $RECORD_ID even with a token (HTTP $code)."
    echo "401 => bad or expired token. 403 => the token's account does not own it."
    exit 1
  fi
  SOURCE="draft (authenticated)"
else
  SOURCE="published record"
fi

# ---- build the file list -----------------------------------------------------
# The two endpoints name their fields differently: published records use
# key/size/checksum/links.self, drafts use filename/filesize/checksum/links.download.
list=$(mktemp)
trap 'rm -f "$meta" "$list"' EXIT

python3 - "$meta" > "$list" <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
for f in d.get("files", []):
    name = f.get("key") or f.get("filename")
    size = f.get("size") or f.get("filesize")
    csum = (f.get("checksum") or "").replace("md5:", "")
    links = f.get("links", {})
    url = links.get("download") or links.get("content") or links.get("self")
    print("\t".join([name, str(size), csum, url]))
PY

if [ -n "$FIGURE" ]; then
  # Keep the rows whose figure column mentions this figure, and the rows no
  # figure claims -- the manifest, the supplementary tables and anything
  # shared, which are small and wanted either way.
  kept=$(mktemp); trap 'rm -f "$meta" "$list" "$kept"' EXIT
  python3 - "$MANIFEST" "$FIGURE" "$list" > "$kept" <<'PYFILTER'
import csv, sys
manifest, figure, listing = sys.argv[1], sys.argv[2], sys.argv[3]
wanted, unclaimed = set(), set()
for row in csv.DictReader(open(manifest)):
    figures = row["figures"]
    if f"Figure {figure}" in figures:
        wanted.add(row["filename"])
    elif not figures or figures == "--":
        unclaimed.add(row["filename"])
keep = wanted | unclaimed
for line in open(listing):
    if line.split("\t")[0] in keep or line.split("\t")[0] == "file_manifest.csv":
        sys.stdout.write(line)
PYFILTER
  mv "$kept" "$list"
fi

n_total=$(wc -l < "$list" | tr -d ' ')
echo "Record $RECORD_ID -- $SOURCE"
[ -n "$FIGURE" ] && echo "Restricted to Figure $FIGURE (per file_manifest.csv)"
echo "Destination: $DEST"
echo "$n_total files"
echo

md5_of() { md5 -q "$1" 2>/dev/null || md5sum "$1" 2>/dev/null | cut -d' ' -f1; }

# Pass 1: decide what actually needs fetching, so the parallel pass is pure I/O.
todo=$(mktemp); trap 'rm -f "$meta" "$list" "$todo"' EXIT
ok=0; failed=0
while IFS=$'\t' read -r name size sum url; do
  target="$DEST/$name"
  if [ -f "$target" ] && [ "$(md5_of "$target")" = "$sum" ]; then
    printf "  current  %9.1f MB  %s\n" "$(echo "scale=1; $size/1048576" | bc)" "$name"
    ok=$((ok + 1))
    continue
  fi
  if [ "$CHECK_ONLY" -eq 1 ]; then
    if [ -f "$target" ]; then st=STALE; else st=absent; fi
    printf "  %-8s %9.1f MB  %s\n" "$st" "$(echo "scale=1; $size/1048576" | bc)" "$name"
    failed=$((failed + 1))
    continue
  fi
  printf '%s\t%s\t%s\t%s\n' "$name" "$size" "$sum" "$url" >> "$todo"
done < "$list"

fetched=0
n_todo=$(wc -l < "$todo" | tr -d ' ')

if [ "$CHECK_ONLY" -eq 0 ] && [ "$n_todo" -gt 0 ]; then
  echo
  echo "  fetching $n_todo file(s), $JOBS at a time"

  # Pass 2: transfer in parallel, in batches of $JOBS. Backgrounded subshells
  # with a plain `wait` keep this portable -- `wait -n` needs bash 4.3 and
  # macOS still ships 3.2.
  #
  # No -C -: resuming onto a stale or truncated file yields a corrupt result
  # that only the checksum would catch.
  fetch_one() {
    local name="$1" url="$2" target="$DEST/$1"
    rm -f "$target"
    if [ -n "${ZENODO_TOKEN:-}" ]; then
      curl -sSL --retry 3 --retry-delay 5 -H "Authorization: Bearer ${ZENODO_TOKEN}" -o "$target" "$url" 2>/dev/null
    else
      curl -sSL --retry 3 --retry-delay 5 -o "$target" "$url" 2>/dev/null
    fi
  }

  batch=0
  while IFS=$'\t' read -r name size sum url; do
    fetch_one "$name" "$url" &
    batch=$((batch + 1))
    if [ "$batch" -ge "$JOBS" ]; then wait; batch=0; fi
  done < "$todo"
  wait

  # Pass 3: verify every transfer against Zenodo's checksum.
  while IFS=$'\t' read -r name size sum url; do
    target="$DEST/$name"
    if [ -f "$target" ] && [ "$(md5_of "$target")" = "$sum" ]; then
      printf "  ok       %9.1f MB  %s\n" "$(echo "scale=1; $size/1048576" | bc)" "$name"
      fetched=$((fetched + 1))
    else
      printf "  FAILED   %9.1f MB  %s (checksum mismatch or transfer error)\n" "$(echo "scale=1; $size/1048576" | bc)" "$name"
      rm -f "$target"
      failed=$((failed + 1))
    fi
  done < "$todo"
fi

echo
if [ "$CHECK_ONLY" -eq 1 ]; then
  echo "  current: $ok   missing or stale: $failed   (of $n_total)"
  [ "$failed" -gt 0 ] && exit 1
else
  echo "  already current: $ok   downloaded: $fetched   failed: $failed"
  [ "$failed" -gt 0 ] && exit 1
fi
echo "  all $n_total files present and checksum-verified in $DEST"
