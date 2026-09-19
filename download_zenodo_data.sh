#!/usr/bin/env bash
#
# Download this paper's processed data from Zenodo into ./data/.
#
#   ./download_zenodo_data.sh              # download everything, skip what is current
#   ./download_zenodo_data.sh --check      # verify what is already there, download nothing
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
[ "${1:-}" = "--check" ] && CHECK_ONLY=1

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

n_total=$(wc -l < "$list" | tr -d ' ')
echo "Record $RECORD_ID -- $SOURCE"
echo "Destination: $DEST"
echo "$n_total files"
echo

md5_of() { md5 -q "$1" 2>/dev/null || md5sum "$1" 2>/dev/null | cut -d' ' -f1; }

ok=0; fetched=0; failed=0
while IFS=$'\t' read -r name size sum url; do
  target="$DEST/$name"

  if [ -f "$target" ] && [ "$(md5_of "$target")" = "$sum" ]; then
    printf "  current  %9.1f MB  %s\n" "$(echo "scale=1; $size/1048576" | bc)" "$name"
    ok=$((ok + 1))
    continue
  fi

  if [ "$CHECK_ONLY" -eq 1 ]; then
    if [ -f "$target" ]; then
      printf "  STALE    %9.1f MB  %s\n" "$(echo "scale=1; $size/1048576" | bc)" "$name"
    else
      printf "  absent   %9.1f MB  %s\n" "$(echo "scale=1; $size/1048576" | bc)" "$name"
    fi
    failed=$((failed + 1))
    continue
  fi

  printf "  fetching %9.1f MB  %s ... " "$(echo "scale=1; $size/1048576" | bc)" "$name"

  # No -C - here: resuming onto a stale or truncated file produces a corrupt
  # result that only the checksum catches. Always fetch whole.
  auth=()
  [ -n "${ZENODO_TOKEN:-}" ] && auth=(-H "Authorization: Bearer $ZENODO_TOKEN")

  rm -f "$target"
  if curl -sSL --retry 3 --retry-delay 5 "${auth[@]}" -o "$target" "$url" \
     && [ "$(md5_of "$target")" = "$sum" ]; then
    echo "ok"
    fetched=$((fetched + 1))
  else
    echo "FAILED (checksum mismatch or transfer error)"
    rm -f "$target"
    failed=$((failed + 1))
  fi
done < "$list"

echo
if [ "$CHECK_ONLY" -eq 1 ]; then
  echo "  current: $ok   missing or stale: $failed   (of $n_total)"
  [ "$failed" -gt 0 ] && exit 1
else
  echo "  already current: $ok   downloaded: $fetched   failed: $failed"
  [ "$failed" -gt 0 ] && exit 1
fi
echo "  all $n_total files present and checksum-verified in $DEST"
