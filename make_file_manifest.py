#!/usr/bin/env python3
"""Generate file_manifest.csv and FILES.md: which Zenodo file feeds which analysis.

Zenodo presents a record as a flat list, so there is no folder structure to
convey which figure a file belongs to. This builds that mapping by reading the
scripts themselves rather than maintaining it by hand, so it cannot drift out
of date the way a hand-written table would.

    ./make_file_manifest.py                 # read the live deposit (needs ZENODO_TOKEN
                                            # while the record is an unpublished draft)
    ./make_file_manifest.py --from FILE     # read a saved deposit JSON instead

Outputs, both alongside this script:
    file_manifest.csv   one row per file: filename,size_bytes,md5,figures,scripts
    FILES.md            the same as a table, for the repo and the Zenodo description
"""
import argparse, csv, json, os, re, sys, urllib.request

RECORD_ID = 17363927
ROOT = os.path.dirname(os.path.abspath(__file__))
SELF_NAME = "file_manifest.csv"

SECTION_FIGURE = {
    "Biological_functions_vs_structure": "Figure 1",
    "Cell_cycle_analysis": "Figure 2",
    "Pharmacologic_inhibition_of_cell_proliferation": "Figure 3",
    "Analysis_iLINCs_data": "Figure 4",
    "Escherichia_coli_stationary_phase": "Figure 5",
}

# Files no script names directly. Archives are referenced by their contents, and
# the supplementary tables accompany the manuscript rather than any one figure.
MANUAL = {
    "Pharmacologic_inhibition_of_cell_proliferation_part1.zip":
        ("Figure 3", "upstream GREIN inputs, read from the unpacked archive"),
    "Pharmacologic_inhibition_of_cell_proliferation_part2.zip":
        ("Figure 3", "upstream GREIN inputs, read from the unpacked archive"),
}


def load_deposit(path=None):
    if path:
        return json.load(open(path))
    token = os.environ.get("ZENODO_TOKEN")
    url = f"https://zenodo.org/api/records/{RECORD_ID}"
    req = urllib.request.Request(url)
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            return json.load(r)
    except Exception:
        if not token:
            sys.exit("Record is not public yet; set ZENODO_TOKEN or pass --from FILE.")
        req = urllib.request.Request(
            f"https://zenodo.org/api/deposit/depositions/{RECORD_ID}",
            headers={"Authorization": f"Bearer {token}"})
        with urllib.request.urlopen(req, timeout=60) as r:
            return json.load(r)


def file_meta(dep):
    """name -> (size_bytes, md5). Both endpoints spell these differently."""
    out = {}
    for f in dep.get("files", []):
        name = f.get("key") or f.get("filename")
        size = f.get("size") or f.get("filesize")
        csum = (f.get("checksum") or "").replace("md5:", "")
        out[name] = (size, csum)
    return out


def scan_scripts(names):
    """Map each deposited filename to the (figure, script) pairs that read it."""
    uses = {n: set() for n in names}
    sections = os.path.join(ROOT, "sections")
    for dirpath, _, filenames in os.walk(sections):
        parts = os.path.relpath(dirpath, sections).split(os.sep)
        section = parts[0] if parts and parts[0] != "." else None
        figure = SECTION_FIGURE.get(section)
        if not figure:
            continue
        for fn in filenames:
            if not fn.endswith(".R"):
                continue
            text = open(os.path.join(dirpath, fn), errors="replace").read()
            for n in names:
                stem = n[:-3] if n.endswith(".gz") else n
                # match the filename as a path component in any quoted string,
                # which catches "path_to/x", "../../data/x" and file.path(dir, "x")
                if re.search(r'["/]' + re.escape(n) + r'"', text) or \
                   re.search(r'["/]' + re.escape(stem) + r'"', text):
                    uses[n].add((figure, fn))
    return uses


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--from", dest="src", help="saved deposit JSON")
    args = ap.parse_args()

    meta = file_meta(load_deposit(args.src))
    # The manifest is itself deposited, so its checksum changes every time this
    # runs. Listing it would publish a hash that is stale on arrival, so skip it.
    meta.pop(SELF_NAME, None)
    sizes = {k: v[0] for k, v in meta.items()}
    uses = scan_scripts(set(sizes))

    rows = []
    for name in sorted(sizes, key=str.lower):
        pairs = sorted(uses.get(name, ()))
        figures = sorted({f for f, _ in pairs})
        scripts = sorted({s for _, s in pairs})
        note = ""
        if not figures and name in MANUAL:
            fig, note = MANUAL[name]
            figures = [fig]
        elif not figures and name.lower().startswith("supplementary"):
            figures = ["--"]
            note = "manuscript supplementary table, not read by any script"
        rows.append({
            "filename": name,
            "size_bytes": sizes[name],
            "md5": meta[name][1],
            "figures": "; ".join(figures),
            "scripts": "; ".join(scripts),
            "note": note,
        })

    csv_path = os.path.join(ROOT, "file_manifest.csv")
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    # grouped markdown
    by_fig = {}
    for r in rows:
        by_fig.setdefault(r["figures"] or "--", []).append(r)

    md = ["# Zenodo file manifest",
          "",
          f"Record [10.5281/zenodo.{RECORD_ID}](https://doi.org/10.5281/zenodo.{RECORD_ID}) "
          f"-- {len(rows)} files, {sum(sizes.values())/1048576:,.0f} MB.",
          "",
          "Zenodo presents a record as a flat file list, so this table is how you tell "
          "which file belongs to which analysis. It is generated from the scripts by "
          "`make_file_manifest.py`, not maintained by hand.",
          "",
          "Fetch everything with `./download_zenodo_data.sh`; files land in `data/` "
          "and every one is verified against its checksum.",
          ""]

    def sortkey(k):
        return (k == "--", k)

    for fig in sorted(by_fig, key=sortkey):
        heading = "Not figure-specific" if fig == "--" else fig
        if ";" in fig:
            heading = f"Shared: {fig}"
        md += [f"## {heading}", "",
               "| File | Size | Read by |", "|---|---:|---|"]
        for r in by_fig[fig]:
            read_by = r["scripts"] or r["note"] or "--"
            mb = int(r["size_bytes"]) / 1048576
            md.append(f"| `{r['filename']}` | {mb:,.2f} MB | {read_by} |")
        md.append("")

    open(os.path.join(ROOT, "FILES.md"), "w").write("\n".join(md))
    print(f"  wrote file_manifest.csv and FILES.md ({len(rows)} files)")
    for fig in sorted(by_fig, key=sortkey):
        print(f"    {fig or '--':38} {len(by_fig[fig])} files")


if __name__ == "__main__":
    main()
