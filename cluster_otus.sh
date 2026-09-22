#!/usr/bin/env bash
# =============================================================
# Cluster the LOT2 ITS sequences into OTUs  ->  LOT2_otu_map.csv
#
# Taxon identity in the source workbook is a BLAST-derived NAME
# STRING, which both over-splits (near-identical sequences under two
# spellings) and over-lumps (one label covering "26 spp."). Clustering
# the sequences themselves gives a reproducible analytical unit.
#
# Sequences live in sheet `Feuil2` of LOT2_samples.xlsx, one per
# isolate, keyed by Hofstetter culture code.
#
# This step is deterministic and only needs re-running when the input
# sequences change, so its output is committed and the R pipeline just
# reads it. That also keeps the pipeline runnable on Windows, where
# vsearch is usually unavailable.
#
# Requires vsearch:  conda create -n lot2 -c bioconda vsearch
# Usage:             bash cluster_otus.sh [path/to/vsearch]
# =============================================================
set -euo pipefail

VSEARCH="${1:-vsearch}"
command -v "$VSEARCH" >/dev/null 2>&1 || {
  echo "vsearch not found. Install it with:" >&2
  echo "  conda create -n lot2 -c bioconda vsearch" >&2
  echo "then pass its path: bash cluster_otus.sh /path/to/vsearch" >&2
  exit 1; }

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# 1. Export sequences from the workbook to FASTA
python3 - "$WORK/its.fasta" <<'PY'
import openpyxl, re, sys
ws = openpyxl.load_workbook('LOT2_samples.xlsx', data_only=True)['Feuil2']
n = 0
with open(sys.argv[1], 'w') as f:
    for r in range(1, ws.max_row + 1):
        code, seq = ws.cell(r, 2).value, ws.cell(r, 4).value
        if code is None or seq is None:
            continue
        s = re.sub(r'[^ACGTN]', '', str(seq).upper())
        if s:
            f.write(">%s\n%s\n" % (str(code).strip(), s))
            n += 1
print("exported %d ITS sequences" % n)
PY

# 2. Greedy centroid clustering at both conventional ITS thresholds.
#    97%   - the long-standing fungal OTU convention
#    98.5% - UNITE species-hypothesis default; kept for sensitivity
for ID in 0.97 0.985; do
  "$VSEARCH" --cluster_size "$WORK/its.fasta" --id "$ID" --iddef 2 \
             --strand both --uc "$WORK/uc_$ID.txt" --quiet
  echo "  id=$ID -> $(grep -c '^C' "$WORK/uc_$ID.txt") clusters"
done

# 3. Collapse the two .uc files into one culture_code -> OTU table
python3 - "$WORK/uc_0.97.txt" "$WORK/uc_0.985.txt" LOT2_otu_map.csv <<'PY'
import csv, sys

def load(path):
    m = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if f[0] in ("S", "H"):
            m[f[8]] = int(f[1])
    return m

a, b = load(sys.argv[1]), load(sys.argv[2])
with open(sys.argv[3], "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["culture_code", "otu_97", "otu_985"])
    for code in sorted(a, key=lambda x: (len(x), x)):
        w.writerow([code, "OTU_%03d" % (a[code] + 1),
                    "OTU985_%03d" % (b.get(code, -1) + 1)])
print("wrote %s (%d isolates, %d OTUs at 97%%)"
      % (sys.argv[3], len(a), len(set(a.values()))))
PY
