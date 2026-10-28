# Fungal Endophyte Community Analysis — LOT2 (Peru)

Culture-based fungal endophytes from two interlocked trees: a **strangler *Ficus*** (leaves; wood from trunk zones 1–5 and canopy branches, zone 6) and its **strangled Lauraceae host** (exposed branches, zone 6 only — the trunk is encased by the *Ficus*).

**All results (every table and plot): [`REPORT_all_results.html`](REPORT_all_results.html)** (Markdown version: `REPORT_all_results.md`). Keep the HTML next to `results/`; images are linked, not embedded.

## Input data

| File | Content |
|---|---|
| `LOT2_samples.xlsx` | One row per isolate (sheet 1): substrate, zone, sampling unit, reconciled orientation, ITS name at **100%** and **99%** BLAST similarity, GenBank hit. Sheet `Feuil2`: ITS sequences. |
| `LOT2_pooled_counts_ITS100.xlsx` | Per-taxon counts + taxonomy, 100% names (V. Hofstetter) |
| `LOT2_pooled_counts_ITS99.xlsx` | Per-taxon counts + taxonomy, 99% names (V. Hofstetter, **in production**) |
| `LOT2_otu_map.csv` | Isolate → 97% OTU, produced by `cluster_otus.sh` |

An optional **`Lifestyle`** column (per isolate in `LOT2_samples.xlsx`, or per taxon in the ITS100/ITS99 workbooks) enables the lifestyle analyses (pathogenic share, strangled Lauraceae vs strangler *Ficus*). Without it they are skipped.

## Taxonomic unit sets (all analyses run in parallel)

| Set | Unit |
|---|---|
| **ITS100** (reference) | BLAST top-score name at 100% ITS similarity |
| **ITS99** | BLAST top-score name at 99% ITS similarity (in production) |
| **OTU97** | 97% ITS OTUs: vsearch greedy centroid clustering (`--cluster_size --id 0.97 --iddef 2 --strand both`) of the `Feuil2` sequences, via `cluster_otus.sh` |

## Analysis rules

- Ranks used: phylum, subphylum, class, order, family, genus, species (genus + epithet). **Superclass and subclass are not used**; each rank is aggregated on its own value, and phylum/subphylum are aligned with the class.
- **Ficus trunk and branch wood form one substrate ("Ficus wood")**. Kept within it: the zone 1–6 height gradient. Kept across substrates: the zone-6 (top canopy) comparison.
- *Incertae sedis* (unresolved) isolates are shown in abundance/pie/relative-abundance charts and excluded from diversity and multivariate tests, rank by rank.
- Two curator corrections are applied on reading (Lauraceae branches all in zone 6; Ficus trunk units S11–S15, one per zone 5→1).
- Only two trees were sampled: host contrasts describe these two individuals, not the two plant taxa.

## Running

```r
source("main.R")   # from the project directory
```

| Script | Purpose |
|---|---|
| `main.R` | Parses inputs, runs every analysis once per unit set, then the comparison and the report |
| `parse_LOT2.R` | Reads the workbooks, name harmonisation, lineage, design checks → `results/data_checks/` |
| `plot_abundance.R`, `plot_abundance_pooled.R`, `plot_pie.R` | Abundance bar charts and pie charts per rank |
| `community_analysis.R` | Alpha/beta diversity, NMDS, PERMANOVA, ANOSIM, PERMDISP, IndVal, zone and orientation analyses |
| `diversity_analysis.R` | Observed vs estimated richness (Chao1/ACE), coverage, per-unit alpha diversity tests |
| `lifestyle_analysis.R` | Lifestyle vs substrate/tree (needs a `Lifestyle` column) |
| `compare_units.R` | ITS100 / ITS99 / OTU97 side by side → `results/comparison/` |
| `generate_report.R` | Writes `REPORT_all_results.md` / `.html` |
| `cluster_otus.sh` | Re-clusters the 97% OTUs (needs vsearch + python openpyxl); only if sequences change |

Outputs: `results/<ITS100|ITS99|OTU97>/{tables,plots}`, `results/comparison/`, `results/data_checks/`.
