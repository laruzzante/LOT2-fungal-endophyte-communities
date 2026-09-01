# Fungal Endophyte Community Analysis — LOT2 (Peru)

> **Analysis environment:** R 4.6.1 | Generated on 2026-09-01

## Overview

This project analyses the fungal endophyte communities isolated from three substrates collected in **Peru** (LOT2 sampling campaign):

| Substrate | Isolates (pooled) | Sample-level replicates (units) |
|-----------|:-:|:-:|
| **Ficus leaves** | 264 | 10 |
| **Ficus wood** | 167 | 12 |
| **Lauraceae leaves** | 219 | 8 |

### Sampling design

Samples were collected at different **tree zones** (heights):

- **Zones 1-5**: Tree trunk (zone 1 = lowest, zone 5 = highest)
- **Zone 6**: Canopy branches (highest zone)

Wood collected from zone 6 corresponds to **branch wood** (no trunk present). Wood from zones 1-5 is **trunk wood**.

| Substrate | Zones sampled | Notes |
|-----------|:---:|---|
| Ficus leaves | 6 only | All leaf samples from canopy |
| Ficus wood | 1-6 | Trunk (zones 1-5) and branch (zone 6) |
| Lauraceae leaves | 5-6 | Predominantly zone 6 |

### Handling of uncertain taxonomy

Entries with `NA`, `"?"`, `"NO"`, or `"incertae sedis"` are pooled into **"Incertae sedis"** before aggregation.

---

## Input Data

- **`LOT2_pooled_counts.xlsx`** (first sheet) — Pooled genotype counts per substrate with full taxonomy (284 genotypes)
- **`LOT2_samples.xlsx`** (first sheet) — Individual isolate records with sampling zone and unit (650 isolates)

## Scripts

| Script | Purpose |
|--------|---------|
| `main.R` | Master script — loads libraries, sources all other scripts in order |
| `parse_LOT2.R` | Reads both Excel files (first sheet only), standardises column names |
| `plot_abundance.R` | Horizontal bar charts of isolate counts by taxonomic level |
| `plot_abundance_pooled.R` | Bar charts with uncertain taxa pooled as "Incertae sedis" |
| `plot_pie.R` | Pie charts of community composition (3 pies per level) |
| `community_analysis.R` | Full community ecology analyses across all comparisons |
| `generate_readme.R` | Generates this README dynamically from analysis results |

```r
# From the project directory:
source("main.R")
```

### R packages used

| Package | Version | Role |
|---------|---------|------|
| readxl | 1.5.0 | Reading Excel input |
| dplyr | 1.2.1 | Data manipulation |
| tidyr | 1.3.2 | Data reshaping |
| ggplot2 | 4.0.3 | Plotting |
| scales | 1.4.0 | Axis formatting |
| vegan | 2.7.5 | Diversity, NMDS, PERMANOVA, ANOSIM, betadisper |
| indicspecies | 1.8.0 | Indicator species analysis (IndVal) |
| ggVennDiagram | 1.5.7 | Venn diagrams |

---

## Output Structure

```
plots/
  abundance/                        # Bar charts (raw)
  abundance_pooled_incertae_sedis/  # Bar charts (uncertain taxa pooled)
  pie_charts/                       # Pie charts
  community/                        # Community ecology plots (NMDS, rarefaction, etc.)
  png/                              # PNG versions for README embedding

tables/                             # Statistical results and summary tables
```

---

# Results

## A. Alpha Diversity (from pooled counts)

| Index | What it measures | Higher means |
|-------|-----------------|--------------|
| **Richness (S)** | Number of distinct taxa | More taxa present |
| **Shannon (H')** | Combines richness and evenness | More diverse |
| **Simpson (1-D)** | Probability two random individuals differ | More diverse |
| **Inverse Simpson** | Effective number of equally-common species | More even |
| **Pielou (J')** | How evenly individuals are distributed | More even |

### At ITS taxon level

| Substrate | S | N | H' | 1-D | Inv. Simp. | J' |
|-----------|:-:|:-:|:-:|:-:|:-:|:-:|
| Lauraceae leaves | 106 | 219 | 4.2745 | 0.9779 | 45.289 | 0.9166 |
| Ficus leaves | 135 | 264 | 4.4705 | 0.9791 | 47.8681 | 0.9114 |
| Ficus wood | 81 | 167 | 4.0073 | 0.9722 | 35.9858 | 0.9119 |

Full data: `tables/alpha_diversity.csv`

![Shannon diversity](plots/png/alpha_shannon_H.png)

---

## B. Community Composition

### B1. Rank-Abundance Curves

![Rank-abundance at genus level](plots/png/rank_abundance_genus.png)

### B2. Relative Abundance

![Relative abundance by phylum](plots/png/rel_abundance_phylum.png)

![Relative abundance by genus](plots/png/rel_abundance_genus.png)

### B3. Venn Diagrams — Shared Taxa

At **genus level** (44 total genera): **7** shared across all 3 substrates, **2** unique to Lauraceae leaves, **10** unique to Ficus leaves, **20** unique to Ficus wood.

![Venn diagram — genus](plots/png/venn_genus.png)

---

# Multivariate Statistical Analyses

Sample-level analyses use **32 samples** and **284 ITS taxa**, with **Bray-Curtis dissimilarity** and **999 permutations**.

---

## C1. Substrate Comparison: Ficus Leaves vs Ficus Wood vs Lauraceae Leaves

### NMDS Ordination

**Stress = 0.0771** (< 0.10 = good; < 0.20 = acceptable)

![NMDS — all substrates](plots/png/substrate_all_nmds.png)

![NMDS with species overlay](plots/png/substrate_all_nmds_species.png)

### PERMANOVA

```
PERMANOVA — Bray-Curtis distance
Formula: community ~ substrate
Permutations: 999

Permutation test for adonis under reduced model
Permutation: free
Number of permutations: 999

adonis2(formula = form, data = meta, permutations = 999, method = "bray")
         Df SumOfSqs      R2      F Pr(>F)    
Model     2    2.792 0.20041 3.6342  0.001 ***
Residual 29   11.140 0.79959                  
Total    31   13.932 1.00000                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

F = 3.6342, R² = 0.20041, p = 0.001

### ANOSIM

R = 0.5604 , p = 0.001 

### Beta-dispersion (PERMDISP)

F = 4.332, p = 0.02

### Pairwise PERMANOVA

| Comparison | F | R² | p (raw) | p (Holm) | Sig. |
|-----------|:-:|:-:|:-:|:-:|:-:|
| Ficus leaves vs Ficus wood | 3.8123 | 0.1536 | 0.001 | 0.003 | ** |
| Ficus leaves vs Lauraceae leaves | 3.3294 | 0.1638 | 0.001 | 0.003 | ** |
| Ficus wood vs Lauraceae leaves | 3.6839 | 0.1555 | 0.001 | 0.003 | ** |

### Rarefaction

![Rarefaction — all substrates](plots/png/substrate_all_rarefaction.png)

### Indicator Species (23 significant, p ≤ 0.05)

Full results: `tables/indicator_species_significant.csv`

---

## C2. Ficus Leaves vs Lauraceae Leaves

PERMANOVA: F = 3.3294, R² = 0.16377, p = 0.001

![NMDS — Ficus vs Lauraceae leaves](plots/png/substrate_leaves_nmds.png)

Full results in `tables/substrate_leaves_*.txt`

---

## E. Zone-Based Analyses

### E1. Ficus Wood Across Zones

Community comparison of Ficus wood isolates collected at different tree zones (1-6).

![NMDS — Ficus wood zones](plots/png/ficus_wood_zones_nmds.png)

Full results in `tables/ficus_wood_zones_*.txt`

### E2. Ficus Wood: Trunk (Zones 1-5) vs Branch (Zone 6)

PERMANOVA: F = 1.9802, R² = 0.15256, p = 0.004

![NMDS — trunk vs branch](plots/png/ficus_wood_trunk_vs_branch_nmds.png)

### E4. Substrate Comparison at Zone 6 (Branch Level Only)

PERMANOVA: F = 3.3279, R² = 0.22443, p = 0.001

![NMDS — substrates at zone 6](plots/png/substrate_zone6_branch_nmds.png)

### E5. Substrate × Position Interaction

Combined factor analysis testing whether community composition differs by substrate-position combinations.

![NMDS — substrate × position](plots/png/substrate_x_position_nmds.png)

---

## Abundance Plots

### `plots/abundance/`

Horizontal grouped bar charts showing absolute isolate counts per taxonomic level, broken down by substrate.

![Abundance by genus (pooled)](plots/png/abundance_by_genus.png)

### `plots/pie_charts/`

![Pie chart by phylum](plots/png/pie_by_phylum.png)

![Pie chart by genus](plots/png/pie_by_genus.png)

---

*Auto-generated on 2026-09-01 by `generate_readme.R`*
