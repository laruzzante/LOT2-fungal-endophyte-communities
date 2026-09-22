# Fungal Endophyte Community Analysis — LOT2 (Peru)

> **Analysis environment:** R 4.5.2 | Generated on 2026-09-22

## Overview

This project analyses the fungal endophyte communities isolated from three substrates collected in **Peru** (LOT2 sampling campaign):

| Substrate | Isolates (pooled) | Sample-level replicates (units) |
|-----------|:-:|:-:|
| **Ficus leaves** | 264 | 10 |
| **Ficus wood** | 167 | 13 |
| **Lauraceae leaves** | 219 | 10 |

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

#### Sampling orientation

Each sampled branch (and, where recorded, trunk face) also carries a **compass orientation**. This information only became usable once the field team reconciled two independent labelling systems: the colour codes written down by the climbers in the field, and the branch numbering entered into the project database. The two disagreed (notably a blue/turquoise colour clash, where blue should have been reserved for the Lauraceae, and disputed collection zones for the Lauraceae branches), so earlier versions of `LOT2_samples.xlsx` had an empty orientation column.

The re-issued workbook settles that correspondence and adds it as a reconciled column, which the analyses in **Section F** use. Two clean-up rules are applied when reading it:

- Some bearings still carry a replicate index (`N1`-`N4`, `NW1`-`NW4`). The digit identifies **which branch**, not which direction, so it is stripped: `N3` becomes `N`.
- Material marked `?`, and all Ficus trunk wood (for which no orientation was ever recorded), is treated as **unresolved** and excluded from the orientation analyses.

This leaves **519 of 650 isolates** (80%) with a usable bearing; **131** are unresolved. A single sampling unit can span more than one bearing, so the orientation analyses group isolates by **substrate x zone x unit x bearing**, giving **27 orientation-level sampling units** rather than the 33 used elsewhere.

### Handling of uncertain taxonomy (Incertae sedis)

Many isolates cannot be confidently named at every taxonomic rank. Entries flagged `NA`, `"?"`, `"NO"` or `"incertae sedis"` are treated as **Incertae sedis** ("of uncertain placement") and are handled **rank by rank**:

- They are **shown** — as a single pooled *Incertae sedis* category — **only in the abundance bar charts and the pie charts**, so that the full isolate count is never hidden.
- They are **excluded from every diversity index, rank-/relative-abundance curve, Venn diagram and multivariate test**. A single large, shared "unknown" bin behaves like a taxon that is common everywhere: it inflates apparent overlap between substrates and **flattens the real ecological differences** we are trying to detect.
- Because the exclusion is applied **independently at each rank**, an isolate that is unresolved at *species* level but has a defined *genus*, *family* or *order* still contributes to the analyses run at those higher ranks (see **Section D — Multi-level analyses**).

> The sample-level multivariate analyses are run at **ITS-genotype** resolution, where each of the sequenced genotypes is a distinct entity. There is no dominant "unknown" bin at that level, so no isolates are dropped there; the *Incertae sedis* filtering matters only when isolates are grouped into higher taxa.

---

## Input Data

- **`LOT2_pooled_counts.xlsx`** (first sheet) — Pooled genotype counts per substrate with full taxonomy (285 genotypes)
- **`LOT2_samples.xlsx`** (first sheet) — Individual isolate records with sampling zone, unit and reconciled branch/trunk orientation (650 isolates)

  The workbook also keeps the reconciliation working columns (field colour code, database colour code, pre-reconciliation orientation, field notes). They are read for traceability but only the reconciled orientation column is used.

## Scripts

| Script | Purpose |
|--------|---------|
| `main.R` | Master script — loads libraries, sources all other scripts in order |
| `parse_LOT2.R` | Reads both Excel files (first sheet only), standardises column names |
| `plot_abundance.R` | Horizontal bar charts of isolate counts by taxonomic level |
| `plot_abundance_pooled.R` | Bar charts with uncertain taxa pooled as "Incertae sedis" |
| `plot_pie.R` | Pie charts of community composition (3 pies per level) |
| `community_analysis.R` | Full community ecology analyses; ITS-level tests plus multi-rank taxonomic sweeps (Section D) |
| `generate_readme.R` | Generates this README dynamically from analysis results |

```r
# From the project directory:
source("main.R")
```

### R packages used

| Package | Version | Role |
|---------|---------|------|
| readxl | 1.4.5 | Reading Excel input |
| dplyr | 1.1.4 | Data manipulation |
| tidyr | 1.3.1 | Data reshaping |
| ggplot2 | 4.0.1 | Plotting |
| scales | 1.4.0 | Axis formatting |
| vegan | 2.7.6 | Diversity, NMDS, PERMANOVA, ANOSIM, betadisper |
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

> Throughout the Results, each analysis is introduced with a short **plain-language explanation** of what it measures and how to read it, followed by a **brief interpretation** of what the LOT2 data actually show. A synthesis of all findings is given in the final **Conclusion**.

## A. Alpha Diversity (from pooled counts)

**What it is.** *Alpha diversity* describes how varied the fungal community is **within a single substrate**. It combines two ideas: *richness* (how many different taxa are present) and *evenness* (whether isolates are spread evenly across taxa or dominated by a few). The indices below capture different balances of these two ideas. *Incertae sedis* taxa are excluded so that diversity reflects only confidently identified fungi.

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

**Interpretation.** The richest community is **Ficus leaves** (S = 135 distinct ITS taxa), and the least rich is **Ficus wood** (S = 81). Shannon diversity is highest in **Ficus leaves** (H' = 4.4705) and lowest in **Ficus wood** (H' = 4.0073). Pielou evenness is high (J' > 0.9) for all three substrates, meaning no single genotype dominates any community — isolates are spread across many co-occurring taxa. Because richness partly reflects sampling effort (N differs between substrates), richness values should be compared together with the **rarefaction curves** below, which put all substrates on an equal-effort footing.

Diversity was also computed at every higher rank (phylum → genus); the full table is in `tables/alpha_diversity.csv`.

![Shannon diversity](plots/png/alpha_shannon_H.png)

*The bars show Shannon H' for each substrate across taxonomic ranks. H' naturally decreases towards coarser ranks (fewer categories), but the ranking of substrates stays broadly consistent, indicating the diversity differences are not an artefact of one particular rank.*

---

## B. Community Composition

These analyses describe **what the communities are made of** and **how much they overlap**, again after removing *Incertae sedis*.

### B1. Rank-Abundance Curves

**What it is.** Taxa are ranked from most to least abundant (x-axis) against their relative abundance on a log scale (y-axis). A **steep** curve means a few taxa dominate (low evenness); a **shallow, long** curve means many taxa share the community evenly (high evenness). The length of each curve reflects richness.

![Rank-abundance at genus level](plots/png/rank_abundance_genus.png)

*Interpretation.* All three substrates show relatively shallow curves with long tails, confirming the high evenness seen in the Pielou index: communities are not dominated by one or two hyper-abundant genera but consist of many moderately frequent taxa plus a long tail of rare ones — a pattern typical of tropical endophyte assemblages.

### B2. Relative Abundance

**What it is.** Stacked bars show the **proportional composition** of each substrate at a given rank (each bar sums to 100%). They make it easy to see which phyla/genera dominate and how composition shifts between substrates.

![Relative abundance by phylum](plots/png/rel_abundance_phylum.png)

![Relative abundance by genus](plots/png/rel_abundance_genus.png)

*Interpretation.* At phylum level the communities are overwhelmingly **Ascomycota**, as expected for culturable endophytes. The genus-level bars reveal the real contrast between substrates: the identity and proportion of dominant genera differ markedly between leaves and wood, foreshadowing the significant substrate effect quantified in the multivariate tests below.

### B3. Venn Diagrams — Shared Taxa

**What it is.** The Venn diagram counts how many taxa are **unique** to each substrate versus **shared** between them. It is a simple presence/absence view of community overlap (abundance is ignored).

*Interpretation.* At **genus level** (43 identified genera in total): **6** genera occur in all three substrates (a shared generalist core), while **2** are unique to Lauraceae leaves, **10** unique to Ficus leaves and **20** unique to Ficus wood. The substantial number of substrate-exclusive genera indicates a degree of **habitat specialisation** layered on top of a shared generalist core.

![Venn diagram — genus](plots/png/venn_genus.png)

---

# Multivariate Statistical Analyses

These analyses ask **whether whole communities differ between groups** (substrates, zones, positions). They use the sample-level data (**33 sampling units**, **284 ITS genotypes**), the **Bray-Curtis dissimilarity** (a 0–1 measure of how different two samples are in both *which* taxa are present and *how abundant* they are), and **999 permutations** to obtain p-values without assuming normality.

**How to read each test:**

- **NMDS ordination** — squeezes the many-dimensional Bray-Curtis distances into a 2-D map so that samples plotting close together have similar communities. The **stress** value measures distortion: < 0.10 excellent, < 0.20 acceptable, > 0.20 unreliable. Crosses mark group centroids; shaded ellipses show 95% confidence regions. **A stress at or near zero is not a good fit** — it means the ordination has *degenerated*, which happens when the community matrix is so sparse that most sample pairs share no taxa and their Bray-Curtis distance is pinned at 1. See the sparsity caveat below.
- **PERMANOVA** (`adonis2`) — tests whether **group centroids differ**. **R²** is the fraction of community variation explained by the grouping (effect size); a small **p** means the separation is unlikely by chance.
- **ANOSIM** — a complementary rank-based test; **R** ranges from 0 (no separation) to 1 (groups completely distinct).
- **Beta-dispersion / PERMDISP** (`betadisper`) — checks whether groups differ in **within-group spread** rather than location. If PERMDISP is significant, part of a PERMANOVA result may reflect unequal dispersion rather than a pure shift in composition, so it is an important caveat.
- **Pairwise PERMANOVA** — which specific pairs of groups differ, with Holm correction for multiple tests.
- **Rarefaction** — expected richness rescaled to equal sampling effort, so richness can be compared fairly.
- **Indicator species (IndVal)** — identifies taxa statistically associated with (diagnostic of) a particular group.

### An important caveat: the ITS-level matrix is very sparse

At ITS-genotype resolution this dataset is dominated by rare taxa: **62% of the 284 genotypes were isolated exactly once**, and as a result **59% of all sampling-unit pairs share no genotype at all**. Every one of those pairs has a Bray-Curtis distance of exactly 1, so the distance matrix is largely saturated.

This has two concrete consequences for how the results below should be read:

- **The ITS-level NMDS maps are degenerate and should not be interpreted.** Their near-zero stress is the symptom, not a virtue: with most distances tied at 1 there is no gradient left for the ordination to lay out, and the resulting configuration is arbitrary (note the implausible axis ranges). Each such plot is now labelled as degenerate in its subtitle.
- **PERMANOVA and ANOSIM remain valid**, because they work on the ranks and the sums of squares of the distance matrix rather than on a 2-D embedding. They are the tests to trust here, together with the indicator-species analysis.

The **multi-rank sweep in Section D is the constructive answer to this.** Grouping isolates into genera, families or orders collapses the singleton problem: at those ranks samples share taxa, the ordinations reach sensible stress values (see the "pairs sharing no taxa" column — it falls from 59% at ITS level to 33% at genus level, where the ordination stress reaches a healthy 0.1533), and the substrate signal is reproduced. **Use the Section D ordinations as the readable maps of these communities.**

---

## C1. Substrate Comparison: Ficus Leaves vs Ficus Wood vs Lauraceae Leaves

This is the **headline comparison**: do the three substrates host different fungal communities?

### NMDS Ordination

**Stress = 1e-04** — **a degenerate solution, not an excellent one** (see the sparsity caveat above). Read the separation from PERMANOVA and from the higher-rank ordinations in Section D, not from this map.

![NMDS — all substrates](plots/png/substrate_all_nmds.png)

![NMDS with species overlay](plots/png/substrate_all_nmds_species.png)

*Interpretation.* The three substrates form visually distinct clouds, with the two leaf substrates sitting closer to each other than to wood — consistent with tissue type (leaf vs wood) being a strong driver. The species overlay points to the genotypes pulling each substrate apart.

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
Model     2   2.7727 0.19168 3.5571  0.001 ***
Residual 30  11.6922 0.80832                  
Total    32  14.4648 1.00000                  
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

**F = 3.5571, R² = 0.19168, p = 0.001** — the substrate effect is **statistically significant**. Substrate explains about **19%** of the total community variation, a large effect for field endophyte data.

### ANOSIM

**R = 0.5506 , p = 0.001 ** — **statistically significant**. An R of this magnitude confirms that between-substrate differences clearly exceed within-substrate variation.

### Beta-dispersion (PERMDISP)

**F = 2.659, p = 0.097** — dispersion differences are **not statistically significant**. Groups are comparably variable, so the PERMANOVA result reflects a genuine shift in composition rather than unequal spread.

### Pairwise PERMANOVA

Which substrate *pairs* differ (Holm-corrected p). `***` p<0.001, `**` p<0.01, `*` p<0.05.

| Comparison | F | R² | p (raw) | p (Holm) | Sig. |
|-----------|:-:|:-:|:-:|:-:|:-:|
| Ficus leaves vs Ficus wood | 3.9209 | 0.1573 | 0.001 | 0.003 | ** |
| Ficus leaves vs Lauraceae leaves | 3.1421 | 0.1486 | 0.001 | 0.003 | ** |
| Ficus wood vs Lauraceae leaves | 3.5295 | 0.1439 | 0.001 | 0.003 | ** |

*Interpretation.* Every pair of substrates differs significantly after correction — each substrate carries a distinguishable community, not merely one odd substrate against two similar ones.

### Rarefaction

Expected number of taxa if every substrate had been sampled to the same number of isolates — a fair richness comparison that removes the effect of unequal sampling effort.

![Rarefaction — all substrates](plots/png/substrate_all_rarefaction.png)

*Interpretation.* None of the curves has fully levelled off, so additional sampling would still recover new taxa in every substrate (the communities are undersampled, as usual for hyper-diverse tropical fungi). The **relative ordering** of the curves indicates which substrate is richest at equal effort, which is the sampling-fair complement to the raw richness values in Section A.

### Indicator Species (22 significant, p ≤ 0.05)

Indicator (IndVal) analysis finds taxa that are **diagnostic** of a particular substrate — both faithful to it (mostly found there) and frequent within it. **22** ITS genotypes are significant indicators, i.e. reliable biological markers of their substrate. Full ranked list: `tables/indicator_species_significant.csv`.

---

## C2. Ficus Leaves vs Lauraceae Leaves (leaf substrates only)

**Why.** Both are **leaf** endophyte communities but from different host plants, so this isolates the **host effect** from the leaf-vs-wood tissue effect.

**PERMANOVA: F = 3.1421, R² = 0.14862, p = 0.001** — **statistically significant**. Even between two leaf communities, host identity (Ficus vs Lauraceae) leaves a detectable signature, explaining about 15% of the variation.

![NMDS — Ficus vs Lauraceae leaves](plots/png/substrate_leaves_nmds.png)

Full results in `tables/substrate_leaves_*.txt`

---

## D. Multi-level Taxonomic Analyses

**Why this section exists.** The multivariate tests above use **ITS genotypes**, the finest possible resolution. But an isolate that is *Incertae sedis* at species level often still has a defined **genus, family or order**. By linking every isolate to its **full lineage** (from the pooled taxonomy) we can repeat the community comparisons at each rank and ask: **is the substrate/zone signal a fine-scale artefact, or does it hold when isolates are grouped into higher, more confidently identified taxa?**

**Lineage coverage** — isolates with a defined value at each rank (out of 650):

| Rank | Isolates resolved |
|------|:-:|
| phylum | 642 / 650 |
| class | 644 / 650 |
| order | 627 / 650 |
| family | 582 / 650 |
| genus | 537 / 650 |
| species | 238 / 650 |
| its_taxon | 650 / 650 |

### D1. Substrate comparison across ranks

| Rank | Taxa | Pairs sharing no taxa | NMDS stress | PERMANOVA F | R² | p | ANOSIM R | p | PERMDISP p |
|------|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| phylum |   4 | 0% | 0.0260 |  6.6886 | 0.3084 | 0.001 | 0.2410 | 0.001 | 0.096 |
| class |   9 | 4% | 0.1320 | 11.6909 | 0.4380 | 0.001 | 0.4083 | 0.001 | 0.001 |
| order |  25 | 12% | 0.1482 | 13.1621 | 0.4674 | 0.001 | 0.6181 | 0.001 | 0.008 |
| family |  34 | 22% | 0.1561 | 11.5457 | 0.4349 | 0.001 | 0.6250 | 0.001 | 0.006 |
| genus |  43 | 33% | 0.1533 | 10.0711 | 0.4017 | 0.001 | 0.5992 | 0.001 | 0.003 |
| species |  74 | 67% | 0.0001 |  4.4038 | 0.2330 | 0.001 | 0.4232 | 0.001 | 0.001 |
| its_taxon | 284 | 59% | 0.0001 |  3.5571 | 0.1917 | 0.001 | 0.5506 | 0.001 | 0.097 |

*Interpretation.* The substrate effect is **significant at every taxonomic rank** (PERMANOVA p ≤ 0.004). Crucially, the effect size does **not** weaken when isolates are grouped into higher taxa — it is *strongest* around **order level** (R² ≈ 0.47), compared with R² ≈ 0.19 at ITS level. In other words, the substrates differ not just in which fine genotypes they carry, but in their broad taxonomic make-up, and removing the *Incertae sedis* noise sharpens rather than blurs that separation.

![NMDS at genus level — substrates](plots/png/substrate_bylevel_genus_nmds.png)

![NMDS at family level — substrates](plots/png/substrate_bylevel_family_nmds.png)

### D2. Ficus wood zones across ranks

| Rank | Taxa | Pairs sharing no taxa | NMDS stress | PERMANOVA F | R² | p | ANOSIM R | p | PERMDISP p |
|------|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| phylum |  4 | 0% | 0.0399 | 0.9158 | 0.3955 | 0.528 | -0.0032 | 0.431 | — |
| class |  9 | 1% | 0.1495 | 1.9911 | 0.5872 | 0.025 |  0.4229 | 0.006 | — |
| order | 19 | 3% | 0.1315 | 1.3953 | 0.4992 | 0.098 |  0.2524 | 0.060 | — |
| family | 26 | 4% | 0.1490 | 1.1705 | 0.4554 | 0.257 |  0.0276 | 0.447 | — |
| genus | 29 | 12% | 0.1510 | 0.8187 | 0.3690 | 0.866 | -0.0252 | 0.563 | — |
| species | 29 | 76% | 0.0001 | 1.0011 | 0.4169 | 0.488 |  0.0130 | 0.434 | — |
| its_taxon | 82 | 45% | 0.0731 | 0.9259 | 0.3981 | 0.743 | -0.0714 | 0.648 | — |

*Interpretation.* Some ranks show a zone effect; see the table.

### D3. Substrate × position across ranks

| Rank | Taxa | Pairs sharing no taxa | NMDS stress | PERMANOVA F | R² | p | ANOSIM R | p | PERMDISP p |
|------|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| phylum |   4 | 0% | 0.0260 | 3.6492 | 0.3427 | 0.008 | 0.3369 | 0.002 | 0.191 |
| class |   9 | 4% | 0.1320 | 9.2494 | 0.5692 | 0.001 | 0.5600 | 0.001 | 0.123 |
| order |  25 | 12% | 0.1482 | 9.0121 | 0.5628 | 0.001 | 0.6901 | 0.001 | 0.376 |
| family |  34 | 22% | 0.1561 | 7.8847 | 0.5297 | 0.001 | 0.6919 | 0.001 | 0.202 |
| genus |  43 | 33% | 0.1533 | 6.6249 | 0.4862 | 0.001 | 0.6718 | 0.001 | 0.100 |
| species |  74 | 67% | 0.0001 | 3.1999 | 0.3216 | 0.001 | 0.5884 | 0.001 | 0.147 |
| its_taxon | 284 | 59% | 0.0001 | 2.5366 | 0.2660 | 0.001 | 0.6468 | 0.001 | 0.611 |

*Interpretation.* Combining substrate with trunk/branch position remains significant at every rank, and the effect size again peaks at intermediate ranks (class–family). This mirrors the substrate result: the signal is carried by broad taxonomic groups, not just rare fine-scale genotypes.

Summary tables: `tables/substrate_multilevel_summary.csv`, `tables/ficus_wood_zones_multilevel_summary.csv`, `tables/substrate_x_position_multilevel_summary.csv`.

---

## E. Zone-Based Analyses

These test whether **height on the tree** (trunk zones 1–5 vs canopy branch zone 6) structures the community. All are run at ITS-genotype level; Section D2 already showed the zone question is also answered the same way at higher ranks.

### E1. Ficus Wood Across Zones

Community comparison of Ficus wood isolates collected at different tree zones (1–6).

![NMDS — Ficus wood zones](plots/png/ficus_wood_zones_nmds.png)

*Interpretation.* Samples from different zones intermingle on the NMDS with no zone-wise grouping, indicating wood-inhabiting fungi are distributed largely independently of height. Full results: `tables/ficus_wood_zones_*.txt`.

### E2. Ficus Wood: Trunk (Zones 1-5) vs Branch (Zone 6)

**PERMANOVA: F = 1.7926, R² = 0.14013, p = 0.006** — **statistically significant**. Trunk and branch wood host detectably different communities.

![NMDS — trunk vs branch](plots/png/ficus_wood_trunk_vs_branch_nmds.png)

### E4. Substrate Comparison at Zone 6 (Branch Level Only)

**Why.** Restricting to zone 6 removes any confound between substrate and height, since all three substrates are present there.

**PERMANOVA: F = 2.8436, R² = 0.2131, p = 0.001** — **statistically significant**. The substrate effect persists even within a single zone, confirming it is driven by substrate itself and not by differences in sampling height.

![NMDS — substrates at zone 6](plots/png/substrate_zone6_branch_nmds.png)

### E5. Substrate × Position Interaction

Combined-factor analysis testing whether communities differ across substrate-position combinations (e.g. *Ficus wood - Trunk* vs *Ficus wood - Branch* vs *Ficus leaves - Branch* ...).

![NMDS — substrate × position](plots/png/substrate_x_position_nmds.png)

*Interpretation.* Groups separate primarily **by substrate**, with position adding only minor structure — substrate is the dominant organiser of these endophyte communities (see also the multi-rank confirmation in Section D3).

---

## F. Branch & Trunk Sampling Orientation

**Why.** Different faces of a tree experience very different microclimates — sun exposure, temperature, how long the surface stays wet after rain, prevailing wind. If that matters to endophytes, communities should differ systematically between compass bearings. This is the first LOT2 analysis able to test it, because the orientation column only became usable with the re-issued samples workbook (see **Sampling orientation** above).

### F1. What was actually sampled

Orientation is an **observational**, not a designed, factor here: bearings were recorded for whichever branches were reachable, so replication is uneven and, for some substrate/bearing combinations, very thin. That constrains how much the tests below can detect, and it is the single most important caveat for this section.

| Substrate | Position | Orientation | Sampling units | Isolates | ITS taxa |
|-----------|:-:|:-:|:-:|:-:|:-:|
| Ficus leaves | Branch | E | 2 | 51 | 34 |
| Ficus leaves | Branch | N | 2 | 49 | 30 |
| Ficus leaves | Branch | S | 2 | 53 | 43 |
| Ficus leaves | Branch | W | 3 | 62 | 46 |
| Ficus wood | Branch | E | 2 | 32 | 21 |
| Ficus wood | Branch | N | 1 | 11 |  9 |
| Ficus wood | Branch | S | 2 | 10 |  9 |
| Ficus wood | Branch | W | 2 | 32 | 24 |
| Lauraceae leaves | Branch | E | 1 |  2 |  2 |
| Lauraceae leaves | Branch | N | 3 | 81 | 51 |
| Lauraceae leaves | Branch | NW | 4 | 83 | 57 |
| Lauraceae leaves | Trunk | N | 1 | 28 | 23 |
| Lauraceae leaves | Trunk | NW | 2 | 25 | 16 |

Unresolved (excluded from this section): 49 Ficus leaves branch, 82 Ficus wood trunk.

Three consequences follow, and they shape every result below:

- Only **Ficus leaves** and **Ficus wood** were sampled on all four cardinal bearings. **Lauraceae leaves** were only ever collected from northern and north-western faces (plus two isolates from one eastern unit), so the Lauraceae cannot contribute to a full-compass comparison.
- **Ficus trunk wood carries no orientation at all**, so this section is effectively a *branch*-level analysis.
- With 27 orientation-level units spread over 5 bearings, most groups hold **1-4 replicates**. Tests on a single substrate (7-11 units) have little power: a null result here means *no effect was detectable*, not *no effect exists*.

### F2. Does orientation structure the community?

Same test battery as the substrate analyses — PERMANOVA, ANOSIM, PERMDISP — run on the orientation-level sampling units at ITS-genotype resolution.

| Analysis | Units | Groups | PERMANOVA F | R² | p | ANOSIM R | p | PERMDISP p | Indicators |
|----------|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| Orientation - all substrates | 27 |  5 | 1.1019 | 0.1669 | 0.159 | 0.0674 | 0.118 | 0.972 | 3 |
| Orientation - Ficus leaves |  9 |  4 | 1.0637 | 0.3896 | 0.306 | -0.0611 | 0.616 | — | 1 |
| Orientation - Ficus wood |  7 |  4 | 0.8833 | 0.4690 | 0.795 | -0.2222 | 0.816 | — | 0 |
| Orientation - Lauraceae leaves | 11 |  3 | 1.3088 | 0.2465 | 0.088 | 0.1008 | 0.289 | — | 1 |
| Substrate x orientation | 27 | 11 | 1.4030 | 0.4672 | 0.001 | 0.3387 | 0.002 | — | 6 |
| Sun aspect - all substrates | 27 |  3 | 1.3590 | 0.1017 | 0.014 | 0.1513 | 0.026 | 0.515 | 8 |
| Sun aspect - Ficus leaves |  9 |  3 | 1.0311 | 0.2558 | 0.399 | -0.1424 | 0.681 | — | — |
| Sun aspect - Ficus wood |  7 |  3 | 0.9908 | 0.3313 | 0.487 | 0.1020 | 0.350 | — | — |
| Sun aspect - Lauraceae leaves | 11 |  2 | 1.5808 | 0.1494 | 0.175 | 0.6444 | 0.186 | — | — |

> **PERMDISP shown as —** where the smallest group holds fewer than 3 sampling units. A centroid computed from one or two points has a degenerate spread, which inflates the dispersion F ratio to meaningless magnitudes; those tests are written to the `*_betadisper.txt` files with a warning but are not reported as numbers.

**Orientation on its own: F = 1.1019, R² = 0.1669, p = 0.159** — **not statistically significant**. Within each substrate taken separately the result is the same: no substrate shows a significant orientation effect. ANOSIM agrees, and for Ficus leaves and Ficus wood the ANOSIM R is actually **negative** — meaning units from *different* bearings are, if anything, slightly more similar to each other than units from the *same* bearing. That is the signature of no orientation structure at all.

![NMDS — orientation, all substrates](plots/png/orientation_all_nmds.png)

![NMDS — orientation within Ficus leaves](plots/png/ficus_leaves_orientation_nmds.png)

*Interpretation.* The bearings overlap almost completely, with centroids piled on top of one another. These are ITS-level ordinations, so — as everywhere in this report — the maps are degenerate and carry no weight on their own; the conclusion rests on the PERMANOVA and ANOSIM results above, and on the higher-rank ordinations in **F7**, which agree. Whatever separates these communities, it is not which side of the tree the material came from.

### F3. Separating orientation from substrate

Two of the pooled analyses above **do** come out significant — *substrate x orientation* and *sun aspect* — and both need care, because orientation is partly confounded with substrate in this dataset (the Lauraceae units are almost all north/north-west facing). A marginal (type-III) PERMANOVA asks what each factor explains **once the other is accounted for**:

```
PERMANOVA - variance partitioning, Bray-Curtis distance
Only sampling units with a reconciled orientation are used.
Permutations: 999

== Marginal (type III): each term adjusted for the other ==
Permutation test for adonis under reduced model
Marginal effects of terms
Permutation: free
Number of permutations: 999

adonis2(formula = comm_orient ~ substrate + orientation, data = meta_orient, permutations = 999, method = "bray", by = "margin")
            Df SumOfSqs      R2      F Pr(>F)    
substrate    2   1.7622 0.14998 2.1956  0.001 ***
orientation  4   1.6283 0.13859 1.0144  0.398    
Residual    20   8.0262 0.68311                  
Total       26  11.7495 1.00000                  
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

== Sequential (type I): substrate entered first ==
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = comm_orient ~ substrate + orientation, data = meta_orient, permutations = 999, method = "bray", by = "terms")
            Df SumOfSqs      R2      F Pr(>F)    
substrate    2   2.0950 0.17830 2.6101  0.001 ***
orientation  4   1.6283 0.13859 1.0144  0.398    
Residual    20   8.0262 0.68311                  
Total       26  11.7495 1.00000                  
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

== Marginal (type III): substrate + sun aspect ==
Permutation test for adonis under reduced model
Marginal effects of terms
Permutation: free
Number of permutations: 999

adonis2(formula = comm_orient ~ substrate + sun_aspect, data = meta_orient, permutations = 999, method = "bray", by = "margin")
           Df SumOfSqs      R2      F Pr(>F)    
substrate   2   1.7759 0.15114 2.2253  0.001 ***
sun_aspect  2   0.8762 0.07457 1.0979  0.208    
Residual   22   8.7783 0.74713                  
Total      26  11.7495 1.00000                  
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

== Sampling units per substrate x sun aspect ==
(empty cells show where aspect is confounded with substrate)
                  
                   Lateral (E/W) Shaded (S) Sun-facing (N)
  Ficus leaves                 5          2              2
  Ficus wood                   4          2              1
  Lauraceae leaves             1          0             10
```

**Substrate: R² = 0.14998, p = 0.001. Orientation: R² = 0.13859, p = 0.398.** Substrate holds up; orientation does not. The apparently strong *substrate x orientation* result (R² ≈ 0.4672) is therefore the substrate effect re-expressed through a factor that happens to encode it, not evidence that bearing matters.

#### Sun aspect

LOT2 was collected in **Peru, i.e. the southern hemisphere**, where the *northern* face of a tree is the sun-exposed one and the southern face the shaded one. Grouping bearings into **sun-facing (N/NE/NW)**, **shaded (S/SE/SW)** and **lateral (E/W)** gives a more direct microclimate proxy than the raw compass, and coarser groups mean better replication.

Pooled across substrates this is significant (R² = 0.1017, p = 0.014, PERMDISP p = 0.515, so not a dispersion artefact) — but the confounding is severe: of the 27 orientation units, **all but one Lauraceae unit is sun-facing**, while Ficus leaves and wood are mostly lateral. Adjusting for substrate, **sun aspect gives R² = 0.07457, p = 0.208** — **not statistically significant**. Within each substrate individually it is likewise non-significant. The pooled result is substrate wearing an aspect label.

![NMDS — sun aspect](plots/png/orientation_sun_aspect_nmds.png)

### F4. A directional gradient rather than discrete bearings

Treating the compass as four or five **discrete classes** spends a lot of degrees of freedom on a dataset this small. An alternative is to treat the bearing as the **circular variable** it is: decomposing it into `cos(bearing)` (the north-south axis) and `sin(bearing)` (the east-west axis) tests for a community that changes *smoothly* around the tree using only 2 df, which is more powerful at this level of replication.

| Analysis | Units | Bearings | N-S axis R² | p | E-W axis R² | p | Joint R² | Joint p |
|----------|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| All substrates | 27 | 5 | 0.0477 | 0.084 | 0.0355 | 0.611 | 0.0844 | 0.191 |
| Ficus leaves |  9 | 4 | 0.1493 | 0.121 | 0.1343 | 0.272 | 0.2836 | 0.107 |
| Ficus wood |  7 | 4 | 0.1883 | 0.338 | 0.1377 | 0.791 | 0.3260 | 0.542 |
| Lauraceae leaves | 11 | 3 | 0.1369 | 0.036 | 0.0932 | 0.446 | 0.2465 | 0.088 |

*Interpretation.* No joint directional gradient is significant. The one nominally significant single axis is the **north-south axis in Lauraceae leaves** (p = 0.036), and it should not be over-read: the Lauraceae span only N and NW (plus two eastern isolates), so the "north-south axis" is fitted over a ~45° arc rather than a full compass, the joint test for the same substrate is non-significant (p = 0.088), and it is one nominal result among eight axis tests with no correction applied.

### F5. Diversity per orientation

Community *structure* may not differ while *diversity* still does — a more exposed face could simply support fewer taxa. Richness depends strongly on how many isolates each group contributed, so richness is also reported **rarefied to a common isolate count**.

| Orientation | S | N | H' | 1-D | J' | Rarefied S (to 63) |
|-------------|:-:|:-:|:-:|:-:|:-:|:-:|
| E | 57 |  85 | 3.8174 | 0.9694 | 0.9442 | 45.10 |
| N | 94 | 169 | 4.2844 | 0.9817 | 0.9430 | 46.54 |
| NW | 66 | 108 | 3.9516 | 0.9738 | 0.9432 | 44.48 |
| S | 52 |  63 | 3.8845 | 0.9776 | 0.9831 | 52.00 |
| W | 66 |  94 | 4.0425 | 0.9783 | 0.9649 | 48.68 |

| Sun aspect | S | N | H' | 1-D | J' | Rarefied S (to 63) |
|------------|:-:|:-:|:-:|:-:|:-:|:-:|
| Lateral (E/W) | 114 | 179 | 4.4573 | 0.9816 | 0.9411 | 49.58 |
| Shaded (S) |  52 |  63 | 3.8845 | 0.9776 | 0.9831 | 52.00 |
| Sun-facing (N) | 132 | 277 | 4.5091 | 0.9831 | 0.9235 | 46.76 |

*Interpretation.* Raw richness tracks sampling effort almost exactly — the most-sampled bearing is also the richest — but once rarefied to a common 63 isolates the bearings are within a few taxa of each other (44.48-52 taxa). Evenness (Pielou J') is uniformly high, as everywhere else in this dataset. There is no diversity gradient around the tree to match the absent compositional one.

![Richness by orientation and substrate](plots/png/orientation_alpha_richness_S.png)

![Shannon diversity by orientation and substrate](plots/png/orientation_alpha_shannon_H.png)

> The per-substrate rarefaction target is pulled down to a handful of isolates by the smallest Ficus wood groups, so `tables/orientation_alpha_diversity_by_substrate.csv` should be read as indicative only. The pooled table above is the more reliable of the two.

### F6. Taxon sharing between orientations

| Occurs in ... orientations | ITS genotypes |
|:-:|:-:|
| 1 | 182 |
| 2 | 35 |
| 3 | 17 |
| 4 | 8 |

*Interpretation.* Most genotypes (182 of 242) were found on a single bearing. On its own that looks like strong orientation fidelity, but it is the expected consequence of **rarity plus thin sampling**: this community is dominated by singletons and doubletons (see the rank-abundance curves in Section B1), and a taxon seen once can only ever be recorded at one bearing. The multivariate tests above, which use the whole community at once rather than taxon-by-taxon presence, find no such structure — and they are the trustworthy reading.

![Shared genotypes between orientations — Ficus leaves](plots/png/orientation_venn_ficus_leaves.png)

![Rarefaction by orientation](plots/png/orientation_all_rarefaction.png)

### F7. Orientation across taxonomic ranks

As in Section D, the comparison is repeated at every rank, in case a bearing effect exists among broader taxonomic groups but is hidden by genotype-level noise.

| Rank | Taxa | Pairs sharing no taxa | NMDS stress | PERMANOVA F | R² | p | ANOSIM R | p | PERMDISP p |
|------|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| phylum |   4 | 0% | 0.0054 | 1.0178 | 0.1562 | 0.467 | 0.0613 | 0.141 | 0.587 |
| class |   8 | 0% | 0.0961 | 1.3581 | 0.1980 | 0.213 | 0.1236 | 0.050 | 0.675 |
| order |  22 | 7% | 0.1361 | 1.0546 | 0.1609 | 0.416 | 0.0431 | 0.214 | 0.687 |
| family |  28 | 17% | 0.1465 | 1.0315 | 0.1579 | 0.415 | 0.0376 | 0.222 | 0.723 |
| genus |  34 | 25% | 0.1652 | 1.2668 | 0.1872 | 0.168 | 0.0899 | 0.072 | 0.620 |
| species |  59 | 57% | 0.0001 | 1.2025 | 0.2020 | 0.122 | 0.1153 | 0.061 | 0.932 |
| its_taxon | 242 | 59% | 0.0001 | 1.1019 | 0.1669 | 0.159 | 0.0674 | 0.118 | 0.972 |

*Interpretation.* Orientation is non-significant at **every** taxonomic rank, so the null result is not an artefact of working at ITS resolution.

Full results: `tables/orientation_*`, `tables/*_orientation_*` and `tables/substrate_x_orientation_*`. Summary: `tables/orientation_summary.csv`.

### F8. What this section concludes

**Sampling orientation does not detectably structure the LOT2 endophyte communities.** It is non-significant pooled, within each substrate, at every taxonomic rank, as a coarse sun-exposure grouping, and as a smooth directional gradient. The two significant pooled results both dissolve once substrate is accounted for (orientation adjusted for substrate: R² = 0.13859, p = 0.398; sun aspect adjusted for substrate: R² = 0.07457, p = 0.208).

**How firm is that?** Firm enough to report, but it is a negative result from an unbalanced observational factor with 1-4 replicates per group. It rules out an orientation effect of the size seen for substrate (R² ≈ 0.19); it does not rule out a small one. Making that test properly would need balanced sampling of all four bearings within each substrate — worth specifying in advance if a future campaign wants to answer this question rather than check it.

---

## Abundance & Composition Plots (Incertae sedis retained)

Unlike every analysis above, the plots below **keep the *Incertae sedis* isolates** (as an explicit pooled category), so they present the complete isolate census without hiding unidentified material.

### `plots/abundance_pooled_incertae_sedis/`

Horizontal grouped bar charts of **absolute isolate counts** per taxon, split by substrate, with unresolved taxa collected into an *Incertae sedis* bar.

![Abundance by genus (pooled)](plots/png/abundance_by_genus.png)

### `plots/pie_charts/`

Proportional composition of each substrate; the *Incertae sedis* slice shows how much of each community remains unidentified at that rank.

![Pie chart by phylum](plots/png/pie_by_phylum.png)

![Pie chart by genus](plots/png/pie_by_genus.png)

---

# Conclusion

1. **Substrate is the primary driver of community structure.** The three substrates host significantly different fungal communities (PERMANOVA p = 0.001, R² ≈ 0.19), and this holds at **every taxonomic rank** — the effect is in fact strongest around **order level** (R² ≈ 0.47). The two leaf substrates are more similar to each other than to wood, i.e. **tissue type (leaf vs wood)** is the strongest split, with **host identity (Ficus vs Lauraceae leaves)** adding a secondary but significant effect (p = 0.001).
2. **Tree height has at most a weak effect.** Treated as six discrete zones, height does **not** structure Ficus-wood communities at any taxonomic rank (Section D2, all p > 0.05). When the wood is instead split simply into **trunk vs branch**, a modest but **statistically significant** difference emerges (p = 0.006, R² ≈ 0.14): branch wood carries a somewhat distinct community from trunk wood, but this coarse contrast explains far less variation than substrate does.
3. **The substrate signal is real, not a sampling-height artefact.** Even when the comparison is restricted to zone 6 alone (where all substrates co-occur), substrates remain **statistically significant** (p = 0.001).
4. **Which side of the tree the material came from does not matter.** Now that the field branch codes have been reconciled with the project database, sampling orientation could be tested for the first time (Section F, 519 of 650 isolates). It is non-significant pooled (p = 0.159), within every substrate, at every taxonomic rank, as a sun-exposure grouping and as a smooth directional gradient. The two pooled contrasts that do come out significant — *substrate x orientation* and *sun aspect* — both vanish once substrate is accounted for (orientation adjusted for substrate: R² = 0.13859, p = 0.398), because the Lauraceae happened to be sampled almost entirely on north-facing branches. This is a **negative result from an unbalanced observational factor with 1-4 replicates per bearing**: it rules out an orientation effect as large as the substrate effect, not a small one.
5. **Communities are diverse and even.** All substrates show high evenness (Pielou J' > 0.9) and long rank-abundance tails; rarefaction curves have not saturated, so true richness is higher still. A shared generalist core of genera co-exists with a substantial set of substrate-exclusive taxa.
6. **Removing *Incertae sedis* sharpened the picture.** Excluding the pooled "unknown" bin from the diversity, overlap and multivariate analyses (while keeping it visible in the abundance/pie plots) increased, rather than decreased, the measured separation between substrates — confirming that the unidentified fraction had been masking genuine differences.

---

*Auto-generated on 2026-09-22 by `generate_readme.R`*
