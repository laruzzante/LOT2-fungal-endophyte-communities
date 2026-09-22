# LOT2 — Review of Analytical Methodology

**Scope.** An independent check of the LOT2 pipeline (`main.R` and the scripts it sources) against its own outputs, asking whether the strategy, parameters and products are fit for the stated goal: characterising fungal endophyte diversity in this ecosystem. Every figure quoted below was recomputed from the current input files.

**Inputs reviewed.** `LOT2_samples.xlsx` (650 isolates), `LOT2_pooled_counts.xlsx` (285 rows), and all 200+ generated tables and plots.

---

## Verdict

The pipeline is well engineered and the multivariate core is sound: the **substrate effect is real** and survives every reasonable challenge I applied, including adjustment for sampling depth. Reproducibility is good — one command regenerates everything.

But three results currently in the report should not be published as they stand: the **pie and abundance charts are contaminated by a spreadsheet totals row**, the **indicator species do not survive multiple-testing correction**, and the **host comparison is pseudoreplicated**. Separately, the study's headline question — how diverse is this system — is answered with descriptive statistics that omit the standard estimators, even though the data support them.

The single largest missed opportunity: all 650 ITS sequences are sitting unused in the workbook.

---

## 1. What works well

- **Substrate is a genuine signal.** PERMANOVA R² = 0.19 (p = 0.001) at ITS level, 0.40 at genus level. It holds after adjusting for sampling depth (marginal R² = 0.27, p = 0.001), and within zone 6 alone, so it is not an artefact of sampling height.
- **Rank-by-rank repetition** (Section D) is a genuinely good design choice. It is what revealed that the ITS-level ordinations are unusable while the genus-level ones are fine.
- **Incertae sedis handled per rank** rather than globally — correct, and better than the common practice of dropping unresolved isolates entirely.
- **Orientation conclusion is correctly argued.** The negative result is properly defended with marginal tests that dissolve the two confounded positives.

---

## 2. Problems, ranked

### P1 — CRITICAL: a totals row is plotted as a taxon

Row 286 of `LOT2_pooled_counts.xlsx` has no culture code and no taxon name, but carries the column sums: **219 / 264 / 167**. `parse_LOT2.R` reads it as a genotype. Downstream, `plot_pie.R` and `plot_abundance_pooled.R` bin it as *Incertae sedis*.

| Chart | *Incertae sedis* as plotted | Actual |
|:-------------|:------------------------|:-------------------|
| Pie by phylum | 50% / 50% / 51% | **0% / 0% / 2%** |
| Pie by genus | 56% / 60% / 60% | **11% / 20% / 20%** |

The README then interprets this slice as "how much of each community remains unidentified at that rank" — off by a factor of 25 at phylum level. It also explains why the two input files disagree by exactly 2× in total (438/528/334 vs 219/264/167).

The multivariate analyses are **unaffected** (they read the sample-level file), and the diversity indices escaped because the row has no taxonomy and is filtered with the *Incertae sedis* bin. Only the pie and abundance charts are wrong.

**Fix:** drop rows with no culture code and no taxon name at parse time, and assert that pooled totals equal sample-level totals.

### P2 — CRITICAL: no indicator species survives correction

`multipatt` p-values are filtered at raw p ≤ 0.05 with no correction, across 238–283 taxa per analysis. Recomputing with Benjamini–Hochberg:

| Analysis | Taxa tested | Raw p ≤ 0.05 | **BH q ≤ 0.05** |
|:-----------|:-----------:|:------------:|:---------------:|
| Substrate | 283 | 21 | **0** |
| Sun aspect | 238 | 4 | **0** |
| Orientation | 242 | 3 | **0** |

At 283 tests, ~14 hits at p ≤ 0.05 are expected by chance alone. The 21 reported for substrate are not distinguishable from noise. Across all analyses the pipeline currently reports **142 "significant" indicator taxa**; none would survive.

**Fix:** apply BH-FDR within each analysis and report q-values. Expect to report zero indicators — which is an honest and interesting result given the singleton-dominated community, not a failure.

### P3 — MAJOR: the host comparison is pseudoreplicated

As recorded, there is **no tree-level replication**. All 264 Ficus leaf isolates share one collection code (`LOT02-67508`) and the sheet's `Tree` column reads "Ficus" throughout. No tree identifier exists anywhere in the data.

- *Ficus* leaves vs *Ficus* wood is a **within-tree tissue contrast** — valid, but describing one tree.
- *Ficus* leaves vs Lauraceae leaves is a **host contrast with n = 1 tree per host**. PERMANOVA treats 20 sampling units as independent replicates of "host", so the p-value tests whether these two *individual trees* differ, not whether the two host species do.

This does not invalidate the numbers; it bounds what they mean. The README's conclusion that "host identity adds a secondary but significant effect" is not supported at the species level.

**Fix:** state the inferential unit explicitly, and rephrase host claims as descriptive of the sampled individuals. If other trees were sampled, add the tree ID and use it as a stratum.

### P4 — MAJOR: alpha diversity is computed with n = 1 per substrate

Section A derives all indices from the pooled file, giving **one value per substrate and no error, no test**. Statements like "the richest community is Ficus leaves" therefore have no statistical support.

The sample-level file already supports a proper test. Per sampling unit:

| Substrate | Units | Mean H′ | SD |
|:----------------|:-----:|:-------:|:----:|
| Ficus leaves | 10 | 2.92 | 0.16 |
| Ficus wood | 13 | 1.99 | 0.49 |
| Lauraceae leaves | 10 | 2.46 | 1.05 |

Kruskal–Wallis: **Shannon p = 0.00049, richness p = 0.00037**. The result the report wants is available and significant — it just is not being computed.

**Fix:** compute alpha diversity per sampling unit and test across substrates.

### P5 — MAJOR: sampling completeness is never quantified

The README says rarefaction curves "have not saturated" but never estimates how far off they are. They are far off:

| Substrate | Observed S | Chao1 (± SE) | Recovered | Good's coverage | Singletons |
|:----------------|:----------:|:------------:|:---------:|:---------------:|:----------:|
| Ficus leaves | 135 | 302 ± 50 | **45%** | 65% | 92 |
| Ficus wood | 82 | 206 ± 49 | **40%** | 67% | 55 |
| Lauraceae leaves | 105 | 186 ± 28 | **56%** | 71% | 63 |

Less than half the estimated richness was recovered, and about a third of isolates belong to genotypes seen exactly once. For a study whose goal is characterising diversity, these are the headline numbers and they are absent.

**Fix:** report Chao1/ACE with SEs and Good's coverage; compare diversity at equal *coverage* rather than equal *count*.

### P6 — MAJOR: sampling depth is confounded with substrate

Ficus wood units yielded roughly half as many isolates as leaf units (9.6 vs 20.2 and 18.8 per unit at genus level). Depth is itself a significant predictor of composition:

- community ~ depth: **R² = 0.216, p = 0.001** (genus level)
- marginal: depth R² = 0.069 (p = 0.003), substrate R² = 0.266 (p = 0.001)

Substrate survives, which is the important point — but part of the raw substrate effect is a depth effect, and this is neither reported nor adjusted for.

**Fix:** include depth as a covariate in the reported model, or standardise by coverage before comparing.

### P7 — MAJOR: dispersion heterogeneity at the ranks the report relies on

At genus level — the ranks the report (correctly) recommends as the readable ones — **PERMDISP p = 0.004**, with Ficus wood far more internally variable (mean distance to centroid 0.538) than either leaf substrate (0.288, 0.322). PERMANOVA cannot distinguish "different centroids" from "different spread", so the genus-level effect size is partly dispersion. Section D's interpretation does not mention it, although the p-value is printed in the summary table.

**Fix:** report PERMDISP alongside every PERMANOVA in the narrative, and interpret wood as *more heterogeneous*, which is itself an ecologically meaningful finding.

### P8 — MODERATE: ITS-level ordinations are degenerate *(already addressed)*

59% of sampling-unit pairs share no genotype, so their Bray–Curtis dissimilarity is tied at exactly 1. NMDS's default weak-tie handling leaves those pairs unconstrained, producing a near-zero stress that is an illusion; refitted with strong ties the stress is **0.249** (unreliable) at ITS level but **0.163** at genus level. The pipeline now reports both and labels the unusable plots.

### P9 — MODERATE: genotype names are not clean OTUs

Taxon identity is a BLAST-derived **name string**, not a clustered unit. **86 of 284 labels (30%)** are ambiguous — e.g. *"Aspergillus brunneus or A. niveoglaucus etc (4 spp.)"*, *"Cladosporium subuliforme or C. xylophilum etc (26 spp.) gen 2"*. Some bundle dozens of species into one "taxon" while nearly identical sequences may sit under two labels.

The two input files also disagree on **18 names in each direction**, including plain typos (*guandongensis* / *guangdongensis*; *prosopidis* / *Prosopidis*; *"P cf.. sporulosa"*). The lineage join silently falls back to matching on the leading genus word for these.

**All 650 ITS sequences are present in sheet `Feuil2`** (306–980 bp, one per isolate, 100% coverage) and are never used.

**Fix:** cluster the sequences into OTUs at 97% (or use ASVs), and use those as the analytical unit. This is the highest-value single change available — it removes the naming problem, and will very likely reduce the singleton fraction that is driving P2, P5 and P8.

### P10 — MODERATE: zone analyses are largely untestable

Ficus wood zones have **1 sampling unit for four of six zones** (zone 1: 1, zone 2: 1, zone 3: 1, zone 4: 2, zone 5: 1, zone 6: 7). "No zone effect" here means "no test was possible", not "no effect exists". The same applies to the per-substrate orientation tests (1–4 units per bearing).

**Fix:** state the replication explicitly wherever a null result is reported; consider dropping the discrete-zone analysis in favour of the trunk-vs-branch contrast, which does have replication.

### P11 — MINOR: no global multiplicity control

The pipeline generates roughly 19 PERMANOVA, 17 ANOSIM and 17 PERMDISP tests plus 120 pairwise rows. Holm correction is applied *within* each pairwise family only. At α = 0.05 across ~50 top-level tests, 2–3 false positives are expected.

**Fix:** designate a small set of pre-specified primary hypotheses; label everything else exploratory.

---

## 3. Proposed methodology

| Step | Current | Proposed |
|:--------------------|:------------------------|:------------------------|
| Analytical unit | BLAST name string | **97% OTUs clustered from `Feuil2` sequences** (VSEARCH/DADA2) |
| Richness | Observed S only | **Chao1/ACE + Good's coverage**, reported with SEs |
| Diversity comparison | Indices on pooled data, n = 1 | **Hill numbers (q = 0,1,2) per sampling unit**, coverage-standardised (`iNEXT`), tested across substrates |
| Rarefaction | Individual-based on pooled group | **Sample-based rarefaction + extrapolation** with CIs |
| Dissimilarity | Bray–Curtis on raw counts | Bray–Curtis on coverage-standardised data; add **Raup–Crick** (null-model based), which is designed for exactly this sparsity |
| Ordination | NMDS at ITS level | **NMDS at genus/family level**, or PCoA; report tie-aware stress |
| Model | `community ~ factor` | `community ~ depth + factor`, with **tree/zone as stratum** where replication allows |
| Indicators | raw p ≤ 0.05 | **BH-FDR q ≤ 0.05** |
| Dispersion | computed, not discussed | reported next to every PERMANOVA |

**Priority order:** P1 and P2 before anything is shown to anyone (both are wrong-as-published). Then P9 (OTU clustering), which improves P5 and P8 at the same time. Then P4 and P5, which are what the diversity question actually asks for. P3 is a framing fix, not an analysis fix, but it changes what the paper can claim.

---

## 4. Questions for the data owner

1. **How many trees were sampled per host?** If more than one, the tree ID needs adding — it determines whether the host comparison is publishable as a species-level result.
2. **Is "Lauraceae" leaves or wood?** `LOT2_pooled_counts.xlsx` labels the column "Lauraceae leaves", but the samples workbook sheet is titled *"66. Fungi-Endo wood (Host)"*, and the 219 Lauraceae isolates match **none** of the three per-substrate sheets. The report's central claim — that leaf-vs-wood tissue type is the strongest split — depends on this.
3. **Is the Ficus growing on the Lauraceae?** "Host" in the sheet title suggests a strangler relationship. If so the two trees are not independent, which matters for P3.
4. **What do the per-substrate sheets represent** (264 / 167 / 20 isolates) relative to the 650 analysed? Sheet 66 in particular contains only 20 isolates that appear nowhere in the analysis.
5. **Should isolates be de-replicated per sampling unit?** Several units contain the same genotype repeatedly; whether these are independent colonisations or one colony subsampled changes the abundance weighting throughout.

---

## 5. Resolution status

All eleven problems have been addressed in the pipeline, and the five provenance questions were answered by the data owner. Their answers are now recorded in the README ("Sampling design: what was sampled, and the assumptions behind it") so collaborators can check the reasoning.

| # | Problem | Status |
|:---|:------------------------|:------------------------|
| P1 | Totals row plotted as a taxon | **Fixed** — `parse_LOT2.R` drops rows with neither culture code nor taxon name; pie/abundance charts now show the true 0–20% *Incertae sedis*, and an integrity check asserts pooled = sample-level totals |
| P2 | Indicator species uncorrected | **Fixed** — Benjamini–Hochberg applied; full p/q tables written to `*_indicator_species_all.csv` |
| P3 | Pseudoreplication | **Documented** — confirmed 1 tree per host, physically interlocked; `tree_id` added, scope caveat now heads the Conclusion |
| P4 | Alpha diversity untestable | **Fixed** — `diversity_analysis.R` computes per-unit diversity with Kruskal–Wallis and pairwise Wilcoxon |
| P5 | Completeness unquantified | **Fixed** — Chao1/ACE, Good's coverage and equal-coverage richness in new Section A0 |
| P6 | Sampling depth confounded | **Fixed** — every PERMANOVA also reported adjusted for depth (marginal/type-III) |
| P7 | Dispersion heterogeneity | **Fixed** — PERMDISP reported alongside PERMANOVA, and the wood-heterogeneity result is now interpreted rather than buried |
| P8 | NMDS degeneracy | **Fixed** — tie-aware stress reported; OTU clustering removed the degeneracy at the primary level |
| P9 | Name strings used as OTUs | **Fixed** — `cluster_otus.sh` clusters the 650 ITS sequences; 97% OTUs are now the analytical unit |
| P10 | Zone analyses untestable | **Documented** — replication and an explicit caution printed in every PERMANOVA report |
| P11 | No global multiplicity control | **Partly fixed** — FDR within IndVal, Holm across alpha indices and pairwise tests. Designating primary vs exploratory hypotheses remains an authorial decision |

### What the fixes changed in the results

- **OTU clustering rescued real signal.** On the raw name labels, *no* indicator taxon survived FDR correction. On 97% OTUs, **14 substrate indicators survive** — fewer tests and consolidated singletons recovered a result the old pipeline could not have reported honestly.
- **Alpha diversity is now testable and significant.** Richness and Shannon differ across substrates (Holm-adjusted p = 0.0015), while **evenness does not** (p = 0.27).
- **An apparent contradiction turned out to be a finding.** Per sampling unit, *Ficus* wood is the *least* diverse substrate; at equal coverage it is the *most* diverse. Leaf units are individually richer, but the wood community is more heterogeneous between units and has a longer tail of rare taxa — the same heterogeneity PERMDISP detects.
- **The substrate effect held up.** It survives adjustment for sampling depth, so it is not an artefact of unequal recovery.
- **One analysis was deleted, not fixed.** `lauraceae_leaves_trunk_vs_branch_*` compared branch against branch, because the Lauraceae trunk was never sampled.

### Still open

**Whether the Lauraceae material is leaves or branch wood.** This is the one unresolved question, and the report's tissue-type conclusion depends on it. If it is branch wood, "the two leaf substrates resemble each other more than either resembles wood" cannot stand. Everything else is unaffected.
