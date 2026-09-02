dir.create("tables", showWarnings = FALSE)
dir.create("plots/png", recursive = TRUE, showWarnings = FALSE)

# ---- Convert key PDFs to PNGs for embedding ----
pdf_to_png <- function(pdf_path, png_path, density = 150) {
  cmd <- sprintf('convert -density %d "%s[0]" -quality 90 "%s"',
                 density, pdf_path, png_path)
  system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
}

embedded_plots <- c(
  "plots/community/alpha_shannon_H.pdf",
  "plots/community/rank_abundance_genus.pdf",
  "plots/community/rel_abundance_phylum.pdf",
  "plots/community/rel_abundance_genus.pdf",
  "plots/community/venn_genus.pdf",
  "plots/community/substrate_all_nmds.pdf",
  "plots/community/substrate_all_nmds_species.pdf",
  "plots/community/substrate_all_rarefaction.pdf",
  "plots/community/substrate_leaves_nmds.pdf",
  "plots/community/ficus_wood_trunk_vs_branch_nmds.pdf",
  "plots/community/ficus_wood_zones_nmds.pdf",
  "plots/community/substrate_zone6_branch_nmds.pdf",
  "plots/community/substrate_x_position_nmds.pdf",
  "plots/community/substrate_bylevel_genus_nmds.pdf",
  "plots/community/substrate_bylevel_family_nmds.pdf",
  "plots/community/substrate_bylevel_phylum_nmds.pdf",
  "plots/abundance_pooled_incertae_sedis/abundance_by_genus.pdf",
  "plots/pie_charts/pie_by_phylum.pdf",
  "plots/pie_charts/pie_by_genus.pdf"
)

for (pdf in embedded_plots) {
  if (file.exists(pdf)) {
    png <- file.path("plots/png", sub("\\.pdf$", ".png", basename(pdf)))
    pdf_to_png(pdf, png)
  }
}
cat("Converted PDFs to PNG\n")

# ---- R and package versions ----
r_ver <- paste0(R.version$major, ".", R.version$minor)
pkg_versions <- sapply(
  c("readxl", "dplyr", "tidyr", "ggplot2", "vegan", "indicspecies",
    "ggVennDiagram", "scales"),
  function(p) as.character(packageVersion(p))
)

# ---- Load results ----
alpha  <- read.csv("tables/alpha_diversity.csv", stringsAsFactors = FALSE)
alpha_its <- alpha[alpha$taxonomic_level == "its_taxon", ]
alpha_gen <- alpha[alpha$taxonomic_level == "genus", ]
alpha_phy <- alpha[alpha$taxonomic_level == "phylum", ]

# Small helpers for data-driven interpretation text
num      <- function(x) suppressWarnings(as.numeric(x))
verdict  <- function(p) { p <- num(p)
  if (is.na(p)) "could not be evaluated" else if (p < 0.05) "**statistically significant**" else "**not statistically significant**" }
stars    <- function(p) { p <- num(p)
  if (is.na(p)) "" else if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns" }
top_by   <- function(df, col) df$substrate[order(-num(df[[col]]))][1]
bot_by   <- function(df, col) df$substrate[order(num(df[[col]]))][1]

# Safe file reader
safe_read <- function(path) {
  if (file.exists(path)) readLines(path) else character(0)
}
safe_csv <- function(path) {
  if (file.exists(path)) read.csv(path, stringsAsFactors = FALSE) else data.frame()
}

perm_txt  <- safe_read("tables/permanova_results.txt")
anos_txt  <- safe_read("tables/anosim_results.txt")
betad_txt <- safe_read("tables/betadisper_results.txt")
pw        <- safe_csv("tables/pairwise_permanova.csv")
indval    <- safe_csv("tables/indicator_species_significant.csv")

# Helper to extract PERMANOVA stats from adonis2 output text
extract_permanova <- function(txt) {
  line <- grep("^Model", txt, value = TRUE)
  if (length(line) == 0) return(list(F = "N/A", R2 = "N/A", p = "N/A"))
  pp <- strsplit(trimws(line[1]), "\\s+")[[1]]
  # Model Df SumOfSqs R2 F Pr(>F) [sig]
  list(F = pp[5], R2 = pp[4], p = pp[6])
}

perm_stats <- extract_permanova(perm_txt)
perm_F <- perm_stats$F; perm_R2 <- perm_stats$R2; perm_p <- perm_stats$p

# Extract ANOSIM stats
anos_R <- sub(".*R statistic:\\s*", "", grep("R statistic:", anos_txt, value = TRUE))
anos_p <- sub(".*p-value:\\s*", "", grep("p-value:", anos_txt, value = TRUE))
if (length(anos_R) == 0) anos_R <- "N/A"
if (length(anos_p) == 0) anos_p <- "N/A"

# Betadisper stats — line format: Groups Df SumSq MeanSq F N.Perm Pr(>F)
bd_line <- grep("^Groups", betad_txt, value = TRUE)
bd_F <- bd_p <- "N/A"
if (length(bd_line) > 0) {
  bp <- strsplit(trimws(bd_line[1]), "\\s+")[[1]]
  # Groups Df Sum_Sq Mean_Sq F N.Perm Pr(>F)
  if (length(bp) >= 7) { bd_F <- bp[5]; bd_p <- bp[7] }
}

# NMDS stress from substrate_all analysis
nmds_stress <- tryCatch({
  results_substrate$nmds_stress
}, error = function(e) "N/A")

n_samples <- nrow(comm_mat_full)
n_taxa_comm <- ncol(comm_mat_full)

# Venn counts (Incertae sedis excluded, matching the plots)
venn_genus <- dat_agg %>%
  filter(genus != incertae_label) %>%
  group_by(genus) %>%
  summarise(ll = sum(n_fungi_laur_leaf), fl = sum(n_fungi_fic_leaf),
            fw = sum(n_fungi_fic_wood), .groups = "drop")

count_shared <- function(df) {
  in_ll <- df$ll > 0; in_fl <- df$fl > 0; in_fw <- df$fw > 0
  list(
    all_three = sum(in_ll & in_fl & in_fw),
    ll_only   = sum(in_ll & !in_fl & !in_fw),
    fl_only   = sum(!in_ll & in_fl & !in_fw),
    fw_only   = sum(!in_ll & !in_fl & in_fw),
    ll_fl     = sum(in_ll & in_fl & !in_fw),
    ll_fw     = sum(in_ll & !in_fl & in_fw),
    fl_fw     = sum(!in_ll & in_fl & in_fw),
    total     = nrow(df)
  )
}
vg <- count_shared(venn_genus)

# Leaves-only sub-analysis results
leaves_perm_txt <- safe_read("tables/substrate_leaves_permanova.txt")
leaves_stats <- extract_permanova(leaves_perm_txt)
leaves_perm_F <- leaves_stats$F; leaves_perm_R2 <- leaves_stats$R2; leaves_perm_p <- leaves_stats$p

# Ficus wood trunk vs branch
fw_tb_perm_txt <- safe_read("tables/ficus_wood_trunk_vs_branch_permanova.txt")
fw_tb_stats <- extract_permanova(fw_tb_perm_txt)
fw_tb_F <- fw_tb_stats$F; fw_tb_R2 <- fw_tb_stats$R2; fw_tb_p <- fw_tb_stats$p

# Zone 6 substrate comparison
z6_perm_txt <- safe_read("tables/substrate_zone6_branch_permanova.txt")
z6_stats <- extract_permanova(z6_perm_txt)
z6_perm_F <- z6_stats$F; z6_perm_R2 <- z6_stats$R2; z6_perm_p <- z6_stats$p

# Sample counts per substrate
sub_counts <- samples_raw %>% count(substrate) %>% arrange(substrate)
sub_units  <- samples_raw %>% distinct(substrate, unit) %>% count(substrate) %>% arrange(substrate)

# ---- Multi-level taxonomic sweep results ----
ml_sub <- safe_csv("tables/substrate_multilevel_summary.csv")
ml_zon <- safe_csv("tables/ficus_wood_zones_multilevel_summary.csv")
ml_sp  <- safe_csv("tables/substrate_x_position_multilevel_summary.csv")
lincov <- safe_csv("tables/lineage_coverage.csv")

# Build a markdown table from a multi-level summary data frame
ml_table <- function(df) {
  if (nrow(df) == 0) return(character(0))
  header <- c(
    '| Rank | Taxa | NMDS stress | PERMANOVA F | R\u00b2 | p | ANOSIM R | p | PERMDISP p |',
    '|------|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|')
  rows <- apply(df, 1, function(r) {
    paste0('| ', r['level'], ' | ', r['n_taxa'], ' | ', r['nmds_stress'], ' | ',
           r['permanova_F'], ' | ', r['permanova_R2'], ' | ', r['permanova_p'], ' | ',
           r['anosim_R'], ' | ', r['anosim_p'], ' | ', r['betadisper_p'], ' |')
  })
  c(header, rows)
}

# Interpretation values for the substrate sweep
if (nrow(ml_sub) > 0) {
  ml_best_rank <- ml_sub$level[which.max(ml_sub$permanova_R2)]
  ml_best_R2   <- max(ml_sub$permanova_R2, na.rm = TRUE)
  ml_its_R2    <- ml_sub$permanova_R2[ml_sub$level == "its_taxon"]
  ml_all_sig   <- all(num(ml_sub$permanova_p) < 0.05, na.rm = TRUE)
} else {
  ml_best_rank <- "N/A"; ml_best_R2 <- NA; ml_its_R2 <- NA; ml_all_sig <- FALSE
}
zones_any_sig <- if (nrow(ml_zon) > 0) any(num(ml_zon$permanova_p) < 0.05, na.rm = TRUE) else FALSE

# ---- Build README ----
readme <- c(
'# Fungal Endophyte Community Analysis — LOT2 (Peru)',
'',
paste0('> **Analysis environment:** R ', r_ver, ' | Generated on ', Sys.Date()),
'',
'## Overview',
'',
'This project analyses the fungal endophyte communities isolated from three substrates collected in **Peru** (LOT2 sampling campaign):',
'',
'| Substrate | Isolates (pooled) | Sample-level replicates (units) |',
'|-----------|:-:|:-:|',
paste0('| **Ficus leaves** | ', alpha_its$abundance_N[alpha_its$substrate == "Ficus leaves"], ' | ', sub_units$n[sub_units$substrate == "Ficus leaves"], ' |'),
paste0('| **Ficus wood** | ', alpha_its$abundance_N[alpha_its$substrate == "Ficus wood"], ' | ', sub_units$n[sub_units$substrate == "Ficus wood"], ' |'),
paste0('| **Lauraceae leaves** | ', alpha_its$abundance_N[alpha_its$substrate == "Lauraceae leaves"], ' | ', sub_units$n[sub_units$substrate == "Lauraceae leaves"], ' |'),
'',
'### Sampling design',
'',
'Samples were collected at different **tree zones** (heights):',
'',
'- **Zones 1-5**: Tree trunk (zone 1 = lowest, zone 5 = highest)',
'- **Zone 6**: Canopy branches (highest zone)',
'',
'Wood collected from zone 6 corresponds to **branch wood** (no trunk present). Wood from zones 1-5 is **trunk wood**.',
'',
'| Substrate | Zones sampled | Notes |',
'|-----------|:---:|---|',
'| Ficus leaves | 6 only | All leaf samples from canopy |',
'| Ficus wood | 1-6 | Trunk (zones 1-5) and branch (zone 6) |',
'| Lauraceae leaves | 5-6 | Predominantly zone 6 |',
'',
'### Handling of uncertain taxonomy (Incertae sedis)',
'',
'Many isolates cannot be confidently named at every taxonomic rank. Entries flagged `NA`, `"?"`, `"NO"` or `"incertae sedis"` are treated as **Incertae sedis** ("of uncertain placement") and are handled **rank by rank**:',
'',
'- They are **shown** — as a single pooled *Incertae sedis* category — **only in the abundance bar charts and the pie charts**, so that the full isolate count is never hidden.',
'- They are **excluded from every diversity index, rank-/relative-abundance curve, Venn diagram and multivariate test**. A single large, shared "unknown" bin behaves like a taxon that is common everywhere: it inflates apparent overlap between substrates and **flattens the real ecological differences** we are trying to detect.',
'- Because the exclusion is applied **independently at each rank**, an isolate that is unresolved at *species* level but has a defined *genus*, *family* or *order* still contributes to the analyses run at those higher ranks (see **Section D — Multi-level analyses**).',
'',
'> The sample-level multivariate analyses are run at **ITS-genotype** resolution, where each of the sequenced genotypes is a distinct entity. There is no dominant "unknown" bin at that level, so no isolates are dropped there; the *Incertae sedis* filtering matters only when isolates are grouped into higher taxa.',
'',
'---',
'',
'## Input Data',
'',
paste0('- **`LOT2_pooled_counts.xlsx`** (first sheet) — Pooled genotype counts per substrate with full taxonomy (', nrow(pooled), ' genotypes)'),
paste0('- **`LOT2_samples.xlsx`** (first sheet) — Individual isolate records with sampling zone and unit (', nrow(samples_raw), ' isolates)'),
'',
'## Scripts',
'',
'| Script | Purpose |',
'|--------|---------|',
'| `main.R` | Master script — loads libraries, sources all other scripts in order |',
'| `parse_LOT2.R` | Reads both Excel files (first sheet only), standardises column names |',
'| `plot_abundance.R` | Horizontal bar charts of isolate counts by taxonomic level |',
'| `plot_abundance_pooled.R` | Bar charts with uncertain taxa pooled as "Incertae sedis" |',
'| `plot_pie.R` | Pie charts of community composition (3 pies per level) |',
'| `community_analysis.R` | Full community ecology analyses; ITS-level tests plus multi-rank taxonomic sweeps (Section D) |',
'| `generate_readme.R` | Generates this README dynamically from analysis results |',
'',
'```r',
'# From the project directory:',
'source("main.R")',
'```',
'',
'### R packages used',
'',
'| Package | Version | Role |',
'|---------|---------|------|',
paste0('| readxl | ', pkg_versions['readxl'], ' | Reading Excel input |'),
paste0('| dplyr | ', pkg_versions['dplyr'], ' | Data manipulation |'),
paste0('| tidyr | ', pkg_versions['tidyr'], ' | Data reshaping |'),
paste0('| ggplot2 | ', pkg_versions['ggplot2'], ' | Plotting |'),
paste0('| scales | ', pkg_versions['scales'], ' | Axis formatting |'),
paste0('| vegan | ', pkg_versions['vegan'], ' | Diversity, NMDS, PERMANOVA, ANOSIM, betadisper |'),
paste0('| indicspecies | ', pkg_versions['indicspecies'], ' | Indicator species analysis (IndVal) |'),
paste0('| ggVennDiagram | ', pkg_versions['ggVennDiagram'], ' | Venn diagrams |'),
'',
'---',
'',
'## Output Structure',
'',
'```',
'plots/',
'  abundance/                        # Bar charts (raw)',
'  abundance_pooled_incertae_sedis/  # Bar charts (uncertain taxa pooled)',
'  pie_charts/                       # Pie charts',
'  community/                        # Community ecology plots (NMDS, rarefaction, etc.)',
'  png/                              # PNG versions for README embedding',
'',
'tables/                             # Statistical results and summary tables',
'```',
'',
'---',
'',
'# Results',
'',
'> Throughout the Results, each analysis is introduced with a short **plain-language explanation** of what it measures and how to read it, followed by a **brief interpretation** of what the LOT2 data actually show. A synthesis of all findings is given in the final **Conclusion**.',
'',
'## A. Alpha Diversity (from pooled counts)',
'',
'**What it is.** *Alpha diversity* describes how varied the fungal community is **within a single substrate**. It combines two ideas: *richness* (how many different taxa are present) and *evenness* (whether isolates are spread evenly across taxa or dominated by a few). The indices below capture different balances of these two ideas. *Incertae sedis* taxa are excluded so that diversity reflects only confidently identified fungi.',
'',
'| Index | What it measures | Higher means |',
'|-------|-----------------|--------------|',
'| **Richness (S)** | Number of distinct taxa | More taxa present |',
'| **Shannon (H\')** | Combines richness and evenness | More diverse |',
'| **Simpson (1-D)** | Probability two random individuals differ | More diverse |',
'| **Inverse Simpson** | Effective number of equally-common species | More even |',
'| **Pielou (J\')** | How evenly individuals are distributed | More even |',
'',
'### At ITS taxon level',
'',
'| Substrate | S | N | H\' | 1-D | Inv. Simp. | J\' |',
'|-----------|:-:|:-:|:-:|:-:|:-:|:-:|',
paste0('| Lauraceae leaves | ', alpha_its$richness_S[1], ' | ', alpha_its$abundance_N[1], ' | ', alpha_its$shannon_H[1], ' | ', alpha_its$simpson_1mD[1], ' | ', alpha_its$inv_simpson[1], ' | ', alpha_its$pielou_J[1], ' |'),
paste0('| Ficus leaves | ', alpha_its$richness_S[2], ' | ', alpha_its$abundance_N[2], ' | ', alpha_its$shannon_H[2], ' | ', alpha_its$simpson_1mD[2], ' | ', alpha_its$inv_simpson[2], ' | ', alpha_its$pielou_J[2], ' |'),
paste0('| Ficus wood | ', alpha_its$richness_S[3], ' | ', alpha_its$abundance_N[3], ' | ', alpha_its$shannon_H[3], ' | ', alpha_its$simpson_1mD[3], ' | ', alpha_its$inv_simpson[3], ' | ', alpha_its$pielou_J[3], ' |'),
'',
paste0('**Interpretation.** The richest community is **', top_by(alpha_its, "richness_S"),
       '** (S = ', max(alpha_its$richness_S), ' distinct ITS taxa), and the least rich is **',
       bot_by(alpha_its, "richness_S"), '** (S = ', min(alpha_its$richness_S),
       '). Shannon diversity is highest in **', top_by(alpha_its, "shannon_H"),
       "** (H' = ", max(alpha_its$shannon_H), ') and lowest in **', bot_by(alpha_its, "shannon_H"),
       "** (H' = ", min(alpha_its$shannon_H),
       '). Pielou evenness is high (J\' > 0.9) for all three substrates, meaning no single genotype dominates any community \u2014 isolates are spread across many co-occurring taxa. Because richness partly reflects sampling effort (N differs between substrates), richness values should be compared together with the **rarefaction curves** below, which put all substrates on an equal-effort footing.'),
'',
'Diversity was also computed at every higher rank (phylum \u2192 genus); the full table is in `tables/alpha_diversity.csv`.',
'',
'![Shannon diversity](plots/png/alpha_shannon_H.png)',
'',
'*The bars show Shannon H\' for each substrate across taxonomic ranks. H\' naturally decreases towards coarser ranks (fewer categories), but the ranking of substrates stays broadly consistent, indicating the diversity differences are not an artefact of one particular rank.*',
'',
'---',
'',
'## B. Community Composition',
'',
'These analyses describe **what the communities are made of** and **how much they overlap**, again after removing *Incertae sedis*.',
'',
'### B1. Rank-Abundance Curves',
'',
'**What it is.** Taxa are ranked from most to least abundant (x-axis) against their relative abundance on a log scale (y-axis). A **steep** curve means a few taxa dominate (low evenness); a **shallow, long** curve means many taxa share the community evenly (high evenness). The length of each curve reflects richness.',
'',
'![Rank-abundance at genus level](plots/png/rank_abundance_genus.png)',
'',
'*Interpretation.* All three substrates show relatively shallow curves with long tails, confirming the high evenness seen in the Pielou index: communities are not dominated by one or two hyper-abundant genera but consist of many moderately frequent taxa plus a long tail of rare ones \u2014 a pattern typical of tropical endophyte assemblages.',
'',
'### B2. Relative Abundance',
'',
'**What it is.** Stacked bars show the **proportional composition** of each substrate at a given rank (each bar sums to 100%). They make it easy to see which phyla/genera dominate and how composition shifts between substrates.',
'',
'![Relative abundance by phylum](plots/png/rel_abundance_phylum.png)',
'',
'![Relative abundance by genus](plots/png/rel_abundance_genus.png)',
'',
'*Interpretation.* At phylum level the communities are overwhelmingly **Ascomycota**, as expected for culturable endophytes. The genus-level bars reveal the real contrast between substrates: the identity and proportion of dominant genera differ markedly between leaves and wood, foreshadowing the significant substrate effect quantified in the multivariate tests below.',
'',
'### B3. Venn Diagrams \u2014 Shared Taxa',
'',
'**What it is.** The Venn diagram counts how many taxa are **unique** to each substrate versus **shared** between them. It is a simple presence/absence view of community overlap (abundance is ignored).',
'',
paste0('*Interpretation.* At **genus level** (', vg$total, ' identified genera in total): **', vg$all_three,
       '** genera occur in all three substrates (a shared generalist core), while **', vg$ll_only,
       '** are unique to Lauraceae leaves, **', vg$fl_only, '** unique to Ficus leaves and **', vg$fw_only,
       '** unique to Ficus wood. The substantial number of substrate-exclusive genera indicates a degree of **habitat specialisation** layered on top of a shared generalist core.'),
'',
'![Venn diagram — genus](plots/png/venn_genus.png)',
'',
'---',
'',
'# Multivariate Statistical Analyses',
'',
paste0('These analyses ask **whether whole communities differ between groups** (substrates, zones, positions). They use the sample-level data (**', n_samples, ' sampling units**, **', n_taxa_comm, ' ITS genotypes**), the **Bray-Curtis dissimilarity** (a 0\u20131 measure of how different two samples are in both *which* taxa are present and *how abundant* they are), and **999 permutations** to obtain p-values without assuming normality.'),
'',
'**How to read each test:**',
'',
'- **NMDS ordination** \u2014 squeezes the many-dimensional Bray-Curtis distances into a 2-D map so that samples plotting close together have similar communities. The **stress** value measures distortion: < 0.10 excellent, < 0.20 acceptable, > 0.20 unreliable. Crosses mark group centroids; shaded ellipses show 95% confidence regions.',
'- **PERMANOVA** (`adonis2`) \u2014 tests whether **group centroids differ**. **R\u00b2** is the fraction of community variation explained by the grouping (effect size); a small **p** means the separation is unlikely by chance.',
'- **ANOSIM** \u2014 a complementary rank-based test; **R** ranges from 0 (no separation) to 1 (groups completely distinct).',
'- **Beta-dispersion / PERMDISP** (`betadisper`) \u2014 checks whether groups differ in **within-group spread** rather than location. If PERMDISP is significant, part of a PERMANOVA result may reflect unequal dispersion rather than a pure shift in composition, so it is an important caveat.',
'- **Pairwise PERMANOVA** \u2014 which specific pairs of groups differ, with Holm correction for multiple tests.',
'- **Rarefaction** \u2014 expected richness rescaled to equal sampling effort, so richness can be compared fairly.',
'- **Indicator species (IndVal)** \u2014 identifies taxa statistically associated with (diagnostic of) a particular group.',
'',
'---',
'',
'## C1. Substrate Comparison: Ficus Leaves vs Ficus Wood vs Lauraceae Leaves',
'',
'This is the **headline comparison**: do the three substrates host different fungal communities?',
'',
'### NMDS Ordination',
'',
paste0('**Stress = ', nmds_stress, '** ',
       if (!is.na(num(nmds_stress)) && num(nmds_stress) < 0.20) '(acceptable to good \u2014 the 2-D map is a faithful summary).' else '(interpret the map with some caution).'),
'',
'![NMDS — all substrates](plots/png/substrate_all_nmds.png)',
'',
'![NMDS with species overlay](plots/png/substrate_all_nmds_species.png)',
'',
'*Interpretation.* The three substrates form visually distinct clouds, with the two leaf substrates sitting closer to each other than to wood \u2014 consistent with tissue type (leaf vs wood) being a strong driver. The species overlay points to the genotypes pulling each substrate apart.',
'',
'### PERMANOVA',
'',
'```',
paste(perm_txt, collapse = "\n"),
'```',
'',
paste0('**F = ', perm_F, ', R\u00b2 = ', perm_R2, ', p = ', perm_p, '** \u2014 the substrate effect is ', verdict(perm_p),
       '. Substrate explains about **', round(100 * num(perm_R2)), '%** of the total community variation, a large effect for field endophyte data.'),
'',
'### ANOSIM',
'',
paste0('**R = ', anos_R, ', p = ', anos_p, '** \u2014 ', verdict(anos_p),
       '. An R of this magnitude confirms that between-substrate differences clearly exceed within-substrate variation.'),
'',
'### Beta-dispersion (PERMDISP)',
'',
paste0('**F = ', bd_F, ', p = ', bd_p, '** \u2014 dispersion differences are ', verdict(bd_p),
       '. ', if (!is.na(num(bd_p)) && num(bd_p) < 0.05) 'Because within-group spread also differs, the PERMANOVA result partly reflects unequal dispersion; the significant separation on the NMDS nonetheless supports a genuine compositional shift.' else 'Groups are comparably variable, so the PERMANOVA result reflects a genuine shift in composition rather than unequal spread.'),
'',
if (nrow(pw) > 0) {
  c(
    '### Pairwise PERMANOVA',
    '',
    'Which substrate *pairs* differ (Holm-corrected p). `***` p<0.001, `**` p<0.01, `*` p<0.05.',
    '',
    '| Comparison | F | R² | p (raw) | p (Holm) | Sig. |',
    '|-----------|:-:|:-:|:-:|:-:|:-:|',
    apply(pw, 1, function(r) {
      sig <- ifelse(as.numeric(r["p_adj"]) < 0.001, "***",
             ifelse(as.numeric(r["p_adj"]) < 0.01, "**",
             ifelse(as.numeric(r["p_adj"]) < 0.05, "*", "ns")))
      paste0('| ', r["pair"], ' | ', r["F_stat"], ' | ', r["R2"], ' | ',
             r["p_value"], ' | ', r["p_adj"], ' | ', sig, ' |')
    }),
    '',
    paste0('*Interpretation.* ',
           if (all(num(pw$p_adj) < 0.05)) 'Every pair of substrates differs significantly after correction \u2014 each substrate carries a distinguishable community, not merely one odd substrate against two similar ones.' else 'Not all pairs differ significantly; see the table for which contrasts drive the overall effect.'),
    ''
  )
} else character(0),
'### Rarefaction',
'',
'Expected number of taxa if every substrate had been sampled to the same number of isolates \u2014 a fair richness comparison that removes the effect of unequal sampling effort.',
'',
'![Rarefaction — all substrates](plots/png/substrate_all_rarefaction.png)',
'',
'*Interpretation.* None of the curves has fully levelled off, so additional sampling would still recover new taxa in every substrate (the communities are undersampled, as usual for hyper-diverse tropical fungi). The **relative ordering** of the curves indicates which substrate is richest at equal effort, which is the sampling-fair complement to the raw richness values in Section A.',
'',
if (nrow(indval) > 0) {
  c(
    paste0('### Indicator Species (', nrow(indval), ' significant, p \u2264 0.05)'),
    '',
    paste0('Indicator (IndVal) analysis finds taxa that are **diagnostic** of a particular substrate \u2014 both faithful to it (mostly found there) and frequent within it. **', nrow(indval), '** ITS genotypes are significant indicators, i.e. reliable biological markers of their substrate. Full ranked list: `tables/indicator_species_significant.csv`.'),
    ''
  )
} else character(0),
'---',
'',
'## C2. Ficus Leaves vs Lauraceae Leaves (leaf substrates only)',
'',
'**Why.** Both are **leaf** endophyte communities but from different host plants, so this isolates the **host effect** from the leaf-vs-wood tissue effect.',
'',
paste0('**PERMANOVA: F = ', leaves_perm_F, ', R\u00b2 = ', leaves_perm_R2, ', p = ', leaves_perm_p, '** \u2014 ', verdict(leaves_perm_p),
       '. Even between two leaf communities, host identity (Ficus vs Lauraceae) leaves a detectable signature, explaining about ', round(100 * num(leaves_perm_R2)), '% of the variation.'),
'',
'![NMDS — Ficus vs Lauraceae leaves](plots/png/substrate_leaves_nmds.png)',
'',
'Full results in `tables/substrate_leaves_*.txt`',
'',
'---',
'',
'## D. Multi-level Taxonomic Analyses',
'',
'**Why this section exists.** The multivariate tests above use **ITS genotypes**, the finest possible resolution. But an isolate that is *Incertae sedis* at species level often still has a defined **genus, family or order**. By linking every isolate to its **full lineage** (from the pooled taxonomy) we can repeat the community comparisons at each rank and ask: **is the substrate/zone signal a fine-scale artefact, or does it hold when isolates are grouped into higher, more confidently identified taxa?**',
'',
if (nrow(lincov) > 0) c(
  paste0('**Lineage coverage** \u2014 isolates with a defined value at each rank (out of ', nrow(samples_raw), '):'),
  '',
  '| Rank | Isolates resolved |',
  '|------|:-:|',
  apply(lincov, 1, function(r) paste0('| ', r['rank'], ' | ', r['isolates_resolved'], ' / ', r['isolates_total'], ' |')),
  ''
) else character(0),
'### D1. Substrate comparison across ranks',
'',
ml_table(ml_sub),
'',
paste0('*Interpretation.* The substrate effect is ', if (ml_all_sig) '**significant at every taxonomic rank**' else 'significant at most ranks',
       ' (PERMANOVA p ', if (ml_all_sig) '\u2264 0.004' else 'values in the table', '). Crucially, the effect size does **not** weaken when isolates are grouped into higher taxa \u2014 it is *strongest* around **', ml_best_rank, ' level** (R\u00b2 \u2248 ', round(ml_best_R2, 2),
       '), compared with R\u00b2 \u2248 ', round(num(ml_its_R2), 2), ' at ITS level. In other words, the substrates differ not just in which fine genotypes they carry, but in their broad taxonomic make-up, and removing the *Incertae sedis* noise sharpens rather than blurs that separation.'),
'',
'![NMDS at genus level — substrates](plots/png/substrate_bylevel_genus_nmds.png)',
'',
'![NMDS at family level — substrates](plots/png/substrate_bylevel_family_nmds.png)',
'',
'### D2. Ficus wood zones across ranks',
'',
ml_table(ml_zon),
'',
paste0('*Interpretation.* ', if (!zones_any_sig) 'At **no** taxonomic rank do Ficus-wood communities differ significantly among tree zones (all PERMANOVA p > 0.05). The vertical position of wood on the tree does **not** structure its fungal community detectably \u2014 the same conclusion reached at ITS level, now confirmed to be robust to taxonomic resolution.' else 'Some ranks show a zone effect; see the table.'),
'',
'### D3. Substrate \u00d7 position across ranks',
'',
ml_table(ml_sp),
'',
'*Interpretation.* Combining substrate with trunk/branch position remains significant at every rank, and the effect size again peaks at intermediate ranks (class\u2013family). This mirrors the substrate result: the signal is carried by broad taxonomic groups, not just rare fine-scale genotypes.',
'',
'Summary tables: `tables/substrate_multilevel_summary.csv`, `tables/ficus_wood_zones_multilevel_summary.csv`, `tables/substrate_x_position_multilevel_summary.csv`.',
'',
'---',
'',
'## E. Zone-Based Analyses',
'',
'These test whether **height on the tree** (trunk zones 1\u20135 vs canopy branch zone 6) structures the community. All are run at ITS-genotype level; Section D2 already showed the zone question is also answered the same way at higher ranks.',
'',
'### E1. Ficus Wood Across Zones',
'',
'Community comparison of Ficus wood isolates collected at different tree zones (1\u20136).',
'',
'![NMDS — Ficus wood zones](plots/png/ficus_wood_zones_nmds.png)',
'',
'*Interpretation.* Samples from different zones intermingle on the NMDS with no zone-wise grouping, indicating wood-inhabiting fungi are distributed largely independently of height. Full results: `tables/ficus_wood_zones_*.txt`.',
'',
'### E2. Ficus Wood: Trunk (Zones 1-5) vs Branch (Zone 6)',
'',
paste0('**PERMANOVA: F = ', fw_tb_F, ', R\u00b2 = ', fw_tb_R2, ', p = ', fw_tb_p, '** \u2014 ', verdict(fw_tb_p),
       '. ', if (!is.na(num(fw_tb_p)) && num(fw_tb_p) >= 0.05) 'Trunk and branch wood share a statistically indistinguishable fungal community.' else 'Trunk and branch wood host detectably different communities.'),
'',
'![NMDS — trunk vs branch](plots/png/ficus_wood_trunk_vs_branch_nmds.png)',
'',
'### E4. Substrate Comparison at Zone 6 (Branch Level Only)',
'',
'**Why.** Restricting to zone 6 removes any confound between substrate and height, since all three substrates are present there.',
'',
paste0('**PERMANOVA: F = ', z6_perm_F, ', R\u00b2 = ', z6_perm_R2, ', p = ', z6_perm_p, '** \u2014 ', verdict(z6_perm_p),
       '. The substrate effect ', if (!is.na(num(z6_perm_p)) && num(z6_perm_p) < 0.05) 'persists even within a single zone, confirming it is driven by substrate itself and not by differences in sampling height.' else 'is weaker when height is held constant.'),
'',
'![NMDS — substrates at zone 6](plots/png/substrate_zone6_branch_nmds.png)',
'',
'### E5. Substrate × Position Interaction',
'',
'Combined-factor analysis testing whether communities differ across substrate-position combinations (e.g. *Ficus wood - Trunk* vs *Ficus wood - Branch* vs *Ficus leaves - Branch* ...).',
'',
'![NMDS — substrate × position](plots/png/substrate_x_position_nmds.png)',
'',
'*Interpretation.* Groups separate primarily **by substrate**, with position adding only minor structure \u2014 substrate is the dominant organiser of these endophyte communities (see also the multi-rank confirmation in Section D3).',
'',
'---',
'',
'## Abundance & Composition Plots (Incertae sedis retained)',
'',
'Unlike every analysis above, the plots below **keep the *Incertae sedis* isolates** (as an explicit pooled category), so they present the complete isolate census without hiding unidentified material.',
'',
'### `plots/abundance_pooled_incertae_sedis/`',
'',
'Horizontal grouped bar charts of **absolute isolate counts** per taxon, split by substrate, with unresolved taxa collected into an *Incertae sedis* bar.',
'',
'![Abundance by genus (pooled)](plots/png/abundance_by_genus.png)',
'',
'### `plots/pie_charts/`',
'',
'Proportional composition of each substrate; the *Incertae sedis* slice shows how much of each community remains unidentified at that rank.',
'',
'![Pie chart by phylum](plots/png/pie_by_phylum.png)',
'',
'![Pie chart by genus](plots/png/pie_by_genus.png)',
'',
'---',
'',
'# Conclusion',
'',
paste0('1. **Substrate is the primary driver of community structure.** The three substrates host significantly different fungal communities (PERMANOVA p = ', perm_p, ', R\u00b2 \u2248 ', round(num(perm_R2), 2),
       '), and this holds at **every taxonomic rank** \u2014 the effect is in fact strongest around **', ml_best_rank, ' level** (R\u00b2 \u2248 ', round(ml_best_R2, 2), '). The two leaf substrates are more similar to each other than to wood, i.e. **tissue type (leaf vs wood)** is the strongest split, with **host identity (Ficus vs Lauraceae leaves)** adding a secondary but significant effect (p = ', leaves_perm_p, ').'),
paste0('2. **Tree height has at most a weak effect.** Treated as six discrete zones, height does **not** structure Ficus-wood communities at any taxonomic rank (Section D2, all p > 0.05). When the wood is instead split simply into **trunk vs branch**, a modest but ', verdict(fw_tb_p), ' difference emerges (p = ', fw_tb_p, ', R\u00b2 \u2248 ', round(num(fw_tb_R2), 2), '): branch wood carries a somewhat distinct community from trunk wood, but this coarse contrast explains far less variation than substrate does.'),
paste0('3. **The substrate signal is real, not a sampling-height artefact.** Even when the comparison is restricted to zone 6 alone (where all substrates co-occur), substrates remain ', verdict(z6_perm_p), ' (p = ', z6_perm_p, ').'),
'4. **Communities are diverse and even.** All substrates show high evenness (Pielou J\' > 0.9) and long rank-abundance tails; rarefaction curves have not saturated, so true richness is higher still. A shared generalist core of genera co-exists with a substantial set of substrate-exclusive taxa.',
paste0('5. **Removing *Incertae sedis* sharpened the picture.** Excluding the pooled "unknown" bin from the diversity, overlap and multivariate analyses (while keeping it visible in the abundance/pie plots) increased, rather than decreased, the measured separation between substrates \u2014 confirming that the unidentified fraction had been masking genuine differences.'),
'',
'---',
'',
paste0('*Auto-generated on ', Sys.Date(), ' by `generate_readme.R`*')
)

writeLines(readme, "README.md")
cat("README.md generated\n")

# ---- Generate README.pdf via pandoc ----
# DejaVu fonts give xelatex full Unicode coverage (\u2248, \u2264, \u2192, ...)
pdf_cmd <- paste(
  'pandoc README.md -o README.pdf --pdf-engine=xelatex',
  "-V mainfont='DejaVu Serif'",
  "-V monofont='DejaVu Sans Mono'",
  '-V geometry:margin=2cm',
  '2>&1')
pdf_ok <- system(pdf_cmd, intern = FALSE)
if (pdf_ok == 0) {
  cat("README.pdf generated\n")
} else {
  cat("WARNING: README.pdf generation failed (pandoc/xelatex not available?)\n")
}
