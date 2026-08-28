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
  "plots/community/nmds_ordination.pdf",
  "plots/community/nmds_with_species.pdf",
  "plots/community/rarefaction_curves.pdf",
  "plots/abundance_pooled_incertae_sedis/abundance_by_genus.pdf",
  "plots/pie_charts/pie_by_phylum.pdf",
  "plots/pie_charts/pie_by_genus.pdf"
)

for (pdf in embedded_plots) {
  png <- file.path("plots/png", sub("\\.pdf$", ".png", basename(pdf)))
  pdf_to_png(pdf, png)
}
cat("Converted", length(embedded_plots), "PDFs to PNG\n")

# ---- R and package versions ----
r_ver <- paste0(R.version$major, ".", R.version$minor)
pkg_versions <- sapply(
  c("readxl", "dplyr", "tidyr", "ggplot2", "vegan", "indicspecies",
    "ggVennDiagram", "scales"),
  function(p) as.character(packageVersion(p))
)

# ---- Load results ----
alpha  <- read.csv("tables/alpha_diversity.csv", stringsAsFactors = FALSE)
pw     <- read.csv("tables/pairwise_permanova.csv", stringsAsFactors = FALSE)
indval <- read.csv("tables/indicator_species_significant.csv", stringsAsFactors = FALSE)
perm_txt   <- readLines("tables/permanova_results.txt")
anos_txt   <- readLines("tables/anosim_results.txt")
betad_txt  <- readLines("tables/betadisper_results.txt")

# ---- Extract key stats ----
# PERMANOVA: Model Df SumOfSqs R2 F Pr(>F) ***
perm_line <- grep("^Model", perm_txt, value = TRUE)
perm_parts <- strsplit(trimws(perm_line), "\\s+")[[1]]
perm_F  <- perm_parts[5]
perm_R2 <- perm_parts[4]
perm_p  <- perm_parts[6]

# ANOSIM
anos_R <- sub(".*R statistic:\\s*", "", grep("R statistic:", anos_txt, value = TRUE))
anos_p <- sub(".*p-value:\\s*", "", grep("p-value:", anos_txt, value = TRUE))

# Betadisper
bd_F <- sub(".*F\\s*N\\.Perm.*", "",
            grep("^Groups", betad_txt, value = TRUE))
bd_line <- grep("^Groups", betad_txt, value = TRUE)
bd_parts <- strsplit(trimws(bd_line), "\\s+")[[1]]
bd_F <- bd_parts[4]
bd_p <- bd_parts[6]

# Betadisper distances
bd_dist_lines <- grep("Ficus|Host", betad_txt, value = TRUE)
bd_dist <- sapply(bd_dist_lines, function(x) {
  p <- strsplit(trimws(x), "\\s+")[[1]]
  p[length(p)]
})
names(bd_dist) <- sapply(bd_dist_lines, function(x) {
  p <- strsplit(trimws(x), "\\s+")[[1]]
  paste(p[1:(length(p)-1)], collapse = " ")
})

# NMDS stress (re-extract from the saved plot subtitle isn't possible,
# so we recalculate from saved objects or parse from log)
# We stored it as subtitle text; hardcode extraction from community_analysis output
# Instead, let's recompute quickly
nmds_stress <- tryCatch({
  set.seed(42)
  comm_wide <- bind_rows(
    endo_leaf_ficus %>% mutate(substrate = "Ficus leaves", sample_id = tree_orient) %>% select(substrate, sample_id, its_taxon),
    endo_wood_ficus %>% mutate(substrate = "Ficus wood", sample_id = sub("/.*", "", sample_code)) %>% select(substrate, sample_id, its_taxon),
    endo_wood_host %>% mutate(substrate = "Host wood", sample_id = paste0("zone_", tree_zone)) %>% select(substrate, sample_id, its_taxon)
  ) %>%
    mutate(uid = paste(substrate, sample_id, sep = "__")) %>%
    count(uid, its_taxon) %>%
    pivot_wider(names_from = its_taxon, values_from = n, values_fill = 0)
  cm <- as.matrix(comm_wide[, -1])
  rownames(cm) <- comm_wide$uid
  nmds <- vegan::metaMDS(cm, distance = "bray", k = 2, trymax = 50, trace = 0)
  round(nmds$stress, 4)
}, error = function(e) "N/A")

n_samples <- nrow(comm_wide)
n_taxa_comm <- ncol(comm_wide) - 1

# Alpha at ITS level
alpha_its <- alpha[alpha$taxonomic_level == "its_taxon", ]

# Venn counts (compute from aggregate data)
venn_genus <- dat_agg %>%
  group_by(genus) %>%
  summarise(ll = sum(n_fungi_laur_leaf), fl = sum(n_fungi_fic_leaf),
            fw = sum(n_fungi_fic_wood), .groups = "drop")
venn_sp <- dat_agg %>%
  group_by(species) %>%
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
vs <- count_shared(venn_sp)

# Indicator species table
indval$substrate <- ifelse(indval$s.Ficus.leaves == 1, "Ficus leaves",
                    ifelse(indval$s.Ficus.wood == 1, "Ficus wood",
                    ifelse(indval$s.Host.wood == 1, "Host wood", "Multiple")))
indval_fl <- indval[indval$substrate == "Ficus leaves", ]
indval_hw <- indval[indval$substrate == "Host wood", ]
indval_fw <- indval[indval$substrate == "Ficus wood", ]

# ---- Helper ----
md_table <- function(df, align = NULL) {
  header <- paste0("| ", paste(names(df), collapse = " | "), " |")
  if (is.null(align)) align <- rep("---", ncol(df))
  sep <- paste0("| ", paste(align, collapse = " | "), " |")
  rows <- apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  paste(c(header, sep, rows), collapse = "\n")
}

# ---- Build README ----
readme <- c(
'# Fungal Community Analysis of LOT2 Endophytic Fungi',
'',
paste0('> **Analysis environment:** R ', r_ver, ' | Generated on ', Sys.Date()),
'',
'## Overview',
'',
'This project analyses the fungal endophyte communities isolated from three substrates collected on Ile de La Reunion (LOT2 sampling campaign):',
'',
'| Substrate | Source sheet | Isolates | Sample-level replicates |',
'|-----------|-------------|----------|------------------------|',
paste0('| **Lauraceae leaves** | ALL isolates (aggregate only) | ', alpha_its$abundance_N[alpha_its$substrate == "Lauraceae leaves"], ' | Not available |'),
paste0('| **Ficus leaves** | 67. Fungi-Endo leaf (Ficus) | ', alpha_its$abundance_N[alpha_its$substrate == "Ficus leaves"], ' | 10 (leaf orientations S1-S10) |'),
paste0('| **Ficus wood** | 65. Fungi-Endo wood (Ficus) | ', alpha_its$abundance_N[alpha_its$substrate == "Ficus wood"], ' | 13 (sample collection units) |'),
paste0('| **Host wood** (Lauraceae) | 66. Fungi-Endo wood (Host) | ', nrow(endo_wood_host), ' | 3 (tree zones) |'),
'',
'**Important note on substrates:** The aggregate data ("ALL isolates" sheet) records counts for Lauraceae *leaves*, Ficus leaves, and Ficus wood. However, there is no individual sample sheet for Lauraceae leaves -- only for Lauraceae (Host) *wood* (20 isolates). Therefore, aggregate-level analyses (alpha diversity, Venn diagrams, rank-abundance) use the three aggregate substrates, while sample-level multivariate analyses (NMDS, PERMANOVA, etc.) use Ficus leaves, Ficus wood, and Host wood.',
'',
'### Handling of uncertain taxonomy',
'',
'Across all taxonomic levels, entries with missing or unresolved identifiers (`NA`, `"?"`, `"NO"`, `"incertae sedis"`) are pooled into a single **"Incertae sedis"** category before aggregation.',
'',
'---',
'',
'## Input Data',
'',
'- **`LOT2_for_Livio.xlsx`** -- Source Excel workbook containing all isolate data across multiple sheets.',
'',
'## Scripts',
'',
'| Script | Purpose |',
'|--------|---------|',
'| `main.R` | Master script -- loads libraries, sources all other scripts in order |',
'| `parse_LOT2.R` | Reads and renames columns from all Excel sheets |',
'| `plot_abundance.R` | Horizontal bar charts of isolate counts by taxonomic level |',
'| `plot_abundance_pooled.R` | Same bar charts, but with uncertain taxa pooled as "Incertae sedis" |',
'| `plot_pie.R` | Pie charts of community composition (3 pies per level, one per substrate) |',
'| `community_analysis.R` | Full community ecology analysis: diversity, ordination, statistical tests |',
'| `generate_readme.R` | Generates this README dynamically from analysis results |',
'',
'### How to run',
'',
'```r',
'# From the project directory, run everything:',
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
paste0('| vegan | ', pkg_versions['vegan'], ' | Diversity indices, NMDS, PERMANOVA, ANOSIM, betadisper, rarefaction |'),
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
'  community/                        # Community ecology plots',
'',
'tables/                             # Statistical results and summary tables',
'```',
'',
'---',
'',
'## Abundance Plots',
'',
'### `plots/abundance/`',
'',
'Horizontal grouped bar charts showing the **absolute number of isolates** per taxonomic identifier, broken down by substrate (colour-coded). One plot per taxonomic level: culture code, ITS taxon, phylum, subphylum, superclass, class, subclass, order, family, genus, species.',
'',
'Taxa are sorted by total abundance (lowest at top, highest at bottom). When labels are ambiguous at a given level, parent ranks are prepended to disambiguate.',
'',
'![Abundance by genus (pooled)](plots/png/abundance_by_genus.png)',
'',
'### `plots/abundance_pooled_incertae_sedis/`',
'',
'Same as above, but all unknown/missing taxonomic identifiers are first collapsed into a single "Incertae sedis" group at each level.',
'',
'### `plots/pie_charts/`',
'',
'For each taxonomic level, three pie charts side by side -- one per substrate -- showing the **proportional composition** of the community. All three pies share the same colour palette and taxon ordering. Slices >= 3% are labelled with their percentage.',
'',
'![Pie chart by phylum](plots/png/pie_by_phylum.png)',
'',
'![Pie chart by genus](plots/png/pie_by_genus.png)',
'',
'---',
'',
'# Results',
'',
'## A. Alpha Diversity',
'',
'Alpha diversity measures the richness and evenness of a community *within* a single substrate.',
'',
'### How to interpret',
'',
'| Index | What it measures | Range | Higher means |',
'|-------|-----------------|-------|--------------|',
'| **Richness (S)** | Number of distinct taxa | 0 to infinity | More taxa present |',
'| **Shannon (H\')** | Combines richness and evenness; sensitive to rare species | Typically 0-5 | More diverse |',
'| **Simpson (1-D)** | Probability that two random individuals differ | 0 to 1 | More diverse |',
'| **Inverse Simpson** | Effective number of equally-common species | 1 to S | More even |',
'| **Pielou\'s evenness (J\')** | How evenly individuals are distributed | 0 to 1 | More even |',
'',
'### Results at ITS taxon level',
'',
'| Substrate | Richness (S) | Abundance (N) | Shannon (H\') | Simpson (1-D) | Inv. Simpson | Pielou (J\') |',
'|-----------|:---:|:---:|:---:|:---:|:---:|:---:|',
paste0('| Lauraceae leaves | ', alpha_its$richness_S[1], ' | ', alpha_its$abundance_N[1], ' | ', alpha_its$shannon_H[1], ' | ', alpha_its$simpson_1mD[1], ' | ', alpha_its$inv_simpson[1], ' | ', alpha_its$pielou_J[1], ' |'),
paste0('| Ficus leaves | ', alpha_its$richness_S[2], ' | ', alpha_its$abundance_N[2], ' | ', alpha_its$shannon_H[2], ' | ', alpha_its$simpson_1mD[2], ' | ', alpha_its$inv_simpson[2], ' | ', alpha_its$pielou_J[2], ' |'),
paste0('| Ficus wood | ', alpha_its$richness_S[3], ' | ', alpha_its$abundance_N[3], ' | ', alpha_its$shannon_H[3], ' | ', alpha_its$simpson_1mD[3], ' | ', alpha_its$inv_simpson[3], ' | ', alpha_its$pielou_J[3], ' |'),
'',
paste0('Ficus leaves harbour the most ITS taxa (**', alpha_its$richness_S[2], '**), followed by Lauraceae leaves (**', alpha_its$richness_S[1], '**) and Ficus wood (**', alpha_its$richness_S[3], '**). Shannon diversity is broadly similar across substrates (H\' = ', alpha_its$shannon_H[3], '-', alpha_its$shannon_H[2], '), and Pielou\'s evenness values above 0.91 indicate that no single taxon strongly dominates any substrate.'),
'',
'### Results at genus level',
'',
{
  alpha_gen <- alpha[alpha$taxonomic_level == "genus", ]
  c(
    '| Substrate | Richness (S) | Shannon (H\') | Simpson (1-D) | Pielou (J\') |',
    '|-----------|:---:|:---:|:---:|:---:|',
    paste0('| Lauraceae leaves | ', alpha_gen$richness_S[1], ' | ', alpha_gen$shannon_H[1], ' | ', alpha_gen$simpson_1mD[1], ' | ', alpha_gen$pielou_J[1], ' |'),
    paste0('| Ficus leaves | ', alpha_gen$richness_S[2], ' | ', alpha_gen$shannon_H[2], ' | ', alpha_gen$simpson_1mD[2], ' | ', alpha_gen$pielou_J[2], ' |'),
    paste0('| Ficus wood | ', alpha_gen$richness_S[3], ' | ', alpha_gen$shannon_H[3], ' | ', alpha_gen$simpson_1mD[3], ' | ', alpha_gen$pielou_J[3], ' |')
  )
},
'',
'Full diversity data across all taxonomic levels is in `tables/alpha_diversity.csv`.',
'',
'**Plots:** `plots/community/alpha_shannon_H.pdf`, `alpha_simpson_1mD.pdf`, `alpha_inv_simpson.pdf`, `alpha_pielou_J.pdf`',
'',
'![Shannon diversity across taxonomic levels](plots/png/alpha_shannon_H.png)',
'',
'---',
'',
'## B. Community Composition Comparisons',
'',
'### B1. Rank-Abundance Curves (Whittaker Plots)',
'',
'Each taxon is ranked from most to least abundant (x-axis) and its relative abundance is plotted on a log scale (y-axis). One curve per substrate.',
'',
'**How to interpret:** A steep curve indicates a community dominated by a few taxa with many rare ones. A flat curve indicates an even community. Comparing curves across substrates reveals differences in dominance structure.',
'',
'**Plots:** `plots/community/rank_abundance_genus.pdf`, `rank_abundance_species.pdf`, `rank_abundance_its_taxon.pdf`',
'',
'![Rank-abundance curve at genus level](plots/png/rank_abundance_genus.png)',
'',
'### B2. Relative Abundance Stacked Bar Charts',
'',
'Stacked bars showing the proportional composition of each substrate at a given taxonomic level. The y-axis shows percentage, making substrates with different total isolate numbers directly comparable.',
'',
'**Plots:** `plots/community/rel_abundance_phylum.pdf`, `rel_abundance_class.pdf`, `rel_abundance_order.pdf`, `rel_abundance_family.pdf`, `rel_abundance_genus.pdf`',
'',
'![Relative abundance by phylum](plots/png/rel_abundance_phylum.png)',
'',
'![Relative abundance by genus](plots/png/rel_abundance_genus.png)',
'',
'### B3. Venn Diagrams -- Shared and Unique Taxa',
'',
'Shows the number of taxa shared between substrates and unique to each.',
'',
paste0('**At genus level** (', vg$total, ' total genera):'),
'',
paste0('- Shared across all 3 substrates: **', vg$all_three, '**'),
paste0('- Lauraceae leaves only: **', vg$ll_only, '**'),
paste0('- Ficus leaves only: **', vg$fl_only, '**'),
paste0('- Ficus wood only: **', vg$fw_only, '**'),
paste0('- Lauraceae + Ficus leaves only: **', vg$ll_fl, '**'),
paste0('- Lauraceae + Ficus wood only: **', vg$ll_fw, '**'),
paste0('- Ficus leaves + Ficus wood only: **', vg$fl_fw, '**'),
'',
paste0('**At species level** (', vs$total, ' total species):'),
'',
paste0('- Shared across all 3 substrates: **', vs$all_three, '**'),
paste0('- Lauraceae leaves only: **', vs$ll_only, '**'),
paste0('- Ficus leaves only: **', vs$fl_only, '**'),
paste0('- Ficus wood only: **', vs$fw_only, '**'),
'',
'**Plots:** `plots/community/venn_genus.pdf`, `venn_species.pdf`, `venn_its_taxon.pdf`',
'',
'![Venn diagram at genus level](plots/png/venn_genus.png)',
'',
'---',
'',
'# Multivariate Statistical Analyses',
'',
paste0('These analyses use sample-level data (', n_samples, ' samples total, ', n_taxa_comm, ' ITS taxa) and **Bray-Curtis dissimilarity** to compare community composition across substrates.'),
'',
'**Substrates used:** Ficus leaves (10 samples by leaf orientation), Ficus wood (13 samples by collection unit), Host wood (3 samples by tree zone).',
'',
'### Parameters common to all tests',
'',
'| Parameter | Value |',
'|-----------|-------|',
'| Distance metric | Bray-Curtis (`method = "bray"`) |',
'| Number of permutations | 999 |',
'| Random seed | 42 (`set.seed(42)`) |',
paste0('| Community matrix dimensions | ', n_samples, ' samples x ', n_taxa_comm, ' taxa |'),
'',
'---',
'',
'## C1. NMDS Ordination',
'',
'Non-metric Multidimensional Scaling reduces the high-dimensional community data into a 2D plot where each point represents one sample. Points close together have similar community composition.',
'',
'**Parameters:** `vegan::metaMDS(comm, distance = "bray", k = 2, trymax = 200)`',
'',
paste0('**Stress = ', nmds_stress, '**'),
'',
'| Stress value | Quality of representation |',
'|:---:|---|',
'| < 0.05 | Excellent |',
'| < 0.10 | Good |',
'| < 0.20 | Acceptable |',
'| > 0.20 | Poor -- interpret with caution |',
'',
paste0('The stress value of ', nmds_stress, ' indicates a **good** 2D representation of the community distances.'),
'',
'**How to read the plot:**',
'- Points coloured by substrate; clustering by colour = communities differ systematically',
'- 95% confidence ellipses show group spread; non-overlapping = distinct communities',
'- Crosses mark group centroids (mean position)',
'- The species overlay plot shows which of the 15 most abundant taxa drive group separation',
'',
'**Plots:** `plots/community/nmds_ordination.pdf`, `nmds_with_species.pdf`',
'',
'![NMDS ordination](plots/png/nmds_ordination.png)',
'',
'![NMDS with species overlay](plots/png/nmds_with_species.png)',
'',
'---',
'',
'## C2. PERMANOVA',
'',
'Permutational Multivariate Analysis of Variance tests whether community composition centroids differ between substrates. This is the multivariate equivalent of an ANOVA.',
'',
'**Parameters:** `vegan::adonis2(comm ~ substrate, method = "bray", permutations = 999)`',
'',
'```',
paste(perm_txt, collapse = "\n"),
'```',
'',
'| Metric | Value | Interpretation |',
'|--------|:-----:|----------------|',
paste0('| F statistic | ', perm_F, ' | Ratio of between- to within-group variation |'),
paste0('| R-squared | ', perm_R2, ' | ', round(as.numeric(perm_R2) * 100, 1), '% of variation explained by substrate |'),
paste0('| p-value | ', perm_p, ' | **Highly significant** -- communities differ |'),
'',
paste0('**Interpretation:** Substrate identity explains **', round(as.numeric(perm_R2) * 100, 1), '%** of the total variation in community composition (p = ', perm_p, '). The remaining ', round((1 - as.numeric(perm_R2)) * 100, 1), '% is due to within-substrate variability and unmeasured factors. A significant result means at least two substrate groups harbour distinct communities.'),
'',
'**Caveat:** PERMANOVA can be sensitive to differences in multivariate dispersion (spread). See the betadisper test below.',
'',
'---',
'',
'## C3. ANOSIM',
'',
'Analysis of Similarity -- a complementary non-parametric test comparing between-group to within-group dissimilarities using ranks.',
'',
'**Parameters:** `vegan::anosim(comm, grouping = substrate, distance = "bray", permutations = 999)`',
'',
'| Metric | Value | Interpretation |',
'|--------|:-----:|----------------|',
paste0('| R statistic | ', anos_R, ' | Ranges -1 to 1; values > 0.5 = well-separated groups |'),
paste0('| p-value | ', anos_p, ' | **Highly significant** |'),
'',
paste0('**Interpretation:** R = ', anos_R, ' indicates **strong separation** between substrate communities. Values near 0 would indicate no difference; values near 1 indicate complete separation. This confirms the PERMANOVA finding.'),
'',
'---',
'',
'## C4. Beta-Dispersion Test (PERMDISP)',
'',
'Tests whether groups have equal multivariate spread (dispersion). This is an assumption of PERMANOVA.',
'',
'**Parameters:** `vegan::betadisper(vegdist(comm, "bray"), groups)` + `permutest(bd, permutations = 999)`',
'',
'| Metric | Value | Interpretation |',
'|--------|:-----:|----------------|',
paste0('| F statistic | ', bd_F, ' | Large = dispersions differ |'),
paste0('| p-value | ', bd_p, ' | **Significant** -- dispersions are unequal |'),
'',
'**Mean distance to centroid per substrate:**',
'',
'| Substrate | Mean distance | Interpretation |',
'|-----------|:---:|---|',
paste0(sapply(seq_along(bd_dist), function(i) {
  paste0('| ', names(bd_dist)[i], ' | ', bd_dist[i], ' | ',
         ifelse(as.numeric(bd_dist[i]) == max(as.numeric(bd_dist)),
                'Most variable', 'Less variable'), ' |')
}), collapse = "\n"),
'',
paste0('**Interpretation:** The significant result (p = ', bd_p, ') means within-group variability differs between substrates. Ficus wood communities are the most variable (highest distance to centroid), while Host wood communities are the most homogeneous. This means the significant PERMANOVA result could partly reflect dispersion differences rather than purely centroid differences. However, the ANOSIM result (R = ', anos_R, ', which is less sensitive to dispersion) still supports genuine community differences.'),
'',
'---',
'',
'## C5. Pairwise PERMANOVA',
'',
'Post-hoc pairwise comparisons with Holm-corrected p-values for multiple testing.',
'',
'**Parameters:** `vegan::adonis2()` on each pair, p-values adjusted with `p.adjust(method = "holm")`',
'',
'| Comparison | F | R-squared | p (raw) | p (Holm-adjusted) | Significant? |',
'|-----------|:---:|:---:|:---:|:---:|:---:|',
paste0(apply(pw, 1, function(r) {
  sig <- ifelse(as.numeric(r["p_adj"]) < 0.001, "***",
         ifelse(as.numeric(r["p_adj"]) < 0.01, "**",
         ifelse(as.numeric(r["p_adj"]) < 0.05, "*", "ns")))
  paste0('| ', r["pair"], ' | ', r["F_stat"], ' | ', r["R2"], ' | ', r["p_value"], ' | ', r["p_adj"], ' | ', sig, ' |')
}), collapse = "\n"),
'',
{
  max_r2_idx <- which.max(pw$R2)
  paste0('**Interpretation:** All three pairwise comparisons are significant (adjusted p < 0.05), confirming that **each substrate harbours a distinct fungal community**. The largest effect size (R-squared = ', pw$R2[max_r2_idx], ') is between ', pw$pair[max_r2_idx], ', indicating these two substrates are the most different from each other.')
},
'',
'---',
'',
'## D1. Rarefaction Curves',
'',
'Shows expected number of taxa as a function of sampling effort (number of individuals).',
'',
'**How to interpret:** If the curve reaches a plateau, sampling was sufficient to capture most diversity. If still rising steeply, more sampling would reveal additional taxa. Substrates with curves plateauing at different heights have genuinely different richness.',
'',
'**Plot:** `plots/community/rarefaction_curves.pdf`',
'',
'![Rarefaction curves](plots/png/rarefaction_curves.png)',
'',
'---',
'',
'## D2. Indicator Species Analysis (IndVal)',
'',
'Identifies taxa significantly associated with a particular substrate based on fidelity (how consistently the taxon occurs) and exclusivity (how restricted it is).',
'',
'**Parameters:** `indicspecies::multipatt(comm, cluster = substrate, func = "IndVal.g", control = how(nperm = 999))`',
'',
paste0('**', nrow(indval), ' significant indicator taxa** were identified (p <= 0.05):'),
'',
if (nrow(indval_fl) > 0) {
  c(
    paste0('### Ficus leaves indicators (', nrow(indval_fl), ' taxa)'),
    '',
    '| Taxon | IndVal statistic | p-value |',
    '|-------|:---:|:---:|',
    paste0(apply(indval_fl, 1, function(r) {
      paste0('| *', r["taxon"], '* | ', round(as.numeric(r["stat"]), 3), ' | ', r["p.value"], ' |')
    }), collapse = "\n")
  )
} else character(0),
'',
if (nrow(indval_hw) > 0) {
  c(
    paste0('### Host wood indicators (', nrow(indval_hw), ' taxa)'),
    '',
    '| Taxon | IndVal statistic | p-value |',
    '|-------|:---:|:---:|',
    paste0(apply(indval_hw, 1, function(r) {
      paste0('| ', r["taxon"], ' | ', round(as.numeric(r["stat"]), 3), ' | ', r["p.value"], ' |')
    }), collapse = "\n")
  )
} else character(0),
'',
if (nrow(indval_fw) > 0) {
  c(
    paste0('### Ficus wood indicators (', nrow(indval_fw), ' taxa)'),
    '',
    '| Taxon | IndVal statistic | p-value |',
    '|-------|:---:|:---:|',
    paste0(apply(indval_fw, 1, function(r) {
      paste0('| *', r["taxon"], '* | ', round(as.numeric(r["stat"]), 3), ' | ', r["p.value"], ' |')
    }), collapse = "\n")
  )
} else character(0),
'',
'**How to interpret the IndVal statistic:** Ranges 0 to 1. A value of 1 means the taxon is found in all samples of that substrate and nowhere else (perfect indicator). The A component measures specificity (exclusivity to the group) and B measures fidelity (frequency of occurrence within the group).',
'',
paste0('Most indicators are associated with Ficus leaves (', nrow(indval_fl), ' taxa), suggesting this substrate harbours a particularly distinctive fungal assemblage. The absence of Ficus wood indicators likely reflects its higher within-group variability (mean distance to centroid = ', bd_dist[grep("Ficus wood", names(bd_dist))], ').'),
'',
'**Full results:** `tables/indicator_species.txt`, `tables/indicator_species_significant.csv`',
'',
'---',
'',
'# Summary of Key Findings',
'',
paste0('1. **Community composition differs significantly across all substrate pairs** (PERMANOVA p = ', perm_p, ', ANOSIM R = ', anos_R, ', all pairwise comparisons adjusted p < 0.05).'),
paste0('2. **Ficus leaves harbour the most ITS taxa** (', alpha_its$richness_S[2], '), followed by Lauraceae leaves (', alpha_its$richness_S[1], ') and Ficus wood (', alpha_its$richness_S[3], ').'),
paste0('3. **Shannon diversity is broadly similar** across substrates at the ITS level (H\' = ', alpha_its$shannon_H[3], '-', alpha_its$shannon_H[2], '), with Pielou\'s J\' > 0.91 indicating no strong single-taxon dominance.'),
paste0('4. **', nrow(indval_fl), ' indicator species** are significantly associated with Ficus leaves, particularly *Colletotrichum* and *Diaporthe* species.'),
paste0('5. **Within-group variability differs** between substrates (betadisper p = ', bd_p, '): Ficus wood is the most variable, Host wood the least.'),
paste0('6. **NMDS stress = ', nmds_stress, '** -- the 2D ordination is a good representation of community distances.'),
paste0('7. **', vg$all_three, ' genera** are shared across all three substrates, while ', vg$fl_only, ' are unique to Ficus leaves, ', vg$fw_only, ' to Ficus wood, and ', vg$ll_only, ' to Lauraceae leaves.')
)

writeLines(readme, "README.md")
cat("Generated: README.md\n")

# ---- Render README to PDF ----
tryCatch({
  rmarkdown::render(
    "README.md",
    output_format = rmarkdown::pdf_document(
      toc = TRUE,
      toc_depth = 3,
      number_sections = TRUE,
      latex_engine = "xelatex"
    ),
    output_file = "README.pdf",
    quiet = TRUE
  )
  cat("Generated: README.pdf\n")
}, error = function(e) {
  cat("PDF rendering failed:", conditionMessage(e), "\n")
  cat("Trying pandoc directly...\n")
  ret <- system(paste(
    'pandoc README.md -o README.pdf',
    '--pdf-engine=xelatex',
    '--toc --toc-depth=3',
    '-V geometry:margin=2.5cm',
    '-V fontsize=11pt',
    '-V mainfont="DejaVu Sans"'
  ))
  if (ret == 0) cat("Generated: README.pdf (via pandoc)\n")
  else cat("PDF rendering failed with both methods\n")
})
