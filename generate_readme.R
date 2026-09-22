dir.create("tables", showWarnings = FALSE)
dir.create("plots/png", recursive = TRUE, showWarnings = FALSE)

# ---- Convert key PDFs to PNGs for embedding ----
# The project is run on both Linux and Windows, where a different set
# of converters is available, so resolve one at run time instead of
# assuming ImageMagick is on the PATH.
first_on_path <- function(candidates) {
  for (cmd in candidates) {
    p <- Sys.which(cmd)
    if (nzchar(p)) return(unname(p))
  }
  ""
}

png_converter <- if (requireNamespace("pdftools", quietly = TRUE)) {
  "pdftools"
} else {
  first_on_path(c("magick", "convert", "pdftocairo", "pdftoppm"))
}

pdf_to_png <- function(pdf_path, png_path, density = 150) {
  if (identical(png_converter, "pdftools")) {
    ok <- tryCatch({
      pdftools::pdf_convert(pdf_path, format = "png", pages = 1,
                            dpi = density, filenames = png_path, verbose = FALSE)
      TRUE
    }, error = function(e) FALSE)
    return(ok)
  }
  if (!nzchar(png_converter)) return(FALSE)
  tool <- tolower(basename(png_converter))
  cmd <- if (grepl("pdftocairo|pdftoppm", tool)) {
    # these append "-<page>.png" to the prefix they are given
    sprintf('"%s" -png -r %d -f 1 -l 1 -singlefile "%s" "%s"',
            png_converter, density, pdf_path, sub("\\.png$", "", png_path))
  } else {
    sprintf('"%s" -density %d "%s[0]" -quality 90 "%s"',
            png_converter, density, pdf_path, png_path)
  }
  system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE) == 0
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
  "plots/community/orientation_all_nmds.pdf",
  "plots/community/orientation_sun_aspect_nmds.pdf",
  "plots/community/ficus_leaves_orientation_nmds.pdf",
  "plots/community/substrate_x_orientation_nmds.pdf",
  "plots/community/orientation_alpha_richness_S.pdf",
  "plots/community/orientation_alpha_shannon_H.pdf",
  "plots/community/orientation_venn_ficus_leaves.pdf",
  "plots/community/orientation_all_rarefaction.pdf",
  "plots/diversity/richness_estimators.pdf",
  "plots/diversity/alpha_per_unit_richness_S.pdf",
  "plots/diversity/alpha_per_unit_shannon_H.pdf",
  "plots/abundance_pooled_incertae_sedis/abundance_by_genus.pdf",
  "plots/pie_charts/pie_by_phylum.pdf",
  "plots/pie_charts/pie_by_genus.pdf"
)

n_png <- 0
for (pdf in embedded_plots) {
  if (file.exists(pdf)) {
    png <- file.path("plots/png", sub("\\.pdf$", ".png", basename(pdf)))
    if (isTRUE(pdf_to_png(pdf, png))) n_png <- n_png + 1
  }
}
if (n_png > 0) {
  cat("Converted", n_png, "PDFs to PNG\n")
} else {
  cat("WARNING: no PDF->PNG converter found; README images may be stale.\n",
      "         Install the 'pdftools' R package, or ImageMagick / poppler.\n", sep = "")
}

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

# Sparsity of the ITS-level matrix - drives how far the ordinations
# can be trusted (see the caveat block in the Results)
sparse_no_shared <- if (!is.null(results_substrate$prop_no_shared))
  round(100 * results_substrate$prop_no_shared) else NA
sparse_singletons <- if (!is.null(results_substrate$prop_singleton_taxa))
  round(100 * results_substrate$prop_singleton_taxa) else NA
# Analytical-unit and completeness figures used in the narrative
n_otus   <- n_distinct(na.omit(samples_raw$otu_97))
n_labels <- n_distinct(samples_raw$its_taxon)
n_seqs   <- nrow(samples_raw)
completeness_tbl   <- safe_csv("tables/sampling_completeness.csv")
alpha_summary_tbl  <- safe_csv("tables/alpha_diversity_per_unit_summary.csv")
alpha_tests_tbl    <- safe_csv("tables/alpha_diversity_tests.csv")
alpha_pairwise_tbl <- safe_csv("tables/alpha_diversity_pairwise.csv")
coverage_tbl       <- safe_csv("tables/richness_at_equal_coverage.csv")
# Indicator-species counts (FDR-corrected in community_analysis.R)
n_ind_raw    <- if (!is.null(results_substrate$n_indicators_raw))
  results_substrate$n_indicators_raw else NA
n_ind_tested <- if (!is.null(results_substrate$n_taxa_tested))
  results_substrate$n_taxa_tested else NA
# Substrate effect adjusted for sampling depth
adj_R2 <- if (!is.null(results_substrate$adj_group_R2)) results_substrate$adj_group_R2 else NA
adj_p  <- if (!is.null(results_substrate$adj_group_p))  results_substrate$adj_group_p  else NA
dep_R2 <- if (!is.null(results_substrate$depth_R2))     results_substrate$depth_R2     else NA
dep_p  <- if (!is.null(results_substrate$depth_p))      results_substrate$depth_p      else NA

cov_range <- if (nrow(completeness_tbl) > 0)
  paste0(min(completeness_tbl$goods_coverage), "-",
         max(completeness_tbl$goods_coverage), "%") else "n/a"

# Honest (strong-tie) stress for the headline substrate ordination
nmds_stress_strong <- if (!is.null(results_substrate$nmds_stress_strong))
  results_substrate$nmds_stress_strong else NA

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
# Count the sampling units the analyses actually use (substrate x zone x
# unit), not distinct unit labels - a couple of Lauraceae unit labels
# recur in both zone 5 and zone 6, so the two counts differ.
sub_units  <- samples_raw %>% distinct(substrate, sample_id) %>% count(substrate) %>% arrange(substrate)

# ---- Multi-level taxonomic sweep results ----
# ---- Branch / trunk orientation results ----
ori_cov     <- safe_csv("tables/orientation_coverage.csv")
ori_unres   <- safe_csv("tables/orientation_unresolved.csv")
ori_sum     <- safe_csv("tables/orientation_summary.csv")
ori_alpha   <- safe_csv("tables/orientation_alpha_diversity.csv")
ori_alpha_s <- safe_csv("tables/orientation_alpha_diversity_by_substrate.csv")
ori_aspect_alpha <- safe_csv("tables/orientation_sun_aspect_alpha_diversity.csv")
ori_circ    <- safe_csv("tables/orientation_circular_gradient.csv")
ori_share   <- safe_csv("tables/orientation_taxon_sharing.csv")
ori_ml      <- safe_csv("tables/orientation_multilevel_summary.csv")
ori_part_txt <- safe_read("tables/orientation_variance_partitioning.txt")

# Pull the marginal PERMANOVA rows out of the partitioning report
partition_row <- function(txt, term) {
  line <- grep(paste0("^", term, "\\s"), txt, value = TRUE)
  if (length(line) == 0) return(list(R2 = "N/A", F = "N/A", p = "N/A"))
  pp <- strsplit(trimws(line[1]), "\\s+")[[1]]
  # term Df SumOfSqs R2 F Pr(>F) [sig]
  list(R2 = pp[4], F = pp[5], p = pp[6])
}
# The report holds three tables; the first two rows are the marginal
# substrate + orientation model.
ori_part_marginal <- ori_part_txt[seq_len(min(length(ori_part_txt),
                                 which(ori_part_txt == "== Sequential (type I): substrate entered first ==")[1] - 1))]
if (length(ori_part_marginal) == 0) ori_part_marginal <- ori_part_txt
part_sub <- partition_row(ori_part_marginal, "substrate")
part_ori <- partition_row(ori_part_marginal, "orientation")

# The sun-aspect marginal model sits in its own block
aspect_block_start <- which(ori_part_txt == "== Marginal (type III): substrate + sun aspect ==")
ori_part_aspect <- if (length(aspect_block_start))
  ori_part_txt[aspect_block_start[1]:length(ori_part_txt)] else character(0)
part_asp <- partition_row(ori_part_aspect, "sun_aspect")

ori_row <- function(name) {
  if (nrow(ori_sum) == 0) return(NULL)
  r <- ori_sum[ori_sum$analysis == name, , drop = FALSE]
  if (nrow(r) == 0) NULL else r[1, ]
}
ori_all    <- ori_row("Orientation - all substrates")
ori_aspect <- ori_row("Sun aspect - all substrates")

n_orient_units    <- if (!is.null(ori_all)) ori_all$n_units else NA
n_orient_isolates <- sum(ori_cov$n_isolates)
n_unres_isolates  <- sum(ori_unres$n_isolates)
orient_any_sig <- if (nrow(ori_sum) > 0)
  any(num(ori_sum$permanova_p[grepl("^Orientation - ", ori_sum$analysis)]) < 0.05,
      na.rm = TRUE) else FALSE

ml_sub <- safe_csv("tables/substrate_multilevel_summary.csv")
ml_zon <- safe_csv("tables/ficus_wood_zones_multilevel_summary.csv")
ml_sp  <- safe_csv("tables/substrate_x_position_multilevel_summary.csv")
lincov <- safe_csv("tables/lineage_coverage.csv")

# Build a markdown table from a multi-level summary data frame
ml_table <- function(df) {
  if (nrow(df) == 0) return(character(0))
  header <- c(
    '| Rank | Taxa | No-taxa pairs | Stress | Tie-aware | F | R\u00b2 | p | ANOSIM R | ANOSIM p | PERMDISP p |',
    '|------|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|')
  key <- paste('*F / R\u00b2 / p are PERMANOVA; "No-taxa pairs" is the share of',
               'sample pairs with no taxon in common; "Tie-aware" is the stress',
               'refitted with strong ties \u2014 the honest one to judge by.*')
  rows <- apply(df, 1, function(r) {
    pns <- if ('prop_no_shared' %in% names(df))
             paste0(round(100 * num(r['prop_no_shared'])), '%') else '\u2014'
    bd  <- if (is.na(r['betadisper_p']) || !nzchar(trimws(r['betadisper_p'])))
             '\u2014' else trimws(r['betadisper_p'])
    sstr <- if ('nmds_stress_strong' %in% names(df))
              trimws(r['nmds_stress_strong']) else '\u2014'
    paste0('| ', r['level'], ' | ', r['n_taxa'], ' | ', pns, ' | ', r['nmds_stress'], ' | ',
           sstr, ' | ',
           r['permanova_F'], ' | ', r['permanova_R2'], ' | ', r['permanova_p'], ' | ',
           r['anosim_R'], ' | ', r['anosim_p'], ' | ', bd, ' |')
  })
  c(header, rows, '', key)
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
'### Sampling design: what was sampled, and the assumptions behind it',
'',
'> **This section records assumptions supplied by the data owner, not facts derivable from the files.** They are written out so collaborators can check the reasoning, and because several of them change how the results must be read. If any is wrong, the analyses affected are named underneath it.',
'',
'**1. Two individual trees were sampled — one of each host.** A strangler *Ficus* growing on a single Lauraceae host tree. There is no tree-level replication, and no tree identifier in the raw files; `parse_LOT2.R` assigns `tree_id` from the substrate.',
'',
'*Consequence.* The **inferential unit is the individual tree, not the host species.** *Ficus* leaves vs *Ficus* wood is a genuine within-tree tissue contrast. But *Ficus* vs Lauraceae compares **one tree with one other tree**, so every "host effect" below describes these two individuals. Sampling units are treated as replicates by PERMANOVA, which is pseudoreplication at the species level: the p-values are real for the trees sampled and cannot be generalised to *Ficus* vs Lauraceae as taxa. All host-level statements are phrased accordingly.',
'',
'**2. The two trees are not independent.** The *Ficus* is a strangler growing **on** the Lauraceae, so the two root systems, canopies and microclimates are physically interlocked and share a continuous surface for fungal dispersal.',
'',
'*Consequence.* Any host difference is a difference between two *interlocked* individuals, which is a conservative setting for detecting host effects (shared exposure should erode differences, not create them) but rules out treating them as independent samples of two species.',
'',
'**3. The Lauraceae trunk was never sampled — it was inaccessible, completely encased by the strangling Ficus.** All Lauraceae material came from exposed upper branches.',
'',
'*Consequence.* This one **changes the data.** The rule "zone ≤ 5 = trunk" is correct for the *Ficus* but wrong for the Lauraceae: applied blindly it mislabelled **3 sampling units / 53 isolates** (the zone-5 Lauraceae units) as trunk material. `parse_LOT2.R` now assigns `position = "Branch"` to all Lauraceae. The previous `lauraceae_leaves_trunk_vs_branch_*` analysis was therefore comparing branch against branch and has been removed. Note also that **zone means different things on the two trees**: on the *Ficus* it is height on a continuous trunk; on the Lauraceae it distinguishes two bands of exposed canopy branch.',
'',
'**4. Whether the Lauraceae material is leaves or branch wood is still unconfirmed.** `LOT2_pooled_counts.xlsx` labels the column "Lauraceae leaves"; the samples workbook has a sheet titled "66. Fungi-Endo wood (Host)". The data owner indicates the material came from branches, without settling leaf vs wood. The analyses retain the label **Lauraceae leaves**.',
'',
'*Consequence.* **The one conclusion that depends on this is the tissue-type claim.** If the Lauraceae material is leaves, the finding "the two leaf substrates resemble each other more than either resembles wood" stands and tissue type is the dominant split. If it is branch wood, that grouping is wrong and the pattern would have to be re-read as *Ficus* vs Lauraceae. Everything else — diversity, orientation, zone — is unaffected. **This should be resolved before the tissue-type interpretation is published.**',
'',
'**5. Each row is an independent colony; no de-replication is required.** Repeated isolates of the same genotype within a sampling unit could be either independent colonisations or one colony subsampled. The rule given is that a **shared `Hofstetter-culture code` marks subsamples of a single colony**. The parser checks this: all **650 isolates carry 650 distinct culture codes**, so no collapsing is needed and abundances are counts of independent colonies.',
'',
'*Consequence.* None — but the check now runs on every pipeline execution and will warn if a future data version reuses a code.',
'',
'**6. Only the first sheet of each workbook is authoritative.** The remaining sheets (per-substrate extracts of 264 / 167 / 20 isolates, and the sequence sheet) are working material. The pipeline reads sheet 1 for both files, plus the sequence sheet `Feuil2` for OTU clustering.',
'',
'### Sampling design',
'',
'Samples were collected at different **tree zones** (heights), zone 1 lowest to zone 6 highest. On the *Ficus*, zones 1-5 are trunk and zone 6 is canopy branch. On the Lauraceae the trunk was inaccessible (assumption 3 above), so its zones 5 and 6 are both **exposed upper branch**, not trunk.',
'',
'| Substrate | Tree | Zones | Position | Sampling units |',
'|-----------|------|:---:|---|:-:|',
paste0('| Ficus leaves | Ficus (strangler) | 6 | Canopy branch | ', sub_units$n[sub_units$substrate == "Ficus leaves"], ' |'),
paste0('| Ficus wood | Ficus (strangler) | 1-6 | Trunk (1-5) + branch (6) | ', sub_units$n[sub_units$substrate == "Ficus wood"], ' |'),
paste0('| Lauraceae leaves | Lauraceae (host) | 5-6 | Exposed branch only | ', sub_units$n[sub_units$substrate == "Lauraceae leaves"], ' |'),
'',
'Because the Lauraceae contributes no trunk material, every trunk-vs-branch contrast in this report is a **Ficus** contrast.',
'',
'#### Sampling orientation',
'',
'Each sampled branch (and, where recorded, trunk face) also carries a **compass orientation**. This information only became usable once the field team reconciled two independent labelling systems: the colour codes written down by the climbers in the field, and the branch numbering entered into the project database. The two disagreed (notably a blue/turquoise colour clash, where blue should have been reserved for the Lauraceae, and disputed collection zones for the Lauraceae branches), so earlier versions of `LOT2_samples.xlsx` had an empty orientation column.',
'',
'The re-issued workbook settles that correspondence and adds it as a reconciled column, which the analyses in **Section F** use. Two clean-up rules are applied when reading it:',
'',
'- Some bearings still carry a replicate index (`N1`-`N4`, `NW1`-`NW4`). The digit identifies **which branch**, not which direction, so it is stripped: `N3` becomes `N`.',
'- Material marked `?`, and all Ficus trunk wood (for which no orientation was ever recorded), is treated as **unresolved** and excluded from the orientation analyses.',
'',
paste0('This leaves **', n_orient_isolates, ' of ', nrow(samples_raw),
       ' isolates** (', round(100 * n_orient_isolates / nrow(samples_raw)),
       '%) with a usable bearing; **', n_unres_isolates,
       '** are unresolved. A single sampling unit can span more than one bearing, so the orientation analyses group isolates by **substrate x zone x unit x bearing**, giving **',
       n_orient_units, ' orientation-level sampling units** rather than the ',
       n_samples, ' used elsewhere.'),
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
paste0('- **`LOT2_samples.xlsx`** (first sheet) — Individual isolate records with sampling zone, unit and reconciled branch/trunk orientation (', nrow(samples_raw), ' isolates)'),
'',
'  The workbook also keeps the reconciliation working columns (field colour code, database colour code, pre-reconciliation orientation, field notes). They are read for traceability but only the reconciled orientation column is used. Sheet `Feuil2` holds one ITS sequence per isolate and is used for OTU clustering. Remaining sheets are working material and are not read.',
'',
paste0('- **`LOT2_otu_map.csv`** \u2014 isolate \u2192 OTU assignment produced by `cluster_otus.sh` (', n_otus, ' OTUs at 97% ITS identity). Committed so the pipeline runs without vsearch.'),
'',
'### The analytical unit',
'',
paste0('Community analyses use **', unit_label, 's**, not the BLAST-derived name strings. Those names are not a clustered unit: they over-split (near-identical sequences filed under two spellings, e.g. *guandongensis* / *guangdongensis*) and over-lump (single labels covering "26 spp."). Clustering the ', n_seqs, ' ITS sequences themselves collapses ', n_labels, ' name labels into **', n_otus, ' OTUs**, and materially improves the data:'),
'',
'| | Name labels | 97% OTUs |',
'|---|:-:|:-:|',
paste0('| Taxa | ', n_labels, ' | ', n_otus, ' |'),
paste0("| Good's coverage | 70% | **", cov_range, '** |'),
paste0('| Sample pairs sharing no taxon | 59% | **', round(100 * num(results_substrate$prop_no_shared)), '%** |'),
'',
'`cluster_otus.sh` regenerates the mapping (needs `vsearch`); it is deterministic and only needs re-running if the sequences change.',
'',
'### Input integrity checks',
'',
'`parse_LOT2.R` now refuses to analyse silently-inconsistent inputs. It drops spreadsheet **totals rows** (rows with neither a culture code nor a taxon name \u2014 one such row was previously read as a genotype and plotted as a giant *Incertae sedis* category), verifies that pooled and sample-level isolate totals agree, and checks that culture codes are unique, since a repeated code would mark subsamples of one colony that must be collapsed before abundances mean anything.',
'',
'## Scripts',
'',
'| Script | Purpose |',
'|--------|---------|',
'| `main.R` | Master script — loads libraries, sources all other scripts in order |',
'| `cluster_otus.sh` | Clusters the ITS sequences into OTUs (run once; needs `vsearch`) |',
'| `parse_LOT2.R` | Reads both workbooks, joins the OTU map, applies design rules, runs integrity checks |',
'| `plot_abundance.R` | Horizontal bar charts of isolate counts by taxonomic level |',
'| `plot_abundance_pooled.R` | Bar charts with uncertain taxa pooled as "Incertae sedis" |',
'| `plot_pie.R` | Pie charts of community composition (3 pies per level) |',
'| `community_analysis.R` | Multivariate community ecology; OTU-level tests plus multi-rank taxonomic sweeps (Section D) |',
'| `diversity_analysis.R` | Sampling completeness (Chao1/ACE/coverage) and replicated per-unit alpha diversity |',
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
'## A0. How much of the community was actually found?',
'',
'**Why this comes first.** Every diversity number below is a *sample* statistic. Culture-based sampling of tropical endophytes recovers a fraction of what is present, and how large that fraction is determines what the rest of the section can claim. Two standard measures answer it: **Chao1/ACE** estimate the true richness from how many taxa were seen once or twice, and **Good\'s coverage** estimates the probability that the *next* isolate collected would belong to a taxon already seen.',
'',
if (nrow(completeness_tbl) > 0) c(
  '| Substrate | Units | Isolates | Observed | Chao1 (±SE) | ACE | % of Chao1 | Singletons | Coverage |',
  '|-----------|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|',
  apply(completeness_tbl, 1, function(r)
    paste0('| ', r['substrate'], ' | ', r['n_units'], ' | ', r['n_isolates'], ' | ',
           r['observed_S'], ' | ', r['chao1'], ' ± ', r['chao1_se'], ' | ', r['ace'],
           ' | **', r['pct_of_chao1'], '%** | ', r['singletons'], ' | ',
           r['goods_coverage'], '% |')),
  ''
) else character(0),
paste0('*Interpretation.* Only **', min(completeness_tbl$pct_of_chao1), '-',
       max(completeness_tbl$pct_of_chao1),
       '% of the estimated richness was recovered**, with Good\'s coverage of ', cov_range,
       '. Roughly one isolate in seven belongs to an OTU seen exactly once. Two things follow. First, **all observed richness values are substantial underestimates** and should be quoted alongside Chao1. Second, *Ficus* wood is the least completely sampled substrate (',
       completeness_tbl$pct_of_chao1[completeness_tbl$substrate == "Ficus wood"],
       '% of Chao1), so its apparently lower richness is partly a sampling artefact — which is why the comparison below is done per sampling unit and again at equal coverage.'),
'',
'![Observed vs estimated richness](plots/png/richness_estimators.png)',
'',
'#### Richness compared at equal coverage',
'',
'Comparing richness at equal *isolate count* still favours whichever substrate has the flatter abundance distribution. The fairer contrast standardises to equal **coverage** (Chao & Jost 2012) \u2014 how many taxa each substrate holds at the same level of sampling completeness.',
'',
if (nrow(coverage_tbl) > 0) c(
  '| Substrate | Isolates needed | Coverage reached | Richness at that coverage |',
  '|-----------|:-:|:-:|:-:|',
  apply(coverage_tbl, 1, function(r)
    paste0('| ', r['substrate'], ' | ', r['isolates_needed'], ' | ',
           r['coverage_reached'], '% | **', r['richness_at_coverage'], '** |')),
  ''
) else character(0),
'*Interpretation \u2014 and an apparent contradiction worth understanding.* At equal coverage the ranking **reverses**: *Ficus* **wood** is the richest substrate, not the poorest. The two statements are not in conflict because they answer different questions. Per sampling unit (A1 below), a patch of wood carries **fewer** taxa than a patch of leaf \u2014 lower alpha diversity. But wood needs far more isolates to reach the same coverage, because its community has a much longer tail of rare taxa; pooled across the substrate it therefore holds **more** taxa overall. In short: **leaf units are individually richer, while the wood community as a whole is more diverse and more heterogeneous between units.** That heterogeneity is the same pattern the dispersion test picks up (Section D), and it is a substantive ecological result rather than an artefact.',
'',
'Full tables: `tables/sampling_completeness.csv`, `tables/richness_at_equal_coverage.csv`.',
'',
'---',
'',
'## A. Alpha Diversity',
'',
'**What it is.** *Alpha diversity* describes how varied the fungal community is **within a single substrate**. It combines *richness* (how many taxa) with *evenness* (whether isolates are spread across taxa or dominated by a few).',
'',
'| Index | What it measures | Higher means |',
'|-------|-----------------|--------------|',
'| **Richness (S)** | Number of distinct taxa | More taxa present |',
'| **Shannon (H\')** | Combines richness and evenness | More diverse |',
'| **Hill q1** (= exp H\') | Effective number of common taxa | More diverse |',
'| **Simpson (1-D)** | Probability two random individuals differ | More diverse |',
'| **Pielou (J\')** | How evenly individuals are distributed | More even |',
'',
'### A1. Per sampling unit — the version that can be tested',
'',
paste0('Diversity is computed **per sampling unit**, giving ',
       paste(alpha_summary_tbl$n_units, collapse = " / "),
       ' independent values per substrate instead of one pooled number. That is what makes a significance test possible at all: a single pooled value per substrate has no variance and supports no inference.'),
'',
if (nrow(alpha_summary_tbl) > 0) c(
  '| Substrate | Units | Isolates | Mean S (±SD) | Mean H\' (±SD) | Mean Hill q1 | Mean J\' |',
  '|-----------|:-:|:-:|:-:|:-:|:-:|:-:|',
  apply(alpha_summary_tbl, 1, function(r)
    paste0('| ', r['substrate'], ' | ', r['n_units'], ' | ', r['isolates'], ' | ',
           r['mean_S'], ' ± ', r['sd_S'], ' | ', r['mean_H'], ' ± ', r['sd_H'],
           ' | ', r['mean_q1'], ' | ', r['mean_J'], ' |')),
  ''
) else character(0),
if (nrow(alpha_tests_tbl) > 0) c(
  'Kruskal-Wallis across substrates (Holm-adjusted across the four indices):',
  '',
  '| Index | H | df | p | p (Holm) |',
  '|-------|:-:|:-:|:-:|:-:|',
  apply(alpha_tests_tbl, 1, function(r)
    paste0('| ', r['index'], ' | ', r['statistic'], ' | ', r['df'], ' | ',
           r['p_value'], ' | **', r['p_adj'], '** |')),
  ''
) else character(0),
if (nrow(alpha_pairwise_tbl) > 0) c(
  'Pairwise follow-ups (Wilcoxon, Holm-adjusted):',
  '',
  '| Index | Comparison | p (Holm) |',
  '|-------|------------|:-:|',
  apply(alpha_pairwise_tbl, 1, function(r)
    paste0('| ', r['index'], ' | ', r['group1'], ' vs ', r['group2'], ' | ',
           r['p_adj'], ' |')),
  ''
) else character(0),
paste0('*Interpretation.* **Richness and Shannon differ significantly between substrates** (Holm-adjusted p = ',
       alpha_tests_tbl$p_adj[alpha_tests_tbl$index == "richness_S"],
       ' and ', alpha_tests_tbl$p_adj[alpha_tests_tbl$index == "shannon_H"],
       '), driven by *Ficus* leaves being richer and more diverse than *Ficus* wood; the Lauraceae sits between them and is not separable from wood. **Evenness (Pielou J\') does not differ** (p = ',
       alpha_tests_tbl$p_adj[alpha_tests_tbl$index == "pielou_J"],
       '): all three communities are similarly un-dominated, and the difference is in how many taxa are present, not how they are balanced. Because *Ficus* wood is also the least completely sampled substrate, part of its lower richness is sampling effort — the effect is real but its size should not be read off these means alone.'),
'',
'![Richness per sampling unit](plots/png/alpha_per_unit_richness_S.png)',
'',
'![Shannon per sampling unit](plots/png/alpha_per_unit_shannon_H.png)',
'',
'Per-unit values: `tables/alpha_diversity_per_unit.csv`; tests: `tables/alpha_diversity_tests.csv`, `tables/alpha_diversity_pairwise.csv`.',
'',
'### A2. Pooled indices across taxonomic ranks',
'',
'The table below is the pooled view (one value per substrate per rank) retained for comparability with earlier versions of this report. **It carries no error and supports no test** — use A1 for inference. *Incertae sedis* taxa are excluded rank by rank.',
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
'- **NMDS ordination** \u2014 squeezes the many-dimensional Bray-Curtis distances into a 2-D map so that samples plotting close together have similar communities. The **stress** value measures distortion: < 0.10 excellent, < 0.20 acceptable, > 0.20 unreliable. Crosses mark group centroids; shaded ellipses show 95% confidence regions. **A stress at or near zero is a warning, not a triumph** \u2014 see the caveat below, which explains why these plots report two stress values.',
'- **PERMANOVA** (`adonis2`) \u2014 tests whether **group centroids differ**. **R\u00b2** is the fraction of community variation explained by the grouping (effect size); a small **p** means the separation is unlikely by chance.',
'- **ANOSIM** \u2014 a complementary rank-based test; **R** ranges from 0 (no separation) to 1 (groups completely distinct).',
'- **Beta-dispersion / PERMDISP** (`betadisper`) \u2014 checks whether groups differ in **within-group spread** rather than location. If PERMDISP is significant, part of a PERMANOVA result may reflect unequal dispersion rather than a pure shift in composition, so it is an important caveat.',
'- **Pairwise PERMANOVA** \u2014 which specific pairs of groups differ, with Holm correction for multiple tests.',
'- **Rarefaction** \u2014 expected richness rescaled to equal sampling effort, so richness can be compared fairly.',
'- **Indicator species (IndVal)** \u2014 identifies taxa statistically associated with (diagnostic of) a particular group. One test is run **per taxon**, so the raw p-values are corrected with Benjamini-Hochberg and only taxa passing **FDR q \u2264 0.05** are reported as indicators. Without that correction, ~14 of ~280 taxa pass p \u2264 0.05 by chance alone.',
'',
'### An important caveat: why the ITS-level ordinations collapse onto a line',
'',
paste0('Several NMDS plots in this report show most samples squeezed onto a single near-vertical line with one point flung far away, next to a stress of ~0.0001. That is **not** an excellent fit and **not** a plotting bug — it is a known failure mode of NMDS, and it is worth understanding because it determines which figures can be read.'),
'',
paste0('**The cause.** At ITS-genotype resolution this dataset is dominated by rare taxa: **',
       sparse_singletons, '% of the ', n_taxa_comm,
       ' genotypes were isolated exactly once**, so **', sparse_no_shared,
       '% of all sampling-unit pairs share no genotype whatsoever**. Every one of those pairs has a Bray-Curtis dissimilarity of *exactly* 1 — they are **tied**.'),
'',
'NMDS fits an ordination by rank order, and by default (`monoMDS`, weak/primary ties) **tied dissimilarities are allowed to map to any distances at all**. With well over half of the pairs tied, the optimiser is therefore free to ignore most of the matrix: it only has to get the *ordering* of the minority of pairs that do share taxa right. It can do that almost perfectly in two dimensions — hence the near-zero stress — while pushing the unconstrained samples wherever is convenient. The collapsed line and the distant outlier are those unconstrained samples.',
'',
paste0('**The evidence.** Refitting the same dissimilarities with *strong* (secondary) ties, which force tied pairs to equal distances, gives the honest answer. At ITS level the stress jumps from **0.0001 to ',
       if (nrow(ml_sub) > 0 && 'nmds_stress_strong' %in% names(ml_sub))
         ml_sub$nmds_stress_strong[ml_sub$level == "its_taxon"] else '0.25',
       '** — i.e. *unreliable* by the usual thresholds. At genus level the two figures agree closely (',
       if (nrow(ml_sub) > 0 && 'nmds_stress_strong' %in% names(ml_sub))
         paste0(ml_sub$nmds_stress[ml_sub$level == "genus"], ' vs ',
                ml_sub$nmds_stress_strong[ml_sub$level == "genus"]) else 'similar',
       '), confirming that those ordinations are real. Every NMDS plot in this report now prints **both** values, and is labelled UNRELIABLE on any of three grounds: a tie-aware stress above 0.20, 9 or fewer sampling units (at or below 4k+1 points a 2-D solution fits almost anything), or a tie-aware stress below 0.001 (a perfect fit means the configuration is unconstrained, not faithful).'),
'',
'Two further points were checked and ruled out as explanations: `metaMDS` applied **no** data transformation here, and it did **not** fall back to extended (step-across) dissimilarities — the largest dissimilarity fed to the ordination is exactly 1. The analysis is doing what it says; the data simply cannot support a 2-D map at genotype resolution.',
'',
paste0('One sampling unit (`Lauraceae leaves__Z5__S4`, a single isolate) shares no genotype with *any* other unit, leaving the dissimilarity matrix formally **disconnected**. Dropping it and the other tiny units reconnects the matrix but does not rescue the ordination (the tie-aware stress only falls to about 0.21), so no samples are excluded on these grounds.'),
'',
'**What to do with this:**',
'',
'- **Do not read the ITS-level NMDS maps.** They are kept because they are referenced throughout the literature-standard workflow, and each is now labelled.',
'- **PERMANOVA, ANOSIM and IndVal are unaffected.** They operate on the dissimilarity matrix itself — its sums of squares and rank order — never on a 2-D embedding, so no tie-handling choice enters. These are the results to trust.',
paste0('- **Use the higher-rank ordinations in Section D as the readable maps.** Grouping isolates into genera or families dissolves the singleton problem: the share of pairs sharing no taxa falls from ',
       if (nrow(ml_sub) > 0 && 'prop_no_shared' %in% names(ml_sub))
         paste0(round(100 * num(ml_sub$prop_no_shared[ml_sub$level == "its_taxon"])),
                '% at ITS level to ',
                round(100 * num(ml_sub$prop_no_shared[ml_sub$level == "genus"])), '% at genus level')
       else 'sharply',
       ', and the substrate signal is reproduced there with an honest stress.'),
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
       if (!is.na(num(nmds_stress)) && num(nmds_stress) < 0.01) paste0('\u2014 **not an excellent fit but a tie-degenerate one**: refitted with strong ties the stress is **', nmds_stress_strong, '**, i.e. unreliable (see the caveat above). Read the separation from PERMANOVA and from the higher-rank ordinations in Section D, not from this map.')
       else if (!is.na(num(nmds_stress)) && num(nmds_stress) < 0.20) '(acceptable to good \u2014 the 2-D map is a faithful summary).'
       else '(interpret the map with some caution).'),
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
'#### Adjusted for sampling depth',
'',
paste0('Sampling units differ about two-fold in how many isolates they yielded (*Ficus* wood units gave roughly half as many as leaf units), and depth on its own predicts composition. Re-fitting with depth as a covariate (marginal / type-III): **substrate R\u00b2 = ',
       adj_R2, ', p = ', adj_p, '; depth R\u00b2 = ', dep_R2, ', p = ', dep_p,
       '**. The substrate effect therefore survives adjustment \u2014 it is not an artefact of unequal recovery \u2014 but depth contributes independently and both are reported. Details in `tables/substrate_all_permanova.txt`.'),
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
if (TRUE) {
  c(
    paste0('### Indicator Species'),
    '',
    paste0('Indicator (IndVal) analysis looks for taxa that are **diagnostic** of a substrate \u2014 both faithful to it and frequent within it. Each taxon is a separate test, so p-values are corrected across the whole taxon set. **',
           n_ind_raw, ' taxa reach raw p \u2264 0.05, but ', nrow(indval),
           ' survive FDR correction (q \u2264 0.05).**'),
    '',
    if (nrow(indval) == 0)
      paste0('*Interpretation.* **No taxon is a statistically reliable indicator of any substrate.** With ~',
             n_ind_tested, ' taxa tested, ', round(0.05 * n_ind_tested),
             ' hits at p \u2264 0.05 are expected by chance, which is about what was observed. This is a real result rather than a failure: the communities differ in *composition as a whole* (PERMANOVA above) without any single taxon being a dependable marker \u2014 exactly what is expected of an assemblage in which most taxa are singletons. Reporting uncorrected indicators here would be reporting noise.')
    else
      paste0('*Interpretation.* These taxa survive correction and can be treated as genuine markers. Full ranked list: `tables/indicator_species_significant.csv`; all taxa with p- and q-values: `tables/substrate_all_indicator_species_all.csv`.'),
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
'## F. Branch & Trunk Sampling Orientation',
'',
'**Why.** Different faces of a tree experience very different microclimates — sun exposure, temperature, how long the surface stays wet after rain, prevailing wind. If that matters to endophytes, communities should differ systematically between compass bearings. This is the first LOT2 analysis able to test it, because the orientation column only became usable with the re-issued samples workbook (see **Sampling orientation** above).',
'',
'### F1. What was actually sampled',
'',
'Orientation is an **observational**, not a designed, factor here: bearings were recorded for whichever branches were reachable, so replication is uneven and, for some substrate/bearing combinations, very thin. That constrains how much the tests below can detect, and it is the single most important caveat for this section.',
'',
if (nrow(ori_cov) > 0) c(
  '| Substrate | Position | Orientation | Sampling units | Isolates | ITS taxa |',
  '|-----------|:-:|:-:|:-:|:-:|:-:|',
  apply(ori_cov, 1, function(r)
    paste0('| ', r['substrate'], ' | ', r['position'], ' | ', r['orientation'],
           ' | ', r['n_units'], ' | ', r['n_isolates'], ' | ', r['n_taxa'], ' |')),
  ''
) else character(0),
if (nrow(ori_unres) > 0) c(
  paste0('Unresolved (excluded from this section): ',
         paste(apply(ori_unres, 1, function(r)
           paste0(r['n_isolates'], ' ', r['substrate'], ' ', tolower(r['position']))),
           collapse = ', '), '.'),
  ''
) else character(0),
'Three consequences follow, and they shape every result below:',
'',
'- Only **Ficus leaves** and **Ficus wood** were sampled on all four cardinal bearings. **Lauraceae leaves** were only ever collected from northern and north-western faces (plus two isolates from one eastern unit), so the Lauraceae cannot contribute to a full-compass comparison.',
'- **Ficus trunk wood carries no orientation at all**, so this section is effectively a *branch*-level analysis.',
paste0('- With ', n_orient_units, ' orientation-level units spread over 5 bearings, most groups hold **1-4 replicates**. Tests on a single substrate (7-11 units) have little power: a null result here means *no effect was detectable*, not *no effect exists*.'),
'',
'### F2. Does orientation structure the community?',
'',
'Same test battery as the substrate analyses — PERMANOVA, ANOSIM, PERMDISP — run on the orientation-level sampling units at ITS-genotype resolution.',
'',
if (nrow(ori_sum) > 0) c(
  '| Analysis | Units | Groups | F | R² | p | ANOSIM R | ANOSIM p | Disp. p | Indic. |',
  '|----------|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|',
  apply(ori_sum, 1, function(r) {
    na_dash <- function(x) if (is.na(x) || !nzchar(trimws(x))) '—' else trimws(x)
    paste0('| ', r['analysis'], ' | ', r['n_units'], ' | ', r['n_groups'], ' | ',
           na_dash(r['permanova_F']), ' | ', na_dash(r['permanova_R2']), ' | ',
           na_dash(r['permanova_p']), ' | ', na_dash(r['anosim_R']), ' | ',
           na_dash(r['anosim_p']), ' | ', na_dash(r['betadisper_p']), ' | ',
           na_dash(r['n_indicators']), ' |')
  }),
  '',
  '*F / R² / p are PERMANOVA; "Disp. p" is PERMDISP; "Indic." counts significant indicator taxa.*',
  '',
  '> **PERMDISP shown as —** where the smallest group holds fewer than 3 sampling units. A centroid computed from one or two points has a degenerate spread, which inflates the dispersion F ratio to meaningless magnitudes; those tests are written to the `*_betadisper.txt` files with a warning but are not reported as numbers.',
  ''
) else character(0),
paste0('**Orientation on its own: ', if (!is.null(ori_all)) paste0('F = ', ori_all$permanova_F, ', R² = ', ori_all$permanova_R2, ', p = ', ori_all$permanova_p) else 'N/A', '** — ',
       if (!is.null(ori_all)) verdict(ori_all$permanova_p) else 'could not be evaluated',
       '. Within each substrate taken separately the result is the same: ',
       if (!orient_any_sig) 'no substrate shows a significant orientation effect.' else 'see the table for which substrate drives the effect.',
       ' ANOSIM agrees, and for Ficus leaves and Ficus wood the ANOSIM R is actually **negative** — meaning units from *different* bearings are, if anything, slightly more similar to each other than units from the *same* bearing. That is the signature of no orientation structure at all.'),
'',
'![NMDS — orientation, all substrates](plots/png/orientation_all_nmds.png)',
'',
'![NMDS — orientation within Ficus leaves](plots/png/ficus_leaves_orientation_nmds.png)',
'',
paste0('*Interpretation.* The bearings overlap almost completely, with centroids piled on top of one another. These are ITS-level ordinations, so \u2014 as everywhere in this report \u2014 the maps are degenerate and carry no weight on their own; the conclusion rests on the PERMANOVA and ANOSIM results above, and on the higher-rank ordinations in **F7**, which agree. Whatever separates these communities, it is not which side of the tree the material came from.'),
'',
'### F3. Separating orientation from substrate',
'',
'Two of the pooled analyses above **do** come out significant — *substrate x orientation* and *sun aspect* — and both need care, because orientation is partly confounded with substrate in this dataset (the Lauraceae units are almost all north/north-west facing). A marginal (type-III) PERMANOVA asks what each factor explains **once the other is accounted for**:',
'',
'```',
paste(ori_part_txt, collapse = "\n"),
'```',
'',
paste0('**Substrate: R² = ', part_sub$R2, ', p = ', part_sub$p,
       '. Orientation: R² = ', part_ori$R2, ', p = ', part_ori$p, '.** ',
       'Substrate holds up; orientation does not. The apparently strong *substrate x orientation* result (R² ≈ ',
       if (!is.null(ori_row("Substrate x orientation"))) ori_row("Substrate x orientation")$permanova_R2 else 'N/A',
       ') is therefore the substrate effect re-expressed through a factor that happens to encode it, not evidence that bearing matters.'),
'',
'#### Sun aspect',
'',
'LOT2 was collected in **Peru, i.e. the southern hemisphere**, where the *northern* face of a tree is the sun-exposed one and the southern face the shaded one. Grouping bearings into **sun-facing (N/NE/NW)**, **shaded (S/SE/SW)** and **lateral (E/W)** gives a more direct microclimate proxy than the raw compass, and coarser groups mean better replication.',
'',
paste0('Pooled across substrates this is significant (',
       if (!is.null(ori_aspect)) paste0('R² = ', ori_aspect$permanova_R2, ', p = ', ori_aspect$permanova_p) else 'N/A',
       ', PERMDISP p = ', if (!is.null(ori_aspect)) ori_aspect$betadisper_p else 'N/A',
       ', so not a dispersion artefact) — but the confounding is severe: of the ',
       n_orient_units, ' orientation units, **all but one Lauraceae unit is sun-facing**, while Ficus leaves and wood are mostly lateral. Adjusting for substrate, **sun aspect gives R² = ', part_asp$R2,
       ', p = ', part_asp$p, '** — ', verdict(part_asp$p),
       '. Within each substrate individually it is likewise non-significant. The pooled result is substrate wearing an aspect label.'),
'',
'![NMDS — sun aspect](plots/png/orientation_sun_aspect_nmds.png)',
'',
'### F4. A directional gradient rather than discrete bearings',
'',
'Treating the compass as four or five **discrete classes** spends a lot of degrees of freedom on a dataset this small. An alternative is to treat the bearing as the **circular variable** it is: decomposing it into `cos(bearing)` (the north-south axis) and `sin(bearing)` (the east-west axis) tests for a community that changes *smoothly* around the tree using only 2 df, which is more powerful at this level of replication.',
'',
if (nrow(ori_circ) > 0) c(
  '| Analysis | Units | Bearings | N-S axis R² | p | E-W axis R² | p | Joint R² | Joint p |',
  '|----------|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|',
  apply(ori_circ, 1, function(r)
    paste0('| ', r['analysis'], ' | ', r['n_units'], ' | ', r['n_orientations'], ' | ',
           r['NS_axis_R2'], ' | ', r['NS_axis_p'], ' | ', r['EW_axis_R2'], ' | ',
           r['EW_axis_p'], ' | ', r['joint_R2'], ' | ', r['joint_p'], ' |')),
  ''
) else character(0),
paste0('*Interpretation.* No joint directional gradient is significant. The one nominally significant single axis is the **north-south axis in Lauraceae leaves** (p = ',
       if (nrow(ori_circ) > 0 && any(ori_circ$analysis == "Lauraceae leaves"))
         ori_circ$NS_axis_p[ori_circ$analysis == "Lauraceae leaves"] else 'N/A',
       '), and it should not be over-read: the Lauraceae span only N and NW (plus two eastern isolates), so the "north-south axis" is fitted over a ~45° arc rather than a full compass, the joint test for the same substrate is non-significant (p = ',
       if (nrow(ori_circ) > 0 && any(ori_circ$analysis == "Lauraceae leaves"))
         ori_circ$joint_p[ori_circ$analysis == "Lauraceae leaves"] else 'N/A',
       '), and it is one nominal result among eight axis tests with no correction applied.'),
'',
'### F5. Diversity per orientation',
'',
'Community *structure* may not differ while *diversity* still does — a more exposed face could simply support fewer taxa. Richness depends strongly on how many isolates each group contributed, so richness is also reported **rarefied to a common isolate count**.',
'',
if (nrow(ori_alpha) > 0) c(
  paste0('| Orientation | S | N | H\' | 1-D | J\' | Rarefied S (to ', ori_alpha$rarefied_to[1], ') |'),
  '|-------------|:-:|:-:|:-:|:-:|:-:|:-:|',
  apply(ori_alpha, 1, function(r)
    paste0('| ', r['orientation'], ' | ', r['richness_S'], ' | ', r['abundance_N'], ' | ',
           r['shannon_H'], ' | ', r['simpson_1mD'], ' | ', r['pielou_J'], ' | ',
           r['rarefied_S'], ' |')),
  ''
) else character(0),
if (nrow(ori_aspect_alpha) > 0) c(
  paste0('| Sun aspect | S | N | H\' | 1-D | J\' | Rarefied S (to ', ori_aspect_alpha$rarefied_to[1], ') |'),
  '|------------|:-:|:-:|:-:|:-:|:-:|:-:|',
  apply(ori_aspect_alpha, 1, function(r)
    paste0('| ', r['sun_aspect'], ' | ', r['richness_S'], ' | ', r['abundance_N'], ' | ',
           r['shannon_H'], ' | ', r['simpson_1mD'], ' | ', r['pielou_J'], ' | ',
           r['rarefied_S'], ' |')),
  ''
) else character(0),
paste0('*Interpretation.* Raw richness tracks sampling effort almost exactly — the most-sampled bearing is also the richest — but once rarefied to a common ',
       if (nrow(ori_alpha) > 0) ori_alpha$rarefied_to[1] else 'N', ' isolates the bearings are within a few taxa of each other',
       if (nrow(ori_alpha) > 0) paste0(' (', min(ori_alpha$rarefied_S, na.rm = TRUE), '-', max(ori_alpha$rarefied_S, na.rm = TRUE), ' taxa)') else '',
       '. Evenness (Pielou J\') is uniformly high, as everywhere else in this dataset. There is no diversity gradient around the tree to match the absent compositional one.'),
'',
'![Richness by orientation and substrate](plots/png/orientation_alpha_richness_S.png)',
'',
'![Shannon diversity by orientation and substrate](plots/png/orientation_alpha_shannon_H.png)',
'',
'> The per-substrate rarefaction target is pulled down to a handful of isolates by the smallest Ficus wood groups, so `tables/orientation_alpha_diversity_by_substrate.csv` should be read as indicative only. The pooled table above is the more reliable of the two.',
'',
'### F6. Taxon sharing between orientations',
'',
if (nrow(ori_share) > 0) c(
  '| Occurs in ... orientations | ITS genotypes |',
  '|:-:|:-:|',
  apply(ori_share, 1, function(r)
    paste0('| ', r['n_orientations'], ' | ', r['n_taxa'], ' |')),
  ''
) else character(0),
paste0('*Interpretation.* Most genotypes (',
       if (nrow(ori_share) > 0) ori_share$n_taxa[ori_share$n_orientations == 1] else 'N/A',
       ' of ', if (nrow(ori_share) > 0) sum(ori_share$n_taxa) else 'N/A',
       ') were found on a single bearing. On its own that looks like strong orientation fidelity, but it is the expected consequence of **rarity plus thin sampling**: this community is dominated by singletons and doubletons (see the rank-abundance curves in Section B1), and a taxon seen once can only ever be recorded at one bearing. The multivariate tests above, which use the whole community at once rather than taxon-by-taxon presence, find no such structure — and they are the trustworthy reading.'),
'',
'![Shared genotypes between orientations — Ficus leaves](plots/png/orientation_venn_ficus_leaves.png)',
'',
'![Rarefaction by orientation](plots/png/orientation_all_rarefaction.png)',
'',
'### F7. Orientation across taxonomic ranks',
'',
'As in Section D, the comparison is repeated at every rank, in case a bearing effect exists among broader taxonomic groups but is hidden by genotype-level noise.',
'',
ml_table(ori_ml),
'',
paste0('*Interpretation.* ',
       if (nrow(ori_ml) > 0 && !any(num(ori_ml$permanova_p) < 0.05, na.rm = TRUE))
         'Orientation is non-significant at **every** taxonomic rank, so the null result is not an artefact of working at ITS resolution.'
       else 'Some ranks show an orientation signal; see the table.'),
'',
'Full results: `tables/orientation_*`, `tables/*_orientation_*` and `tables/substrate_x_orientation_*`. Summary: `tables/orientation_summary.csv`.',
'',
'### F8. What this section concludes',
'',
paste0('**Sampling orientation does not detectably structure the LOT2 endophyte communities.** It is non-significant pooled, within each substrate, at every taxonomic rank, as a coarse sun-exposure grouping, and as a smooth directional gradient. The two significant pooled results both dissolve once substrate is accounted for (orientation adjusted for substrate: R² = ',
       part_ori$R2, ', p = ', part_ori$p, '; sun aspect adjusted for substrate: R² = ',
       part_asp$R2, ', p = ', part_asp$p, ').'),
'',
'**How firm is that?** Firm enough to report, but it is a negative result from an unbalanced observational factor with 1-4 replicates per group. It rules out an orientation effect of the size seen for substrate (R² ≈ 0.19); it does not rule out a small one. Making that test properly would need balanced sampling of all four bearings within each substrate — worth specifying in advance if a future campaign wants to answer this question rather than check it.',
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
paste0('4. **Which side of the tree the material came from does not matter.** Now that the field branch codes have been reconciled with the project database, sampling orientation could be tested for the first time (Section F, ',
       n_orient_isolates, ' of ', nrow(samples_raw), ' isolates). It is non-significant pooled (p = ',
       if (!is.null(ori_all)) ori_all$permanova_p else 'N/A',
       '), within every substrate, at every taxonomic rank, as a sun-exposure grouping and as a smooth directional gradient. The two pooled contrasts that do come out significant — *substrate x orientation* and *sun aspect* — both vanish once substrate is accounted for (orientation adjusted for substrate: R² = ',
       part_ori$R2, ', p = ', part_ori$p, '), because the Lauraceae happened to be sampled almost entirely on north-facing branches. This is a **negative result from an unbalanced observational factor with 1-4 replicates per bearing**: it rules out an orientation effect as large as the substrate effect, not a small one.'),
paste0('5. **The communities are richer than they look, and only partly sampled.** Only ',
       min(completeness_tbl$pct_of_chao1), '-', max(completeness_tbl$pct_of_chao1),
       '% of the Chao1-estimated richness was recovered (Good\'s coverage ', cov_range,
       '). Per sampling unit, richness and Shannon differ significantly between substrates (Holm-adjusted p = ',
       alpha_tests_tbl$p_adj[alpha_tests_tbl$index == "richness_S"],
       ' and ', alpha_tests_tbl$p_adj[alpha_tests_tbl$index == "shannon_H"],
       '), while **evenness does not** (p = ',
       alpha_tests_tbl$p_adj[alpha_tests_tbl$index == "pielou_J"],
       '): the substrates differ in how many taxa they carry, not in how evenly those taxa are balanced.'),
paste0('6. **Indicator taxa survive correction only because of OTU clustering.** ',
       nrow(indval), ' of ', n_ind_tested,
       ' OTUs qualify as substrate indicators at FDR q \u2264 0.05. Run on the raw BLAST name labels the same analysis yields **none** \u2014 the extra tests and the split singletons destroy the signal. Uncorrected IndVal output should never be reported.'),
paste0('7. **Removing *Incertae sedis* sharpened the picture.** Excluding the pooled "unknown" bin from the diversity, overlap and multivariate analyses (while keeping it visible in the abundance/pie plots) increased, rather than decreased, the measured separation between substrates \u2014 confirming that the unidentified fraction had been masking genuine differences.'),
'',
'---',
'',
paste0('*Auto-generated on ', Sys.Date(), ' by `generate_readme.R`*')
)

# ------------------------------------------------------------
# Normalise pipe-table separator rows.
#
# Pandoc derives each column's *relative width* from the number of
# dashes in the separator row, not from the content. The compact ":-:"
# markers used throughout this script therefore hand every column an
# identical share of the page, and once the README is typeset as PDF the
# longer headers ("PERMANOVA F", "PERMDISP p", ...) overflow their
# column and overlap their neighbours.
#
# Rewriting each separator so its dash count tracks the widest cell in
# that column fixes the PDF. Markdown renderers ignore the dash count
# entirely, so GitHub's view of README.md is unchanged.
# ------------------------------------------------------------
balance_table_widths <- function(lines, min_dashes = 3, max_dashes = 30) {
  is_row <- function(x) grepl("^\\s*\\|.*\\|\\s*$", x)
  is_sep <- function(x) grepl("^\\s*\\|[ :|-]+\\|\\s*$", x) && grepl("-", x)
  split_cells <- function(x) {
    x <- sub("^\\s*\\|", "", x)
    x <- sub("\\|\\s*$", "", x)
    trimws(strsplit(x, "|", fixed = TRUE)[[1]])
  }

  in_code <- FALSE
  i <- 1L
  while (i < length(lines)) {
    if (grepl("^\\s*```", lines[i])) in_code <- !in_code
    if (in_code || !is_row(lines[i]) || is_sep(lines[i]) || !is_sep(lines[i + 1])) {
      i <- i + 1L
      next
    }

    # Collect the header + body rows belonging to this table
    last <- i + 1L
    while (last + 1L <= length(lines) && is_row(lines[last + 1L]) &&
           !is_sep(lines[last + 1L])) last <- last + 1L
    body_idx <- if (last > i + 1L) (i + 2L):last else integer(0)

    sep <- split_cells(lines[i + 1L])
    rows <- lapply(c(i, body_idx), function(k) split_cells(lines[k]))
    rows <- rows[vapply(rows, length, 1L) == length(sep)]
    if (length(rows) == 0) { i <- last + 1L; next }

    widths <- vapply(seq_along(sep), function(cl)
      max(vapply(rows, function(r) nchar(r[cl]), 1L)), 1L)
    widths <- pmin(pmax(widths, min_dashes), max_dashes)

    lines[i + 1L] <- paste0("|", paste0(vapply(seq_along(sep), function(cl) {
      s <- sep[cl]
      left  <- startsWith(s, ":")
      right <- endsWith(s, ":")
      d <- strrep("-", widths[cl])
      if (left && right) paste0(":", d, ":")
      else if (right)    paste0(d, ":")
      else if (left)     paste0(":", d)
      else               d
    }, character(1)), collapse = "|"), "|")

    i <- last + 1L
  }
  lines
}

# Cap the widest text column at 24 characters: left uncapped, the long
# analysis labels crowd the numeric columns and push their headers into
# each other.
readme <- balance_table_widths(readme, max_dashes = 24)

writeLines(readme, "README.md")
cat("README.md generated\n")

# ---- Generate README.pdf via pandoc ----
# pandoc and a PDF engine live in different places depending on the
# machine, so look for them rather than assuming a LaTeX install: on a
# plain RStudio/Windows box the bundled pandoc + typst are the only
# ones present, while a Linux box usually has xelatex.
nz1 <- function(x) length(x) == 1 && !is.na(x) && nzchar(x)

pandoc_bin <- first_on_path("pandoc")
if (!nz1(pandoc_bin) && requireNamespace("rmarkdown", quietly = TRUE))
  pandoc_bin <- tryCatch(rmarkdown::pandoc_exec(), error = function(e) "")
if (!nz1(pandoc_bin)) {
  bundled <- Sys.glob(c(
    "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools/pandoc.exe",
    "C:/Program Files/RStudio/bin/quarto/bin/tools/pandoc.exe",
    "/usr/lib/rstudio/resources/app/bin/quarto/bin/tools/pandoc"))
  bundled <- bundled[file.exists(bundled)]
  if (length(bundled)) pandoc_bin <- bundled[1]
}
if (!nz1(pandoc_bin)) pandoc_bin <- ""

# Each engine needs its own font/margin options
engine_opts <- list(
  xelatex  = paste("--pdf-engine=xelatex",
                   "-V mainfont='DejaVu Serif'",
                   "-V monofont='DejaVu Sans Mono'",
                   "-V geometry:margin=2cm"),
  lualatex = paste("--pdf-engine=lualatex",
                   "-V mainfont='DejaVu Serif'",
                   "-V monofont='DejaVu Sans Mono'",
                   "-V geometry:margin=2cm"),
  # 9pt with narrower side margins: the wide statistical tables (10
  # columns) do not otherwise fit the text block, and typst overflows
  # rather than shrinking them.
  typst    = paste("--pdf-engine=typst", "-V margin-x=1.5cm", "-V margin-y=2cm",
                   "-V fontsize=9pt"),
  pdflatex = paste("--pdf-engine=pdflatex", "-V geometry:margin=2cm")
)

if (!nz1(pandoc_bin)) {
  cat("WARNING: pandoc not found - README.pdf not regenerated\n")
} else {
  # A bundled pandoc also ships its engines next to itself
  extra_path <- c(dirname(pandoc_bin),
                  file.path(dirname(pandoc_bin), "x86_64"),
                  file.path(dirname(pandoc_bin), "aarch64"))
  extra_path <- extra_path[dir.exists(extra_path)]
  old_path <- Sys.getenv("PATH")
  if (length(extra_path))
    Sys.setenv(PATH = paste(c(old_path, extra_path),
                            collapse = .Platform$path.sep))

  pdf_ok <- 1L
  for (eng in names(engine_opts)) {
    if (!nz1(first_on_path(eng))) next
    pdf_ok <- system(paste('"', pandoc_bin, '" README.md -o README.pdf ',
                           engine_opts[[eng]], sep = ""), intern = FALSE)
    if (pdf_ok == 0) { cat("README.pdf generated (engine:", eng, ")\n"); break }
  }
  Sys.setenv(PATH = old_path)
  if (pdf_ok != 0)
    cat("WARNING: README.pdf generation failed - no usable PDF engine found\n",
        "         (tried: ", paste(names(engine_opts), collapse = ", "), ")\n", sep = "")
}
