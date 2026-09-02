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

# Venn counts
venn_genus <- dat_agg %>%
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
'### Handling of uncertain taxonomy',
'',
'Entries with `NA`, `"?"`, `"NO"`, or `"incertae sedis"` are pooled into **"Incertae sedis"** before aggregation.',
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
'| `community_analysis.R` | Full community ecology analyses across all comparisons |',
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
'## A. Alpha Diversity (from pooled counts)',
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
'Full data: `tables/alpha_diversity.csv`',
'',
'![Shannon diversity](plots/png/alpha_shannon_H.png)',
'',
'---',
'',
'## B. Community Composition',
'',
'### B1. Rank-Abundance Curves',
'',
'![Rank-abundance at genus level](plots/png/rank_abundance_genus.png)',
'',
'### B2. Relative Abundance',
'',
'![Relative abundance by phylum](plots/png/rel_abundance_phylum.png)',
'',
'![Relative abundance by genus](plots/png/rel_abundance_genus.png)',
'',
'### B3. Venn Diagrams — Shared Taxa',
'',
paste0('At **genus level** (', vg$total, ' total genera): **', vg$all_three, '** shared across all 3 substrates, **', vg$ll_only, '** unique to Lauraceae leaves, **', vg$fl_only, '** unique to Ficus leaves, **', vg$fw_only, '** unique to Ficus wood.'),
'',
'![Venn diagram — genus](plots/png/venn_genus.png)',
'',
'---',
'',
'# Multivariate Statistical Analyses',
'',
paste0('Sample-level analyses use **', n_samples, ' samples** and **', n_taxa_comm, ' ITS taxa**, with **Bray-Curtis dissimilarity** and **999 permutations**.'),
'',
'---',
'',
'## C1. Substrate Comparison: Ficus Leaves vs Ficus Wood vs Lauraceae Leaves',
'',
'### NMDS Ordination',
'',
paste0('**Stress = ', nmds_stress, '** (< 0.10 = good; < 0.20 = acceptable)'),
'',
'![NMDS — all substrates](plots/png/substrate_all_nmds.png)',
'',
'![NMDS with species overlay](plots/png/substrate_all_nmds_species.png)',
'',
'### PERMANOVA',
'',
'```',
paste(perm_txt, collapse = "\n"),
'```',
'',
paste0('F = ', perm_F, ', R² = ', perm_R2, ', p = ', perm_p),
'',
'### ANOSIM',
'',
paste0('R = ', anos_R, ', p = ', anos_p),
'',
'### Beta-dispersion (PERMDISP)',
'',
paste0('F = ', bd_F, ', p = ', bd_p),
'',
if (nrow(pw) > 0) {
  c(
    '### Pairwise PERMANOVA',
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
    ''
  )
} else character(0),
'### Rarefaction',
'',
'![Rarefaction — all substrates](plots/png/substrate_all_rarefaction.png)',
'',
if (nrow(indval) > 0) {
  c(
    paste0('### Indicator Species (', nrow(indval), ' significant, p ≤ 0.05)'),
    '',
    'Full results: `tables/indicator_species_significant.csv`',
    ''
  )
} else character(0),
'---',
'',
'## C2. Ficus Leaves vs Lauraceae Leaves',
'',
paste0('PERMANOVA: F = ', leaves_perm_F, ', R² = ', leaves_perm_R2, ', p = ', leaves_perm_p),
'',
'![NMDS — Ficus vs Lauraceae leaves](plots/png/substrate_leaves_nmds.png)',
'',
'Full results in `tables/substrate_leaves_*.txt`',
'',
'---',
'',
'## E. Zone-Based Analyses',
'',
'### E1. Ficus Wood Across Zones',
'',
'Community comparison of Ficus wood isolates collected at different tree zones (1-6).',
'',
'![NMDS — Ficus wood zones](plots/png/ficus_wood_zones_nmds.png)',
'',
'Full results in `tables/ficus_wood_zones_*.txt`',
'',
'### E2. Ficus Wood: Trunk (Zones 1-5) vs Branch (Zone 6)',
'',
paste0('PERMANOVA: F = ', fw_tb_F, ', R² = ', fw_tb_R2, ', p = ', fw_tb_p),
'',
'![NMDS — trunk vs branch](plots/png/ficus_wood_trunk_vs_branch_nmds.png)',
'',
'### E4. Substrate Comparison at Zone 6 (Branch Level Only)',
'',
paste0('PERMANOVA: F = ', z6_perm_F, ', R² = ', z6_perm_R2, ', p = ', z6_perm_p),
'',
'![NMDS — substrates at zone 6](plots/png/substrate_zone6_branch_nmds.png)',
'',
'### E5. Substrate × Position Interaction',
'',
'Combined factor analysis testing whether community composition differs by substrate-position combinations.',
'',
'![NMDS — substrate × position](plots/png/substrate_x_position_nmds.png)',
'',
'---',
'',
'## Abundance Plots',
'',
'### `plots/abundance/`',
'',
'Horizontal grouped bar charts showing absolute isolate counts per taxonomic level, broken down by substrate.',
'',
'![Abundance by genus (pooled)](plots/png/abundance_by_genus.png)',
'',
'### `plots/pie_charts/`',
'',
'![Pie chart by phylum](plots/png/pie_by_phylum.png)',
'',
'![Pie chart by genus](plots/png/pie_by_genus.png)',
'',
'---',
'',
paste0('*Auto-generated on ', Sys.Date(), ' by `generate_readme.R`*')
)

writeLines(readme, "README.md")
cat("README.md generated\n")

# ---- Generate README.pdf via pandoc ----
pdf_ok <- system("pandoc README.md -o README.pdf --pdf-engine=xelatex 2>&1",
                 intern = FALSE)
if (pdf_ok == 0) {
  cat("README.pdf generated\n")
} else {
  cat("WARNING: README.pdf generation failed (pandoc/xelatex not available?)\n")
}
