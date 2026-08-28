library(vegan)
library(indicspecies)
library(ggVennDiagram)

dir.create("plots/community",   recursive = TRUE, showWarnings = FALSE)
dir.create("tables",            recursive = TRUE, showWarnings = FALSE)

tax_hierarchy <- c("phylum", "subphylum", "superclass", "class",
                   "subclass", "order", "family", "genus", "species")
uncertain_values <- c("?", "NO", "incertae sedis")
incertae_label   <- "Incertae sedis"

dat_agg <- all_isolates %>%
  mutate(across(all_of(tax_hierarchy),
                ~ ifelse(is.na(.) | . %in% uncertain_values, incertae_label, .)),
         its_taxon = ifelse(is.na(its_taxon) | its_taxon %in% uncertain_values,
                            incertae_label, its_taxon),
         across(c(n_fungi_laur_leaf, n_fungi_fic_leaf, n_fungi_fic_wood),
                ~ replace_na(., 0)))

# ============================================================
# A. ALPHA DIVERSITY (aggregate, per substrate)
# ============================================================
cat("\n====== A. Alpha diversity ======\n")

diversity_results <- list()

for (level in c("its_taxon", tax_hierarchy)) {
  if (level == "its_taxon") {
    agg <- dat_agg %>%
      group_by(its_taxon) %>%
      summarise(Lauraceae_leaves = sum(n_fungi_laur_leaf),
                Ficus_leaves     = sum(n_fungi_fic_leaf),
                Ficus_wood       = sum(n_fungi_fic_wood), .groups = "drop")
  } else {
    idx <- which(tax_hierarchy == level)
    agg <- dat_agg %>%
      group_by(across(all_of(tax_hierarchy[1:idx]))) %>%
      summarise(Lauraceae_leaves = sum(n_fungi_laur_leaf),
                Ficus_leaves     = sum(n_fungi_fic_leaf),
                Ficus_wood       = sum(n_fungi_fic_wood), .groups = "drop")
  }

  mat <- as.matrix(agg[, c("Lauraceae_leaves", "Ficus_leaves", "Ficus_wood")])
  comm <- t(mat)  # substrates as rows, taxa as columns

  S     <- specnumber(comm)
  N     <- rowSums(comm)
  H     <- diversity(comm, index = "shannon")
  D     <- diversity(comm, index = "simpson")
  invD  <- diversity(comm, index = "invsimpson")
  J     <- H / log(S)

  div_df <- data.frame(
    taxonomic_level = level,
    substrate       = c("Lauraceae leaves", "Ficus leaves", "Ficus wood"),
    richness_S      = S,
    abundance_N     = N,
    shannon_H       = round(H, 4),
    simpson_1mD     = round(D, 4),
    inv_simpson     = round(invD, 4),
    pielou_J        = round(J, 4),
    row.names       = NULL
  )
  diversity_results[[level]] <- div_df
}

diversity_table <- bind_rows(diversity_results)
write.csv(diversity_table, "tables/alpha_diversity.csv", row.names = FALSE)
cat("Saved: tables/alpha_diversity.csv\n")

# Plot diversity indices
for (idx_name in c("shannon_H", "simpson_1mD", "inv_simpson", "pielou_J")) {
  p <- ggplot(diversity_table,
              aes(x = taxonomic_level, y = .data[[idx_name]],
                  fill = substrate)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = c("Lauraceae leaves" = "#2E86AB",
                                 "Ficus leaves"     = "#A23B72",
                                 "Ficus wood"       = "#F18F01")) +
    scale_x_discrete(limits = c("its_taxon", tax_hierarchy)) +
    labs(title = paste(idx_name, "across substrates"),
         x = "Taxonomic level", y = idx_name, fill = "Substrate") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title  = element_text(face = "bold"),
          legend.position = "top")
  ggsave(file.path("plots/community", paste0("alpha_", idx_name, ".pdf")),
         plot = p, width = 10, height = 6)
}
cat("Saved: alpha diversity plots\n")

# ============================================================
# B1. RANK-ABUNDANCE CURVES (genus level)
# ============================================================
cat("\n====== B1. Rank-abundance curves ======\n")

for (level in c("genus", "species", "its_taxon")) {
  if (level == "its_taxon") {
    agg <- dat_agg %>%
      group_by(its_taxon) %>%
      summarise(Lauraceae_leaves = sum(n_fungi_laur_leaf),
                Ficus_leaves     = sum(n_fungi_fic_leaf),
                Ficus_wood       = sum(n_fungi_fic_wood), .groups = "drop")
  } else {
    idx <- which(tax_hierarchy == level)
    agg <- dat_agg %>%
      group_by(across(all_of(tax_hierarchy[1:idx]))) %>%
      summarise(Lauraceae_leaves = sum(n_fungi_laur_leaf),
                Ficus_leaves     = sum(n_fungi_fic_leaf),
                Ficus_wood       = sum(n_fungi_fic_wood), .groups = "drop")
  }

  ra_list <- list()
  for (sub in c("Lauraceae_leaves", "Ficus_leaves", "Ficus_wood")) {
    counts <- sort(agg[[sub]][agg[[sub]] > 0], decreasing = TRUE)
    ra_list[[sub]] <- data.frame(
      substrate  = gsub("_", " ", sub),
      rank       = seq_along(counts),
      abundance  = counts,
      rel_abund  = counts / sum(counts)
    )
  }
  ra_data <- bind_rows(ra_list)

  p <- ggplot(ra_data, aes(x = rank, y = rel_abund, colour = substrate)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    scale_colour_manual(values = c("Lauraceae leaves" = "#2E86AB",
                                   "Ficus leaves"     = "#A23B72",
                                   "Ficus wood"       = "#F18F01")) +
    scale_y_log10() +
    labs(title = paste("Rank-abundance curve —", level),
         x = "Rank", y = "Relative abundance (log scale)",
         colour = "Substrate") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), legend.position = "top")

  ggsave(file.path("plots/community", paste0("rank_abundance_", level, ".pdf")),
         plot = p, width = 10, height = 6)
}
cat("Saved: rank-abundance curves\n")

# ============================================================
# B2. RELATIVE ABUNDANCE STACKED BAR CHARTS
# ============================================================
cat("\n====== B2. Relative abundance stacked bars ======\n")

for (level in c("phylum", "class", "order", "family", "genus")) {
  idx <- which(tax_hierarchy == level)
  agg <- dat_agg %>%
    group_by(across(all_of(tax_hierarchy[1:idx]))) %>%
    summarise(Lauraceae_leaves = sum(n_fungi_laur_leaf),
              Ficus_leaves     = sum(n_fungi_fic_leaf),
              Ficus_wood       = sum(n_fungi_fic_wood), .groups = "drop") %>%
    mutate(label = .data[[level]])

  plot_data <- agg %>%
    select(label, Lauraceae_leaves, Ficus_leaves, Ficus_wood) %>%
    group_by(label) %>%
    summarise(across(everything(), sum), .groups = "drop") %>%
    pivot_longer(-label, names_to = "substrate", values_to = "count") %>%
    mutate(substrate = gsub("_", " ", substrate)) %>%
    group_by(substrate) %>%
    mutate(rel = count / sum(count)) %>%
    ungroup()

  # Order taxa by total abundance
  tax_order <- plot_data %>%
    group_by(label) %>% summarise(tot = sum(count), .groups = "drop") %>%
    arrange(desc(tot)) %>% pull(label)
  plot_data$label <- factor(plot_data$label, levels = rev(tax_order))

  n <- length(tax_order)
  pal <- if (n <= 12) {
    scales::hue_pal()(n)
  } else {
    colorRampPalette(c(
      "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628",
      "#F781BF","#999999","#66C2A5","#FC8D62","#8DA0CB","#E78AC3",
      "#A6D854","#FFD92F","#1B9E77","#D95F02","#7570B3","#E7298A"
    ))(n)
  }

  p <- ggplot(plot_data,
              aes(x = substrate, y = rel, fill = label)) +
    geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = setNames(pal, tax_order), name = level) +
    scale_y_continuous(labels = scales::percent) +
    labs(title = paste("Relative abundance by", level),
         x = NULL, y = "Relative abundance") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), legend.position = "right")

  ggsave(file.path("plots/community",
                    paste0("rel_abundance_", level, ".pdf")),
         plot = p, width = 10, height = 7)
}
cat("Saved: relative abundance stacked bar plots\n")

# ============================================================
# B3. VENN DIAGRAMS (genus & species)
# ============================================================
cat("\n====== B3. Venn diagrams ======\n")

for (level in c("genus", "species", "its_taxon")) {
  if (level == "its_taxon") {
    agg <- dat_agg %>%
      group_by(its_taxon) %>%
      summarise(Lauraceae_leaves = sum(n_fungi_laur_leaf),
                Ficus_leaves     = sum(n_fungi_fic_leaf),
                Ficus_wood       = sum(n_fungi_fic_wood), .groups = "drop") %>%
      rename(taxon = its_taxon)
  } else {
    idx <- which(tax_hierarchy == level)
    agg <- dat_agg %>%
      group_by(across(all_of(tax_hierarchy[1:idx]))) %>%
      summarise(Lauraceae_leaves = sum(n_fungi_laur_leaf),
                Ficus_leaves     = sum(n_fungi_fic_leaf),
                Ficus_wood       = sum(n_fungi_fic_wood), .groups = "drop") %>%
      mutate(taxon = .data[[level]])
  }

  taxa_lists <- list(
    `Lauraceae leaves` = agg$taxon[agg$Lauraceae_leaves > 0],
    `Ficus leaves`     = agg$taxon[agg$Ficus_leaves > 0],
    `Ficus wood`       = agg$taxon[agg$Ficus_wood > 0]
  )

  p <- ggVennDiagram(taxa_lists, label_alpha = 0) +
    scale_fill_gradient(low = "#F4FAFE", high = "#2E86AB") +
    labs(title = paste("Shared taxa —", level)) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  ggsave(file.path("plots/community", paste0("venn_", level, ".pdf")),
         plot = p, width = 8, height = 7)
}
cat("Saved: Venn diagrams\n")

# ============================================================
# C. SAMPLE-LEVEL MULTIVARIATE ANALYSES
# ============================================================
cat("\n====== C. Building sample × taxon matrix ======\n")

# Build sample-level community matrix from individual sheets
# Ficus leaves: replicate = tree_orient (S1-S10)
fl <- endo_leaf_ficus %>%
  mutate(substrate = "Ficus leaves", sample_id = tree_orient) %>%
  select(substrate, sample_id, its_taxon)

# Ficus wood: replicate = sample prefix
fw <- endo_wood_ficus %>%
  mutate(substrate = "Ficus wood",
         sample_id = sub("/.*", "", sample_code)) %>%
  select(substrate, sample_id, its_taxon)

# Host wood (Lauraceae): replicate = tree_zone
hw <- endo_wood_host %>%
  mutate(substrate = "Host wood",
         sample_id = paste0("zone_", tree_zone)) %>%
  select(substrate, sample_id, its_taxon)

all_samples <- bind_rows(fl, fw, hw) %>%
  mutate(uid = paste(substrate, sample_id, sep = "__"))

# Community matrix: samples × taxa
comm_wide <- all_samples %>%
  count(uid, its_taxon) %>%
  pivot_wider(names_from = its_taxon, values_from = n, values_fill = 0) %>%
  as.data.frame()
rownames(comm_wide) <- comm_wide$uid
comm_mat <- as.matrix(comm_wide[, -1])

# Metadata
meta <- all_samples %>%
  distinct(uid, substrate, sample_id) %>%
  arrange(match(uid, rownames(comm_mat)))

cat("Community matrix:", nrow(comm_mat), "samples x", ncol(comm_mat), "taxa\n")
cat("Substrates:", paste(table(meta$substrate), names(table(meta$substrate)),
                         sep = "×", collapse = ", "), "\n")

# ---- C1. NMDS ----
cat("\n====== C1. NMDS ======\n")
set.seed(42)
nmds <- metaMDS(comm_mat, distance = "bray", k = 2, trymax = 200)
cat("Stress:", round(nmds$stress, 4), "\n")

nmds_scores <- as.data.frame(scores(nmds, display = "sites"))
nmds_scores$substrate <- meta$substrate
nmds_scores$sample_id <- meta$sample_id

substrate_colours_3 <- c("Ficus leaves" = "#A23B72",
                          "Ficus wood"   = "#F18F01",
                          "Host wood"    = "#2E86AB")

# Compute centroids
centroids <- nmds_scores %>%
  group_by(substrate) %>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2), .groups = "drop")

p_nmds <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, colour = substrate)) +
  stat_ellipse(aes(fill = substrate), geom = "polygon",
               alpha = 0.15, level = 0.95, linetype = 2) +
  geom_point(size = 3) +
  geom_point(data = centroids, shape = 4, size = 5, stroke = 1.5) +
  scale_colour_manual(values = substrate_colours_3) +
  scale_fill_manual(values = substrate_colours_3) +
  labs(title = "NMDS ordination of fungal communities",
       subtitle = paste("Bray-Curtis | Stress =", round(nmds$stress, 3)),
       colour = "Substrate", fill = "Substrate") +
  theme_minimal(base_size = 12) +
  theme(plot.title    = element_text(face = "bold"),
        legend.position = "top")

ggsave("plots/community/nmds_ordination.pdf", plot = p_nmds, width = 9, height = 7)
cat("Saved: NMDS ordination plot\n")

# Also plot with species scores (top taxa)
sp_scores <- as.data.frame(scores(nmds, display = "species"))
sp_scores$taxon <- rownames(sp_scores)
# Show only most abundant taxa
taxon_totals <- colSums(comm_mat)
top_taxa <- names(sort(taxon_totals, decreasing = TRUE))[1:min(15, ncol(comm_mat))]
sp_top <- sp_scores[sp_scores$taxon %in% top_taxa, ]

p_nmds_sp <- p_nmds +
  geom_text(data = sp_top, aes(x = NMDS1, y = NMDS2, label = taxon),
            inherit.aes = FALSE, size = 2.5, alpha = 0.7, fontface = "italic") +
  labs(title = "NMDS with top 15 taxa overlaid")

ggsave("plots/community/nmds_with_species.pdf", plot = p_nmds_sp, width = 11, height = 8)
cat("Saved: NMDS with species overlay\n")

# ---- C2. PERMANOVA ----
cat("\n====== C2. PERMANOVA (adonis2) ======\n")
set.seed(42)
perm <- adonis2(comm_mat ~ substrate, data = meta, method = "bray",
                permutations = 999)
print(perm)

sink("tables/permanova_results.txt")
cat("PERMANOVA — Bray-Curtis distance\n")
cat("Formula: community ~ substrate\n")
cat("Permutations: 999\n\n")
print(perm)
sink()
cat("Saved: tables/permanova_results.txt\n")

# ---- C3. ANOSIM ----
cat("\n====== C3. ANOSIM ======\n")
set.seed(42)
anos <- anosim(comm_mat, grouping = meta$substrate, distance = "bray",
               permutations = 999)
cat("ANOSIM R:", round(anos$statistic, 4),
    "  p-value:", anos$signif, "\n")

sink("tables/anosim_results.txt", append = FALSE)
cat("ANOSIM — Bray-Curtis distance\n")
cat("Permutations: 999\n\n")
cat("R statistic:", round(anos$statistic, 4), "\n")
cat("p-value:", anos$signif, "\n\n")
print(summary(anos))
sink()
cat("Saved: tables/anosim_results.txt\n")

# ---- C4. Beta-dispersion ----
cat("\n====== C4. Beta-dispersion (betadisper) ======\n")
bc_dist <- vegdist(comm_mat, method = "bray")
bd <- betadisper(bc_dist, meta$substrate)
bd_perm <- permutest(bd, permutations = 999)
cat("Betadisper F:", round(bd_perm$tab$F[1], 4),
    "  p-value:", round(bd_perm$tab$`Pr(>F)`[1], 4), "\n")

sink("tables/betadisper_results.txt", append = FALSE)
cat("Beta-dispersion test (betadisper + permutest)\n")
cat("Tests homogeneity of multivariate dispersions\n\n")
print(bd_perm)
cat("\nGroup mean distances to centroid:\n")
print(data.frame(substrate = levels(bd$group),
                 mean_dist = round(tapply(bd$distances, bd$group, mean), 4)))
sink()
cat("Saved: tables/betadisper_results.txt\n")

# ---- C5. Pairwise PERMANOVA ----
cat("\n====== C5. Pairwise PERMANOVA ======\n")
substrates <- unique(meta$substrate)
pairs <- combn(substrates, 2, simplify = FALSE)
pw_results <- list()
for (pr in pairs) {
  sel <- meta$substrate %in% pr
  set.seed(42)
  pw <- adonis2(comm_mat[sel, ] ~ substrate, data = meta[sel, ],
                method = "bray", permutations = 999)
  pw_results[[paste(pr, collapse = " vs ")]] <- data.frame(
    pair     = paste(pr, collapse = " vs "),
    F_stat   = round(pw$F[1], 4),
    R2       = round(pw$R2[1], 4),
    p_value  = pw$`Pr(>F)`[1]
  )
}
pw_table <- bind_rows(pw_results)
pw_table$p_adj <- p.adjust(pw_table$p_value, method = "holm")
print(pw_table)

write.csv(pw_table, "tables/pairwise_permanova.csv", row.names = FALSE)
cat("Saved: tables/pairwise_permanova.csv\n")

# ============================================================
# D1. RAREFACTION CURVES
# ============================================================
cat("\n====== D1. Rarefaction curves ======\n")

# Per-substrate rarefaction (pool samples within substrate)
comm_by_sub <- list()
for (sub in substrates) {
  rows <- meta$substrate == sub
  comm_by_sub[[sub]] <- colSums(comm_mat[rows, , drop = FALSE])
}
comm_sub_mat <- do.call(rbind, comm_by_sub)

min_n <- min(rowSums(comm_sub_mat))
rarefy_data <- list()
for (sub in substrates) {
  n_total <- sum(comm_sub_mat[sub, ])
  steps <- unique(c(seq(1, n_total, length.out = 50), n_total))
  rare <- rarefy(comm_sub_mat[sub, , drop = FALSE],
                 sample = floor(steps))
  rarefy_data[[sub]] <- data.frame(
    substrate = sub,
    n_individuals = floor(steps),
    expected_species = as.numeric(rare)
  )
}
rare_df <- bind_rows(rarefy_data)

p_rare <- ggplot(rare_df, aes(x = n_individuals, y = expected_species,
                               colour = substrate)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = substrate_colours_3) +
  labs(title = "Rarefaction curves",
       x = "Number of individuals", y = "Expected number of taxa (ITS)",
       colour = "Substrate") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")

ggsave("plots/community/rarefaction_curves.pdf", plot = p_rare, width = 9, height = 6)
cat("Saved: rarefaction curves\n")

# ============================================================
# D2. INDICATOR SPECIES ANALYSIS
# ============================================================
cat("\n====== D2. Indicator species analysis ======\n")
set.seed(42)
indval <- multipatt(comm_mat, cluster = meta$substrate,
                    func = "IndVal.g", control = how(nperm = 999))

indval_summary <- capture.output(summary(indval, indvalcomp = TRUE))
writeLines(indval_summary, "tables/indicator_species.txt")
cat("Saved: tables/indicator_species.txt\n")

# Extract significant indicators
sig <- indval$sign[indval$sign$p.value <= 0.05, , drop = FALSE]
if (nrow(sig) > 0) {
  sig$taxon <- rownames(sig)
  sig <- sig %>% arrange(p.value)
  write.csv(sig, "tables/indicator_species_significant.csv", row.names = FALSE)
  cat("Significant indicators (p<=0.05):", nrow(sig), "\n")
  cat("Saved: tables/indicator_species_significant.csv\n")
} else {
  cat("No significant indicator species found (p<=0.05)\n")
}

cat("\n====== All community analyses complete ======\n")
