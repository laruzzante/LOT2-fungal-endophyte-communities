# =============================================================
# Side-by-side comparison of the three taxonomic unit sets
# (ITS100 reference, ITS99 in production, OTU97 automatic)
# -> results/comparison/{tables,plots}
# =============================================================
cmp_dir <- "results/comparison"
dir.create(file.path(cmp_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cmp_dir, "plots"),  recursive = TRUE, showWarnings = FALSE)
set_levels <- names(unit_sets)
set_names  <- setNames(paste0(set_levels, " - ", sapply(unit_sets, `[[`, "label")), set_levels)

combine <- function(file) {
  out <- bind_rows(lapply(set_levels, function(u) {
    f <- file.path("results", u, "tables", file)
    if (file.exists(f)) cbind(unit_set = u, read.csv(f, check.names = FALSE, stringsAsFactors = FALSE))
  }))
  if (nrow(out)) {
    out$unit_set <- factor(out$unit_set, levels = set_levels)
    write.csv(out, file.path(cmp_dir, "tables", sub("\\.csv$", "_all_unit_sets.csv", file)),
              row.names = FALSE)
  }
  out
}

# ---- Unit-set overview ----
overview <- bind_rows(lapply(set_levels, function(u) {
  s <- unit_sets[[u]]$samples
  n_per_taxon <- table(s$its_taxon)
  data.frame(unit_set = u, definition = unit_sets[[u]]$label, isolates = nrow(s),
             taxa = length(n_per_taxon),
             singleton_taxa = sum(n_per_taxon == 1),
             pct_singleton_taxa = round(100 * mean(n_per_taxon == 1), 1),
             taxa_Lauraceae_leaves = n_distinct(s$its_taxon[s$substrate == "Lauraceae leaves"]),
             taxa_Ficus_leaves     = n_distinct(s$its_taxon[s$substrate == "Ficus leaves"]),
             taxa_Ficus_wood       = n_distinct(s$its_taxon[s$substrate == "Ficus wood"]))
}))
write.csv(overview, file.path(cmp_dir, "tables", "unit_set_overview.csv"), row.names = FALSE)

# ---- Observed vs estimated total richness ----
compl <- combine("sampling_completeness.csv")
if (nrow(compl)) {
  d <- compl %>%
    select(unit_set, substrate, Observed = observed_S, Chao1 = chao1, ACE = ace, chao1_se) %>%
    pivot_longer(c(Observed, Chao1, ACE), names_to = "estimator", values_to = "richness") %>%
    mutate(estimator = factor(estimator, levels = c("Observed", "Chao1", "ACE")),
           se = ifelse(estimator == "Chao1", chao1_se, NA))
  p <- ggplot(d, aes(substrate, richness, fill = estimator)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = richness - se, ymax = richness + se),
                  position = position_dodge(width = 0.8), width = 0.15, na.rm = TRUE) +
    geom_text(aes(label = round(richness)), position = position_dodge(width = 0.8),
              vjust = -0.6, size = 3) +
    facet_wrap(~ unit_set, nrow = 1, labeller = as_labeller(set_names)) +
    scale_fill_manual(values = c(Observed = "#2E86AB", Chao1 = "#F18F01", ACE = "#A23B72")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(title = "Observed vs estimated total richness, per unit set",
         subtitle = "Error bar = 1 SE of Chao1",
         x = NULL, y = "Number of taxa", fill = NULL) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), legend.position = "top",
          axis.text.x = element_text(angle = 20, hjust = 1))
  ggsave(file.path(cmp_dir, "plots", "richness_observed_vs_estimated_all_unit_sets.pdf"),
         p, width = 15, height = 6.5)

  p <- ggplot(compl, aes(substrate, goods_coverage, fill = unit_set)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = paste0(goods_coverage, "%")), position = position_dodge(width = 0.8),
              vjust = -0.6, size = 3) +
    labs(title = "Good's coverage per unit set", x = NULL, y = "Good's coverage (%)",
         fill = "Unit set") +
    theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(cmp_dir, "plots", "goods_coverage_all_unit_sets.pdf"), p, width = 10, height = 6)
}

combine("richness_at_equal_coverage.csv")
combine("alpha_diversity_per_unit_summary.csv")
combine("alpha_diversity_tests.csv")
combine("alpha_diversity_pairwise.csv")
combine("lineage_coverage.csv")
combine("substrate_multilevel_summary.csv")
combine("substrate_pairwise_multilevel.csv")
combine("ficus_wood_zones_multilevel_summary.csv")
combine("orientation_multilevel_summary.csv")
combine("orientation_circular_gradient.csv")
combine("lifestyle_fisher_trees.csv")
combine("lifestyle_by_substrate.csv")

# ---- Per-unit alpha diversity, side by side ----
apu <- combine("alpha_diversity_per_unit.csv")
if (nrow(apu)) {
  for (idx in c("richness_S", "shannon_H")) {
    p <- ggplot(apu, aes(substrate, .data[[idx]], fill = substrate)) +
      geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.55) +
      geom_jitter(width = 0.12, height = 0, size = 1.8) +
      facet_wrap(~ unit_set, nrow = 1, labeller = as_labeller(set_names)) +
      scale_fill_manual(values = substrate_colours_3, guide = "none") +
      labs(title = paste("Per-sampling-unit", idx, "- all unit sets"), x = NULL, y = idx) +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold"), axis.text.x = element_text(angle = 20, hjust = 1))
    ggsave(file.path(cmp_dir, "plots", paste0("alpha_per_unit_", idx, "_all_unit_sets.pdf")),
           p, width = 14, height = 6)
  }
}

# ---- Multivariate headline tests, side by side ----
mv <- combine("multivariate_summary.csv")
if (nrow(mv)) {
  mv$analysis <- factor(mv$analysis, levels = rev(unique(mv$analysis)))
  p <- ggplot(mv, aes(permanova_R2, analysis, colour = unit_set, shape = permanova_p < 0.05)) +
    geom_point(size = 3.5, position = position_dodge(width = 0.6)) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "PERMANOVA p < 0.05") +
    labs(title = "PERMANOVA effect size (R\u00b2) per analysis and unit set",
         x = "R\u00b2", y = NULL, colour = "Unit set") +
    theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(cmp_dir, "plots", "permanova_R2_all_unit_sets.pdf"), p, width = 11, height = 7)
}

cat("Saved: results/comparison\n")
