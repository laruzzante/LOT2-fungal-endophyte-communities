# =============================================================
# Fungal lifestyle vs substrate / host tree
#
# The Lauraceae is being strangled by the Ficus, so a host effect is
# expected: does the weakened Lauraceae carry more pathogenic fungi than
# its strangler? Needs a 'Lifestyle' column (per isolate in
# LOT2_samples.xlsx, or per taxon in the ITS100 / ITS99 workbooks).
# Categories containing "pathogen" or "parasit" count as pathogenic.
# =============================================================
if (!have_lifestyle) {
  writeLines(c("Lifestyle analyses NOT RUN.",
               "No 'Lifestyle' column was found in LOT2_samples.xlsx or in the ITS100/ITS99 workbooks.",
               "Add one (one value per isolate, or per taxon) and re-run main.R."),
             "tables/lifestyle_NOT_RUN.txt")
  cat("Lifestyle analyses skipped (no Lifestyle column)\n")
} else {
  cat("\n====== Lifestyle analyses ======\n")
  dir.create("plots/lifestyle", recursive = TRUE, showWarnings = FALSE)

  ls_dat <- samples_raw %>%
    mutate(lifestyle  = ifelse(is.na(lifestyle), "Unknown", lifestyle),
           pathogenic = grepl("pathogen|parasit", lifestyle, ignore.case = TRUE),
           tree = ifelse(substrate == "Lauraceae leaves", "Lauraceae (strangled host)",
                         "Ficus (strangler)"))

  # 1. Composition per substrate and per tree
  compo <- function(d, grp) d %>%
    group_by(.data[[grp]], lifestyle) %>%
    summarise(n_isolates = n(), n_taxa = n_distinct(its_taxon), .groups = "drop_last") %>%
    mutate(pct_isolates = round(100 * n_isolates / sum(n_isolates), 1)) %>% ungroup()
  write.csv(compo(ls_dat, "substrate"), "tables/lifestyle_by_substrate.csv", row.names = FALSE)
  write.csv(compo(ls_dat, "tree"), "tables/lifestyle_by_tree.csv", row.names = FALSE)
  write.csv(compo(ls_dat %>% filter(zone == 6), "substrate"),
            "tables/lifestyle_by_substrate_zone6.csv", row.names = FALSE)

  p <- ggplot(compo(ls_dat, "substrate"), aes(substrate, n_isolates, fill = lifestyle)) +
    geom_col(position = "fill", width = 0.7, colour = "white") +
    scale_y_continuous(labels = scales::percent) +
    labs(title = "Fungal lifestyle composition per substrate",
         subtitle = paste("Share of isolates | unit set:", unit_label),
         x = NULL, y = "Share of isolates", fill = "Lifestyle") +
    theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"))
  ggsave("plots/lifestyle/lifestyle_composition_by_substrate.pdf", p, width = 10, height = 6)

  # 2. Is lifestyle associated with substrate? (known lifestyles only)
  known <- ls_dat %>% filter(lifestyle != "Unknown")
  tab <- table(known$substrate, known$lifestyle)
  set.seed(42)
  chi <- chisq.test(tab, simulate.p.value = TRUE, B = 9999)

  # 3. Each lifestyle in each pair of substrates / trees (Fisher, BH-adjusted)
  fisher_pairs <- function(d, grp) {
    out <- list()
    for (pr in combn(sort(unique(d[[grp]])), 2, simplify = FALSE))
      for (ct in c("Pathogenic (any)", sort(unique(d$lifestyle)))) {
        dd <- d[d[[grp]] %in% pr, ]
        hit <- if (ct == "Pathogenic (any)") dd$pathogenic else dd$lifestyle == ct
        m <- table(factor(dd[[grp]], levels = pr), factor(hit, levels = c(TRUE, FALSE)))
        ft <- fisher.test(m)
        out[[length(out) + 1]] <- data.frame(
          comparison = paste(pr, collapse = " vs "), lifestyle = ct,
          pct_group1 = round(100 * m[1, 1] / sum(m[1, ]), 1),
          pct_group2 = round(100 * m[2, 1] / sum(m[2, ]), 1),
          odds_ratio = round(unname(ft$estimate), 3), p_value = signif(ft$p.value, 4))
      }
    r <- bind_rows(out)
    r$p_adj_BH <- signif(p.adjust(r$p_value, "BH"), 4)
    r
  }
  write.csv(fisher_pairs(ls_dat, "substrate"), "tables/lifestyle_fisher_substrates.csv", row.names = FALSE)
  write.csv(fisher_pairs(ls_dat, "tree"), "tables/lifestyle_fisher_trees.csv", row.names = FALSE)

  # 4. Replicated version: share of pathogenic isolates per sampling unit
  unit_share <- ls_dat %>% group_by(sample_id, substrate, zone, tree) %>%
    summarise(n_isolates = n(), pct_pathogenic = 100 * mean(pathogenic),
              pct_unknown = 100 * mean(lifestyle == "Unknown"), .groups = "drop")
  write.csv(unit_share, "tables/lifestyle_pathogenic_share_per_unit.csv", row.names = FALSE)
  kw <- kruskal.test(pct_pathogenic ~ factor(substrate), data = unit_share)
  pw <- pairwise.wilcox.test(unit_share$pct_pathogenic, unit_share$substrate,
                             p.adjust.method = "holm", exact = FALSE)

  p <- ggplot(unit_share, aes(substrate, pct_pathogenic, fill = substrate)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.55) +
    geom_jitter(width = 0.12, height = 0, size = 2.2) +
    scale_fill_manual(values = substrate_colours_3, guide = "none") +
    labs(title = "Share of pathogenic isolates per sampling unit",
         subtitle = paste0("Kruskal-Wallis p = ", signif(kw$p.value, 3), " | unit set: ", unit_label),
         x = NULL, y = "% pathogenic isolates") +
    theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"))
  ggsave("plots/lifestyle/pathogenic_share_per_unit.pdf", p, width = 8, height = 6)

  # 5. Lifestyle composition per sampling unit (PERMANOVA, Bray-Curtis)
  ls_mat <- ls_dat %>% count(sample_id, lifestyle) %>%
    pivot_wider(names_from = lifestyle, values_from = n, values_fill = 0) %>% as.data.frame()
  rownames(ls_mat) <- ls_mat$sample_id
  ls_mat <- as.matrix(ls_mat[, -1, drop = FALSE])
  ls_meta <- unit_share[match(rownames(ls_mat), unit_share$sample_id), ]
  set.seed(42)
  perm_all <- adonis2(ls_mat ~ substrate, data = ls_meta, method = "bray", permutations = 999)
  z6 <- ls_meta$zone == 6
  set.seed(42)
  perm_z6 <- adonis2(ls_mat[z6, ] ~ substrate, data = ls_meta[z6, ], method = "bray", permutations = 999)

  sink("tables/lifestyle_tests.txt")
  cat("Lifestyle analyses | unit set:", unit_label, "\n")
  cat("Pathogenic = lifestyle containing 'pathogen' or 'parasit'\n\n")
  cat("== Substrate x lifestyle contingency (known lifestyles, isolates) ==\n"); print(tab)
  cat("\nChi-square, Monte-Carlo p (9999 sim.):\n"); print(chi)
  cat("\n== % pathogenic isolates per sampling unit: Kruskal-Wallis ==\n"); print(kw)
  cat("\nPairwise Wilcoxon (Holm):\n"); print(pw$p.value)
  cat("\n== Lifestyle composition per sampling unit ~ substrate (PERMANOVA) ==\n"); print(perm_all)
  cat("\n== Same, zone 6 (top canopy) only ==\n"); print(perm_z6)
  sink()
  cat("Saved: lifestyle tables and plots\n")
}
