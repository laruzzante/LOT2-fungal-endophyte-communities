library(vegan)
library(indicspecies)
library(ggVennDiagram)

dir.create("plots/community",   recursive = TRUE, showWarnings = FALSE)
dir.create("tables",            recursive = TRUE, showWarnings = FALSE)

tax_hierarchy <- c("phylum", "subphylum", "superclass", "class",
                   "subclass", "order", "family", "genus", "species")
uncertain_values <- c("?", "NO", "incertae sedis")
incertae_label   <- "Incertae sedis"

substrate_colours_3 <- c("Ficus leaves"     = "#A23B72",
                          "Ficus wood"       = "#F18F01",
                          "Lauraceae leaves" = "#2E86AB")

# Clean pooled data
dat_agg <- pooled %>%
  mutate(across(all_of(tax_hierarchy),
                ~ ifelse(is.na(.) | . %in% uncertain_values, incertae_label, .)),
         its_taxon = ifelse(is.na(its_taxon) | its_taxon %in% uncertain_values,
                            incertae_label, its_taxon),
         across(c(n_fungi_laur_leaf, n_fungi_fic_leaf, n_fungi_fic_wood),
                ~ replace_na(., 0)))

# ------------------------------------------------------------
# Per-isolate lineage: link each sample isolate to its full
# taxonomy so that community analyses can be repeated at every
# taxonomic rank (an isolate that is Incertae sedis at species
# level may still carry a defined genus / family / order ...).
# Incertae sedis / unresolved values are set to NA at the rank
# where they are unresolved, so they are simply dropped from
# every rank-specific analysis (they are kept only in the
# dedicated abundance bar charts and pie charts).
# ------------------------------------------------------------
clean_tax <- function(x) ifelse(is.na(x) | x %in% uncertain_values, NA_character_, x)

# One lineage row per ITS genotype (from the pooled taxonomy table)
lineage_its <- pooled %>%
  distinct(its_taxon, .keep_all = TRUE) %>%
  transmute(its_taxon, across(all_of(tax_hierarchy), clean_tax))

# Genus -> higher-rank lineage (majority vote); used as a fallback
# for the few sample ITS names that do not match the pooled list.
lineage_genus <- pooled %>%
  mutate(genus = clean_tax(genus)) %>%
  filter(!is.na(genus)) %>%
  group_by(genus) %>%
  summarise(across(c(phylum, subphylum, superclass, class, subclass, order, family),
                   ~ { v <- clean_tax(.x); v <- v[!is.na(v)]
                       if (length(v)) names(sort(table(v), decreasing = TRUE))[1] else NA_character_ }),
            .groups = "drop")

samples_lineage <- samples_raw %>% left_join(lineage_its, by = "its_taxon")

# Fill lineage from the leading genus token where the exact ITS
# name was not found in the pooled taxonomy (name variants).
need_fix <- is.na(samples_lineage$phylum) & is.na(samples_lineage$genus)
if (any(need_fix)) {
  token <- sub(" .*", "", samples_lineage$its_taxon[need_fix])
  fb <- lineage_genus[match(token, lineage_genus$genus), ]
  for (col in c("phylum","subphylum","superclass","class","subclass","order","family")) {
    samples_lineage[[col]][need_fix] <- fb[[col]]
  }
  samples_lineage$genus[need_fix] <- ifelse(token %in% lineage_genus$genus, token, NA_character_)
}

lineage_levels <- c("phylum", "class", "order", "family", "genus", "species", "its_taxon")

cat("Isolate lineage coverage per rank:\n")
lineage_coverage <- sapply(lineage_levels, function(lv) {
  v <- if (lv == "its_taxon") samples_lineage$its_taxon else clean_tax(samples_lineage[[lv]])
  sum(!is.na(v))
})
for (lv in lineage_levels)
  cat(sprintf("  %-9s %d / %d\n", lv, lineage_coverage[[lv]], nrow(samples_lineage)))
write.csv(data.frame(rank = lineage_levels,
                     isolates_resolved = as.integer(lineage_coverage),
                     isolates_total = nrow(samples_lineage)),
          "tables/lineage_coverage.csv", row.names = FALSE)

# ============================================================
# A. ALPHA DIVERSITY (from pooled data, per substrate)
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

  # Exclude the pooled Incertae sedis node at this rank
  lvl_col <- if (level == "its_taxon") "its_taxon" else level
  agg <- agg[agg[[lvl_col]] != incertae_label, , drop = FALSE]

  mat <- as.matrix(agg[, c("Lauraceae_leaves", "Ficus_leaves", "Ficus_wood")])
  comm <- t(mat)

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

for (idx_name in c("shannon_H", "simpson_1mD", "inv_simpson", "pielou_J")) {
  p <- ggplot(diversity_table,
              aes(x = taxonomic_level, y = .data[[idx_name]], fill = substrate)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = substrate_colours_3) +
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
# B1. RANK-ABUNDANCE CURVES
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

  # Exclude the pooled Incertae sedis node at this rank
  lvl_col <- if (level == "its_taxon") "its_taxon" else level
  agg <- agg[agg[[lvl_col]] != incertae_label, , drop = FALSE]

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
    scale_colour_manual(values = substrate_colours_3) +
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

  # Exclude the pooled Incertae sedis node at this rank
  agg <- agg[agg$label != incertae_label, , drop = FALSE]

  plot_data <- agg %>%
    select(label, Lauraceae_leaves, Ficus_leaves, Ficus_wood) %>%
    group_by(label) %>%
    summarise(across(everything(), sum), .groups = "drop") %>%
    pivot_longer(-label, names_to = "substrate", values_to = "count") %>%
    mutate(substrate = gsub("_", " ", substrate)) %>%
    group_by(substrate) %>%
    mutate(rel = count / sum(count)) %>%
    ungroup()

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
# B3. VENN DIAGRAMS
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

  # Exclude the pooled Incertae sedis node at this rank
  agg <- agg[agg$taxon != incertae_label, , drop = FALSE]

  taxa_lists <- list(
    `Lauraceae leaves` = agg$taxon[agg$Lauraceae_leaves > 0],
    `Ficus leaves`     = agg$taxon[agg$Ficus_leaves > 0],
    `Ficus wood`       = agg$taxon[agg$Ficus_wood > 0]
  )

  p <- ggVennDiagram(taxa_lists, label_alpha = 0) +
    scale_fill_gradient(low = "#F4FAFE", high = "#2E86AB") +
    labs(title = paste("Shared taxa —", level)) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.margin = margin(15, 25, 15, 25))

  ggsave(file.path("plots/community", paste0("venn_", level, ".pdf")),
         plot = p, width = 10, height = 7)
}
cat("Saved: Venn diagrams\n")

# ============================================================
# HELPER: Build community matrix from sample-level data
# ============================================================
build_comm_matrix <- function(data, group_col = "sample_id") {
  comm_wide <- data %>%
    count(.data[[group_col]], its_taxon) %>%
    pivot_wider(names_from = its_taxon, values_from = n, values_fill = 0) %>%
    as.data.frame()
  rownames(comm_wide) <- comm_wide[[group_col]]
  as.matrix(comm_wide[, -1])
}

# ============================================================
# HELPER: Build community matrix at an arbitrary taxonomic rank
# from the lineage-annotated sample data. Isolates that are
# unresolved (Incertae sedis) at the requested rank are dropped.
# ============================================================
build_comm_matrix_level <- function(data, level, group_col = "sample_id") {
  val <- if (level == "its_taxon") data$its_taxon else clean_tax(data[[level]])
  d <- data.frame(grp = data[[group_col]], taxon = val, stringsAsFactors = FALSE)
  d <- d[!is.na(d$taxon), , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  cm <- d %>% count(grp, taxon) %>%
    pivot_wider(names_from = taxon, values_from = n, values_fill = 0) %>%
    as.data.frame()
  rownames(cm) <- cm$grp
  as.matrix(cm[, -1, drop = FALSE])
}

# ============================================================
# HELPER: Repeat the multivariate community comparison across
# every taxonomic rank for one grouping factor. Returns a tidy
# summary (PERMANOVA / ANOSIM / PERMDISP + NMDS stress per rank)
# and optionally writes one NMDS ordination plot per rank.
# ============================================================
run_level_sweep <- function(sel, group_var, colours, plot_prefix = NULL,
                            plotdir = "plots/community",
                            meta = meta_full, id_col = "sample_id",
                            data = samples_lineage) {
  meta_s <- meta[sel, , drop = FALSE]
  data_s <- data[data[[id_col]] %in% meta_s[[id_col]], , drop = FALSE]
  out <- list()
  for (lv in lineage_levels) {
    cm <- build_comm_matrix_level(data_s, lv, id_col)
    if (is.null(cm)) next
    cm <- cm[, colSums(cm) > 0, drop = FALSE]
    keep <- rowSums(cm) > 0
    cm <- cm[keep, , drop = FALSE]
    m <- meta_s[match(rownames(cm), meta_s[[id_col]]), , drop = FALSE]
    g <- m[[group_var]]
    if (nrow(cm) < 3 || length(unique(g)) < 2 || ncol(cm) < 2) next

    set.seed(42)
    nm <- tryCatch(metaMDS(cm, distance = "bray", k = 2, trymax = 100, trace = 0),
                   error = function(e) NULL)
    stress <- if (!is.null(nm)) round(nm$stress, 4) else NA_real_
    set.seed(42)
    pm <- tryCatch(adonis2(cm ~ g, data = data.frame(g = g),
                           method = "bray", permutations = 999),
                   error = function(e) NULL)
    set.seed(42)
    an <- tryCatch(anosim(cm, g, distance = "bray", permutations = 999),
                   error = function(e) NULL)
    # Only interpretable when every group has >= 3 samples (see the
    # dispersion block in run_multivariate)
    bd <- if (min(table(g)) >= 3)
            tryCatch(permutest(betadisper(vegdist(cm, "bray"), g), permutations = 999),
                     error = function(e) NULL)
          else NULL

    out[[lv]] <- data.frame(
      level        = lv,
      n_samples    = nrow(cm),
      n_taxa       = ncol(cm),
      # Fraction of sample pairs with no taxon in common: above ~0.5 the
      # Bray-Curtis matrix is saturated and the NMDS degenerates.
      prop_no_shared = round(mean(no.shared(cm)), 4),
      nmds_stress  = stress,
      permanova_F  = if (!is.null(pm)) round(pm$F[1], 4)  else NA_real_,
      permanova_R2 = if (!is.null(pm)) round(pm$R2[1], 4) else NA_real_,
      permanova_p  = if (!is.null(pm)) pm$`Pr(>F)`[1]     else NA_real_,
      anosim_R     = if (!is.null(an)) round(an$statistic, 4) else NA_real_,
      anosim_p     = if (!is.null(an)) an$signif          else NA_real_,
      betadisper_p = if (!is.null(bd)) round(bd$tab$`Pr(>F)`[1], 4) else NA_real_,
      stringsAsFactors = FALSE
    )

    degenerate <- !is.na(stress) && stress < 0.01 && nrow(cm) > 4

    if (!is.null(plot_prefix) && !is.null(nm)) {
      sc <- as.data.frame(scores(nm, display = "sites")); sc$grp <- g
      cent <- sc %>% group_by(grp) %>%
        summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2), .groups = "drop")
      pl <- ggplot(sc, aes(NMDS1, NMDS2, colour = grp)) +
        geom_point(size = 3) +
        geom_point(data = cent, shape = 4, size = 5, stroke = 1.5) +
        scale_colour_manual(values = colours)

      # 95% confidence ellipses, matching the run_multivariate plots.
      # stat_ellipse needs at least 4 points per group: with fewer it
      # emits "Too few points to calculate an ellipse" and draws nothing,
      # so small groups are filtered out here rather than passed in. The
      # fill legend is suppressed because the colour legend already
      # names every group.
      ell_n <- table(sc$grp)
      ell_grps <- names(ell_n[ell_n >= 4])
      if (length(ell_grps) > 0) {
        pl <- pl +
          stat_ellipse(data = sc[sc$grp %in% ell_grps, , drop = FALSE],
                       aes(fill = grp), geom = "polygon",
                       alpha = 0.15, level = 0.95, linetype = 2) +
          scale_fill_manual(values = colours, guide = "none")
      }

      pl <- pl +
        labs(title = paste0("NMDS at ", lv, " level \u2014 ", plot_prefix),
             subtitle = paste0("Bray-Curtis | Stress = ", stress,
                               " | ", ncol(cm), " taxa | ",
                               sprintf("%.0f%%", 100 * mean(no.shared(cm))),
                               " of sample pairs share no taxa",
                               if (degenerate)
                                 "\nDEGENERATE ORDINATION - too sparse to interpret"
                               else ""),
             colour = group_var) +
        theme_minimal(base_size = 12) +
        theme(plot.title = element_text(face = "bold"), legend.position = "top")
      ggsave(file.path(plotdir, paste0(plot_prefix, "_", lv, "_nmds.pdf")),
             plot = pl, width = 9, height = 6.5)
    }
  }
  bind_rows(out)
}

# ============================================================
# HELPER: Run full multivariate analysis suite
# ============================================================
run_multivariate <- function(comm_mat, meta, group_var, label, outdir,
                             colours = NULL, do_indval = TRUE) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  prefix <- file.path(outdir, label)

  groups <- meta[[group_var]]
  n_groups <- length(unique(groups))
  n_samples <- nrow(comm_mat)

  cat("\n--- ", label, ": ", n_samples, " samples, ", ncol(comm_mat),
      " taxa, ", n_groups, " groups ---\n", sep = "")

  if (n_samples < 3 || n_groups < 2) {
    cat("  Skipping: insufficient samples or groups\n")
    return(NULL)
  }

  # Remove empty taxa columns
  comm_mat <- comm_mat[, colSums(comm_mat) > 0, drop = FALSE]

  # Remove samples with zero total abundance
  keep <- rowSums(comm_mat) > 0
  if (sum(keep) < 3) {
    cat("  Skipping: too few non-empty samples\n")
    return(NULL)
  }
  comm_mat <- comm_mat[keep, , drop = FALSE]
  meta <- meta[keep, , drop = FALSE]
  groups <- meta[[group_var]]

  results <- list(n_samples = nrow(comm_mat), n_taxa = ncol(comm_mat),
                  n_groups = length(unique(groups)))

  if (is.null(colours)) {
    n <- length(unique(groups))
    colours <- setNames(scales::hue_pal()(n), sort(unique(groups)))
  }

  # ---- NMDS ----
  set.seed(42)
  nmds <- tryCatch(
    metaMDS(comm_mat, distance = "bray", k = 2, trymax = 200, trace = 0),
    error = function(e) { cat("  NMDS failed:", e$message, "\n"); NULL }
  )

  # How sparse is this matrix? When most sample pairs share no taxa at
  # all, Bray-Curtis saturates at 1 and the ordination has almost no
  # gradient left to recover: metaMDS then returns a (near) zero stress
  # that means "degenerate solution", not "excellent fit".
  prop_no_shared <- mean(no.shared(comm_mat))
  results$prop_no_shared <- round(prop_no_shared, 4)
  results$prop_singleton_taxa <- round(mean(colSums(comm_mat) == 1), 4)

  if (!is.null(nmds)) {
    results$nmds_stress <- round(nmds$stress, 4)
    results$nmds_degenerate <- nmds$stress < 0.01 && nrow(comm_mat) > 4
    cat("  NMDS stress:", results$nmds_stress,
        if (isTRUE(results$nmds_degenerate)) " (DEGENERATE - see sparsity)" else "", "\n")
    cat("  Sample pairs sharing no taxa:",
        sprintf("%.1f%%", 100 * prop_no_shared),
        " singleton taxa:", sprintf("%.1f%%", 100 * results$prop_singleton_taxa), "\n")

    nmds_scores <- as.data.frame(scores(nmds, display = "sites"))
    nmds_scores[[group_var]] <- groups

    centroids <- nmds_scores %>%
      group_by(.data[[group_var]]) %>%
      summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2), .groups = "drop")

    p_nmds <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2,
                                       colour = .data[[group_var]])) +
      geom_point(size = 3) +
      geom_point(data = centroids, shape = 4, size = 5, stroke = 1.5) +
      scale_colour_manual(values = colours) +
      labs(title = paste("NMDS —", label),
           subtitle = paste0(
             "Bray-Curtis | Stress = ", signif(nmds$stress, 3),
             " | ", sprintf("%.0f%%", 100 * prop_no_shared),
             " of sample pairs share no taxa",
             if (isTRUE(results$nmds_degenerate))
               "\nDEGENERATE ORDINATION - too sparse to interpret; use the higher-rank plots"
             else ""),
           colour = group_var) +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold"),
            legend.position = "top",
            legend.margin = margin(0, 0, 5, 0),
            plot.margin = margin(10, 15, 10, 15)) +
      guides(colour = guide_legend(nrow = 2))

    grp_n <- table(groups)
    if (any(grp_n >= 3)) {
      ellipse_data <- nmds_scores[groups %in% names(grp_n[grp_n >= 3]), ]
      p_nmds <- p_nmds +
        stat_ellipse(data = ellipse_data,
                     aes(fill = .data[[group_var]]),
                     geom = "polygon", alpha = 0.15, level = 0.95, linetype = 2) +
        scale_fill_manual(values = colours) +
        guides(fill = guide_legend(nrow = 2))
    }

    ggsave(paste0(prefix, "_nmds.pdf"), plot = p_nmds, width = 10, height = 7.5)

    sp_scores <- as.data.frame(scores(nmds, display = "species"))
    sp_scores$taxon <- rownames(sp_scores)
    top_n <- min(15, ncol(comm_mat))
    top_taxa <- names(sort(colSums(comm_mat), decreasing = TRUE))[1:top_n]
    sp_top <- sp_scores[sp_scores$taxon %in% top_taxa, ]

    p_nmds_sp <- p_nmds +
      geom_text(data = sp_top, aes(x = NMDS1, y = NMDS2, label = taxon),
                inherit.aes = FALSE, size = 2.5, alpha = 0.7, fontface = "italic") +
      labs(title = paste("NMDS with top", top_n, "taxa —", label))
    ggsave(paste0(prefix, "_nmds_species.pdf"), plot = p_nmds_sp, width = 12, height = 8.5)
  }

  # ---- PERMANOVA ----
  set.seed(42)
  form <- as.formula(paste("comm_mat ~", group_var))
  perm <- tryCatch(
    adonis2(form, data = meta, method = "bray", permutations = 999),
    error = function(e) { cat("  PERMANOVA failed:", e$message, "\n"); NULL }
  )

  if (!is.null(perm)) {
    results$permanova <- perm
    cat("  PERMANOVA F:", round(perm$F[1], 4), " R2:", round(perm$R2[1], 4),
        " p:", perm$`Pr(>F)`[1], "\n")

    sink(paste0(prefix, "_permanova.txt"))
    cat("PERMANOVA — Bray-Curtis distance\n")
    cat("Formula: community ~", group_var, "\n")
    cat("Permutations: 999\n\n")
    print(perm)
    sink()
  }

  # ---- ANOSIM ----
  set.seed(42)
  anos <- tryCatch(
    anosim(comm_mat, grouping = groups, distance = "bray", permutations = 999),
    error = function(e) { cat("  ANOSIM failed:", e$message, "\n"); NULL }
  )

  if (!is.null(anos)) {
    results$anosim_R <- round(anos$statistic, 4)
    results$anosim_p <- anos$signif
    cat("  ANOSIM R:", results$anosim_R, " p:", results$anosim_p, "\n")

    sink(paste0(prefix, "_anosim.txt"))
    cat("ANOSIM — Bray-Curtis distance\n")
    cat("Permutations: 999\n\n")
    cat("R statistic:", round(anos$statistic, 4), "\n")
    cat("p-value:", anos$signif, "\n\n")
    print(summary(anos))
    sink()
  }

  # ---- Beta-dispersion ----
  bc_dist <- vegdist(comm_mat, method = "bray")
  bd <- tryCatch(
    betadisper(bc_dist, groups),
    error = function(e) { cat("  betadisper failed:", e$message, "\n"); NULL }
  )

  if (!is.null(bd)) {
    bd_perm <- permutest(bd, permutations = 999)
    # A group of 1-2 samples has a degenerate distance-to-centroid
    # (zero, or the same value for both points), which can drive the
    # permutest F to absurd magnitudes. The test is only interpretable
    # when every group holds at least 3 samples.
    min_grp_n <- min(table(groups))
    results$betadisper_min_group_n <- as.integer(min_grp_n)
    results$betadisper_F <- round(bd_perm$tab$F[1], 4)
    results$betadisper_p <- round(bd_perm$tab$`Pr(>F)`[1], 4)

    if (min_grp_n < 3) {
      cat("  Betadisper not interpretable (smallest group has", min_grp_n,
          "sample(s)) - reported as NA\n")
      results$betadisper_note <- paste0(
        "not interpretable: smallest group has ", min_grp_n, " sample(s)")
      results$betadisper_F <- NA_real_
      results$betadisper_p <- NA_real_
    } else {
      cat("  Betadisper F:", results$betadisper_F, " p:", results$betadisper_p, "\n")
    }

    sink(paste0(prefix, "_betadisper.txt"))
    cat("Beta-dispersion test (betadisper + permutest)\n\n")
    if (min_grp_n < 3) {
      cat("WARNING: the smallest group holds only ", min_grp_n,
          " sample(s). Distances to a centroid computed from 1-2 points are\n",
          "degenerate, so the F ratio and p-value below are NOT interpretable\n",
          "and are reported as NA in the summary tables. They are printed here\n",
          "for transparency only.\n\n", sep = "")
    }
    print(bd_perm)
    cat("\nGroup mean distances to centroid:\n")
    print(data.frame(group = levels(bd$group),
                     n = as.integer(table(bd$group)),
                     mean_dist = round(tapply(bd$distances, bd$group, mean), 4)))
    sink()
  }

  # ---- Pairwise PERMANOVA ----
  grp_levels <- unique(groups)
  if (length(grp_levels) >= 2) {
    pairs <- combn(grp_levels, 2, simplify = FALSE)
    pw_results <- list()
    for (pr in pairs) {
      sel <- groups %in% pr
      if (sum(sel) >= 3) {
        set.seed(42)
        pw <- tryCatch(
          adonis2(comm_mat[sel, ] ~ grp, data = data.frame(grp = groups[sel]),
                  method = "bray", permutations = 999),
          error = function(e) NULL
        )
        if (!is.null(pw)) {
          pw_results[[paste(pr, collapse = " vs ")]] <- data.frame(
            pair    = paste(pr, collapse = " vs "),
            F_stat  = round(pw$F[1], 4),
            R2      = round(pw$R2[1], 4),
            p_value = pw$`Pr(>F)`[1]
          )
        }
      }
    }
    if (length(pw_results) > 0) {
      pw_table <- bind_rows(pw_results)
      pw_table$p_adj <- p.adjust(pw_table$p_value, method = "holm")
      results$pairwise <- pw_table
      write.csv(pw_table, paste0(prefix, "_pairwise_permanova.csv"), row.names = FALSE)
      cat("  Pairwise PERMANOVA saved\n")
    }
  }

  # ---- Rarefaction ----
  grp_levels <- unique(groups)
  comm_by_grp <- list()
  for (g in grp_levels) {
    rows <- groups == g
    comm_by_grp[[g]] <- colSums(comm_mat[rows, , drop = FALSE])
  }
  comm_grp_mat <- do.call(rbind, comm_by_grp)

  rarefy_data <- list()
  for (g in grp_levels) {
    n_total <- sum(comm_grp_mat[g, ])
    if (n_total < 2) next
    steps <- unique(c(seq(1, n_total, length.out = 50), n_total))
    rare <- rarefy(comm_grp_mat[g, , drop = FALSE], sample = floor(steps))
    rarefy_data[[g]] <- data.frame(
      group = g,
      n_individuals = floor(steps),
      expected_species = as.numeric(rare)
    )
  }
  if (length(rarefy_data) > 0) {
    rare_df <- bind_rows(rarefy_data)
    p_rare <- ggplot(rare_df, aes(x = n_individuals, y = expected_species,
                                   colour = group)) +
      geom_line(linewidth = 1) +
      scale_colour_manual(values = colours) +
      labs(title = paste("Rarefaction curves —", label),
           x = "Number of individuals", y = "Expected number of taxa (ITS)",
           colour = group_var) +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold"), legend.position = "top")
    ggsave(paste0(prefix, "_rarefaction.pdf"), plot = p_rare, width = 9, height = 6)
    cat("  Rarefaction curves saved\n")
  }

  # ---- Indicator species ----
  if (do_indval && n_groups >= 2 && n_samples >= 5) {
    set.seed(42)
    indval <- tryCatch(
      multipatt(comm_mat, cluster = groups, func = "IndVal.g",
                control = how(nperm = 999)),
      error = function(e) { cat("  IndVal failed:", e$message, "\n"); NULL }
    )

    if (!is.null(indval)) {
      indval_summary <- capture.output(summary(indval, indvalcomp = TRUE))
      writeLines(indval_summary, paste0(prefix, "_indicator_species.txt"))

      sig <- indval$sign[indval$sign$p.value <= 0.05, , drop = FALSE]
      if (nrow(sig) > 0) {
        sig$taxon <- rownames(sig)
        sig <- sig %>% arrange(p.value)
        write.csv(sig, paste0(prefix, "_indicator_species_significant.csv"),
                  row.names = FALSE)
        results$n_indicators <- nrow(sig)
        cat("  Significant indicators:", nrow(sig), "\n")
      } else {
        results$n_indicators <- 0
        cat("  No significant indicators\n")
      }
    }
  }

  results
}

# ============================================================
# C. BUILD COMMUNITY MATRIX FROM SAMPLES DATA
# ============================================================
cat("\n====== C. Building sample-level community matrix ======\n")

comm_mat_full <- build_comm_matrix(samples_raw, "sample_id")

meta_full <- samples_raw %>%
  distinct(sample_id, substrate, zone, unit, position) %>%
  as.data.frame()
rownames(meta_full) <- meta_full$sample_id
meta_full <- meta_full[rownames(comm_mat_full), ]

cat("Full community matrix:", nrow(comm_mat_full), "samples x",
    ncol(comm_mat_full), "taxa\n")

# ============================================================
# C1. SUBSTRATE COMPARISON: All 3 substrates
# ============================================================
cat("\n====== C1. All substrates comparison ======\n")

results_substrate <- run_multivariate(
  comm_mat_full, meta_full, "substrate",
  label = "substrate_all",
  outdir = "tables",
  colours = substrate_colours_3
)

# Copy key results to standard table locations
if (!is.null(results_substrate$permanova)) {
  sink("tables/permanova_results.txt")
  cat("PERMANOVA — Bray-Curtis distance\n")
  cat("Formula: community ~ substrate\nPermutations: 999\n\n")
  print(results_substrate$permanova)
  sink()
}

# ============================================================
# C2. SUBSTRATE SUBANALYSIS: Ficus leaves vs Lauraceae leaves
# ============================================================
cat("\n====== C2. Ficus leaves vs Lauraceae leaves ======\n")

sel_leaves <- meta_full$substrate %in% c("Ficus leaves", "Lauraceae leaves")
results_leaves <- run_multivariate(
  comm_mat_full[sel_leaves, ], meta_full[sel_leaves, ], "substrate",
  label = "substrate_leaves",
  outdir = "tables",
  colours = substrate_colours_3[c("Ficus leaves", "Lauraceae leaves")]
)

# ============================================================
# E. ZONE-BASED ANALYSES
# ============================================================
cat("\n====== E. Zone-based analyses ======\n")

zone_colours <- c("1" = "#440154", "2" = "#3B528B", "3" = "#21918C",
                   "4" = "#5EC962", "5" = "#ADDC30", "6" = "#FDE725")
position_colours <- c("Trunk" = "#8B4513", "Branch" = "#228B22")

# ---- E1. Ficus wood: all zones ----
cat("\n--- E1. Ficus wood across zones ---\n")
sel_fw <- meta_full$substrate == "Ficus wood"
meta_fw <- meta_full[sel_fw, ]
meta_fw$zone_f <- as.character(meta_fw$zone)
comm_fw <- comm_mat_full[sel_fw, ]

results_fw_zones <- run_multivariate(
  comm_fw, meta_fw, "zone_f",
  label = "ficus_wood_zones",
  outdir = "tables",
  colours = zone_colours
)

# ---- E2. Ficus wood: trunk (zones 1-5) vs branch (zone 6) ----
cat("\n--- E2. Ficus wood trunk vs branch ---\n")

results_fw_trunkbranch <- run_multivariate(
  comm_fw, meta_fw, "position",
  label = "ficus_wood_trunk_vs_branch",
  outdir = "tables",
  colours = position_colours
)

# ---- E3. Trunk vs branch per substrate (where data permits) ----
cat("\n--- E3. Trunk vs branch per substrate ---\n")

subst_with_both <- meta_full %>%
  group_by(substrate) %>%
  summarise(has_trunk = any(position == "Trunk"),
            has_branch = any(position == "Branch"), .groups = "drop") %>%
  filter(has_trunk & has_branch) %>%
  pull(substrate)

results_trunk_branch_per_sub <- list()
if (length(subst_with_both) > 0) {
  for (sub in subst_with_both) {
    sel_sub <- meta_full$substrate == sub
    cat("\n  Trunk vs branch for:", sub, "\n")
    results_trunk_branch_per_sub[[sub]] <- run_multivariate(
      comm_mat_full[sel_sub, ], meta_full[sel_sub, ], "position",
      label = paste0(gsub(" ", "_", tolower(sub)), "_trunk_vs_branch"),
      outdir = "tables",
      colours = position_colours
    )
  }
}

# ---- E4. All substrates at zone 6 (branch level only) ----
cat("\n--- E4. Substrate comparison at zone 6 (branch level) ---\n")
sel_z6 <- meta_full$zone == 6
results_z6 <- run_multivariate(
  comm_mat_full[sel_z6, ], meta_full[sel_z6, ], "substrate",
  label = "substrate_zone6_branch",
  outdir = "tables",
  colours = substrate_colours_3
)

# ---- E5. Substrate x position interaction (combined factor) ----
cat("\n--- E5. Substrate x position interaction ---\n")
meta_full$sub_pos <- paste(meta_full$substrate, meta_full$position, sep = " - ")
n_combos <- length(unique(meta_full$sub_pos))
combo_colours <- setNames(
  scales::hue_pal()(n_combos),
  sort(unique(meta_full$sub_pos))
)

results_subpos <- run_multivariate(
  comm_mat_full, meta_full, "sub_pos",
  label = "substrate_x_position",
  outdir = "tables",
  colours = combo_colours
)

# ============================================================
# D. MULTI-LEVEL TAXONOMIC COMMUNITY ANALYSES
# Repeat the key multivariate comparisons at every taxonomic
# rank, using each isolate's full lineage. This shows whether
# the substrate / zone signal detected at ITS-genotype level
# persists when isolates are grouped into higher taxa.
# ============================================================
cat("\n====== D. Multi-level taxonomic community analyses ======\n")

meta_full$zone_f <- as.character(meta_full$zone)

cat("\n--- D1. Substrate comparison across taxonomic ranks ---\n")
sweep_substrate <- run_level_sweep(
  rep(TRUE, nrow(meta_full)), "substrate",
  substrate_colours_3, plot_prefix = "substrate_bylevel")
write.csv(sweep_substrate, "tables/substrate_multilevel_summary.csv", row.names = FALSE)
print(sweep_substrate)

cat("\n--- D2. Ficus wood zones across taxonomic ranks ---\n")
sweep_zones <- run_level_sweep(
  meta_full$substrate == "Ficus wood", "zone_f",
  zone_colours, plot_prefix = "ficus_wood_zones_bylevel")
write.csv(sweep_zones, "tables/ficus_wood_zones_multilevel_summary.csv", row.names = FALSE)
print(sweep_zones)

cat("\n--- D3. Substrate x position across taxonomic ranks ---\n")
sweep_subpos <- run_level_sweep(
  rep(TRUE, nrow(meta_full)), "sub_pos",
  combo_colours, plot_prefix = NULL)
write.csv(sweep_subpos, "tables/substrate_x_position_multilevel_summary.csv", row.names = FALSE)
print(sweep_subpos)

# ============================================================
# F. ZONE GRADIENT — FICUS WOOD
# ============================================================
cat("\n====== F. Ficus wood zone gradient test ======\n")

if (!is.null(results_fw_zones)) {
  meta_fw_num <- meta_fw
  meta_fw_num$zone_num <- meta_fw_num$zone
  set.seed(42)
  perm_zone_ord <- tryCatch(
    adonis2(comm_fw ~ zone_num, data = meta_fw_num, method = "bray",
            permutations = 999),
    error = function(e) NULL
  )
  if (!is.null(perm_zone_ord)) {
    sink("tables/ficus_wood_zones_permanova_continuous.txt")
    cat("PERMANOVA — Ficus wood community vs zone height (continuous)\n\n")
    print(perm_zone_ord)
    sink()
    cat("  Zone gradient PERMANOVA saved\n")
  }
}

# ============================================================
# G. BRANCH / TRUNK SAMPLING ORIENTATION
#
# The re-issued samples workbook reconciled the field branch
# numbering with the project database, which finally makes the
# compass orientation of each sampled branch / trunk face usable.
#
# Orientation is recorded per sampled face, and a single sampling
# unit can span more than one face, so these analyses use the
# finer `sample_id_orient` grouping (substrate x zone x unit x
# bearing). Isolates whose orientation could not be reconciled
# ("?" in the field sheet, plus all Ficus trunk wood) are dropped.
# ============================================================
cat("\n====== G. Branch/trunk sampling orientation ======\n")

orientation_colours <- c("N"  = "#1F4E79", "NE" = "#2E86AB", "E"  = "#4DAF4A",
                         "SE" = "#A6D854", "S"  = "#F18F01", "SW" = "#E7298A",
                         "W"  = "#A23B72", "NW" = "#7570B3")
aspect_colours <- c("Sun-facing (N)" = "#E8A33D",
                    "Shaded (S)"     = "#3B528B",
                    "Lateral (E/W)"  = "#21918C")

samples_orient <- samples_raw %>% filter(!is.na(orientation))

cat("Isolates with a resolved orientation:", nrow(samples_orient),
    "/", nrow(samples_raw), "\n")

if (nrow(samples_orient) == 0) {
  cat("No orientation data available - skipping Section G\n")
} else {

# ---- G0. Sampling coverage by orientation ----
orient_coverage <- samples_orient %>%
  group_by(substrate, position, orientation) %>%
  summarise(n_units    = n_distinct(sample_id_orient),
            n_isolates = n(),
            n_taxa     = n_distinct(its_taxon),
            .groups = "drop") %>%
  arrange(substrate, position, orientation)
write.csv(orient_coverage, "tables/orientation_coverage.csv", row.names = FALSE)
cat("\nSampling coverage per orientation:\n")
print(as.data.frame(orient_coverage))

unresolved <- samples_raw %>%
  filter(is.na(orientation)) %>%
  count(substrate, position, name = "n_isolates")
write.csv(unresolved, "tables/orientation_unresolved.csv", row.names = FALSE)
cat("\nIsolates without a usable orientation:\n")
print(as.data.frame(unresolved))

# ---- Orientation-level community matrix ----
comm_orient <- build_comm_matrix(samples_orient, "sample_id_orient")

meta_orient <- samples_orient %>%
  distinct(sample_id_orient, substrate, zone, unit, position,
           orientation, orientation_deg, sun_aspect) %>%
  as.data.frame()
rownames(meta_orient) <- meta_orient$sample_id_orient
meta_orient <- meta_orient[rownames(comm_orient), ]
meta_orient$zone_f  <- as.character(meta_orient$zone)
meta_orient$sub_ori <- paste(meta_orient$substrate, meta_orient$orientation, sep = " - ")

cat("\nOrientation community matrix:", nrow(comm_orient), "sampling units x",
    ncol(comm_orient), "ITS genotypes\n")

# ---- G1. Orientation across all substrates ----
cat("\n--- G1. Orientation, all substrates pooled ---\n")
results_orient_all <- run_multivariate(
  comm_orient, meta_orient, "orientation",
  label = "orientation_all",
  outdir = "tables",
  colours = orientation_colours)

if (!is.null(results_orient_all$permanova)) {
  sink("tables/orientation_permanova_results.txt")
  cat("PERMANOVA - Bray-Curtis distance\n")
  cat("Formula: community ~ orientation (all substrates pooled)\n")
  cat("Permutations: 999\n\n")
  print(results_orient_all$permanova)
  sink()
}

# ---- G2. Orientation within each substrate ----
cat("\n--- G2. Orientation within each substrate ---\n")

orient_by_substrate <- meta_orient %>%
  group_by(substrate) %>%
  summarise(n_units = n(), n_orient = n_distinct(orientation), .groups = "drop")
print(as.data.frame(orient_by_substrate))

results_orient_sub <- list()
for (sub in orient_by_substrate$substrate[orient_by_substrate$n_orient >= 2 &
                                          orient_by_substrate$n_units >= 3]) {
  cat("\n  Orientation within:", sub, "\n")
  sel_sub <- meta_orient$substrate == sub
  results_orient_sub[[sub]] <- run_multivariate(
    comm_orient[sel_sub, , drop = FALSE], meta_orient[sel_sub, , drop = FALSE],
    "orientation",
    label = paste0(gsub(" ", "_", tolower(sub)), "_orientation"),
    outdir = "tables",
    colours = orientation_colours)
}

# ---- G3. Substrate x orientation combined factor ----
cat("\n--- G3. Substrate x orientation ---\n")
n_so <- length(unique(meta_orient$sub_ori))
sub_ori_colours <- setNames(scales::hue_pal()(n_so), sort(unique(meta_orient$sub_ori)))
results_sub_ori <- run_multivariate(
  comm_orient, meta_orient, "sub_ori",
  label = "substrate_x_orientation",
  outdir = "tables",
  colours = sub_ori_colours)

# ---- G4. Sun aspect (southern-hemisphere exposure proxy) ----
# LOT2 was collected in Peru: the northern face of a tree is the
# sun-exposed one, the southern face the shaded one, E/W lateral.
cat("\n--- G4. Sun aspect (sun-facing / lateral / shaded) ---\n")
results_aspect <- run_multivariate(
  comm_orient, meta_orient, "sun_aspect",
  label = "orientation_sun_aspect",
  outdir = "tables",
  colours = aspect_colours)

results_aspect_sub <- list()
for (sub in unique(meta_orient$substrate)) {
  sel_sub <- meta_orient$substrate == sub
  if (length(unique(meta_orient$sun_aspect[sel_sub])) < 2) next
  cat("\n  Sun aspect within:", sub, "\n")
  results_aspect_sub[[sub]] <- run_multivariate(
    comm_orient[sel_sub, , drop = FALSE], meta_orient[sel_sub, , drop = FALSE],
    "sun_aspect",
    label = paste0(gsub(" ", "_", tolower(sub)), "_sun_aspect"),
    outdir = "tables",
    colours = aspect_colours,
    do_indval = FALSE)
}

# ---- G5. Variance partitioning: substrate vs orientation ----
# Marginal (type-III) PERMANOVA: how much community variation does
# orientation explain once substrate is accounted for, and vice versa?
cat("\n--- G5. Substrate + orientation variance partitioning ---\n")
set.seed(42)
perm_partition <- tryCatch(
  adonis2(comm_orient ~ substrate + orientation, data = meta_orient,
          method = "bray", permutations = 999, by = "margin"),
  error = function(e) { cat("  failed:", e$message, "\n"); NULL })

set.seed(42)
perm_sequential <- tryCatch(
  adonis2(comm_orient ~ substrate + orientation, data = meta_orient,
          method = "bray", permutations = 999, by = "terms"),
  error = function(e) NULL)

# The same question for the coarser sun-aspect grouping. This one
# matters: Lauraceae was only ever sampled on its northern/north-western
# faces, so aspect is partly confounded with substrate and the pooled
# aspect test on its own cannot separate the two.
set.seed(42)
perm_partition_aspect <- tryCatch(
  adonis2(comm_orient ~ substrate + sun_aspect, data = meta_orient,
          method = "bray", permutations = 999, by = "margin"),
  error = function(e) NULL)

aspect_confound <- table(meta_orient$substrate, meta_orient$sun_aspect)

if (!is.null(perm_partition)) {
  sink("tables/orientation_variance_partitioning.txt")
  cat("PERMANOVA - variance partitioning, Bray-Curtis distance\n")
  cat("Only sampling units with a reconciled orientation are used.\n")
  cat("Permutations: 999\n\n")
  cat("== Marginal (type III): each term adjusted for the other ==\n")
  print(perm_partition)
  if (!is.null(perm_sequential)) {
    cat("\n== Sequential (type I): substrate entered first ==\n")
    print(perm_sequential)
  }
  if (!is.null(perm_partition_aspect)) {
    cat("\n== Marginal (type III): substrate + sun aspect ==\n")
    print(perm_partition_aspect)
  }
  cat("\n== Sampling units per substrate x sun aspect ==\n")
  cat("(empty cells show where aspect is confounded with substrate)\n")
  print(aspect_confound)
  sink()
  cat("  Saved: tables/orientation_variance_partitioning.txt\n")
  print(perm_partition)
  if (!is.null(perm_partition_aspect)) {
    cat("\nSun aspect adjusted for substrate:\n")
    print(perm_partition_aspect)
    cat("\nSampling units per substrate x sun aspect:\n")
    print(aspect_confound)
  }
}

# ---- G6. Directional gradient test (circular) ----
# Treating the compass bearing as a circular variable rather than as
# discrete classes: a community that changes smoothly around the tree
# is captured by the sine/cosine pair with only 2 df, which is more
# powerful than a 4-level factor at this level of replication.
cat("\n--- G6. Circular (directional gradient) test ---\n")

circular_test <- function(cm, md, label) {
  if (nrow(cm) < 4 || length(unique(md$orientation)) < 3) return(NULL)
  rad <- md$orientation_deg * pi / 180
  df <- data.frame(cos_b = cos(rad), sin_b = sin(rad))
  set.seed(42)
  fit <- tryCatch(adonis2(cm ~ cos_b + sin_b, data = df, method = "bray",
                          permutations = 999, by = "margin"),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  # adonis2() with the default `by` collapses both terms into a single
  # "Model" row - that row is the joint 2-df directional test.
  set.seed(42)
  fit_joint <- tryCatch(adonis2(cm ~ cos_b + sin_b, data = df, method = "bray",
                                permutations = 999),
                        error = function(e) NULL)
  joint_row <- if (!is.null(fit_joint)) match("Model", rownames(fit_joint)) else NA_integer_
  data.frame(
    analysis      = label,
    n_units       = nrow(cm),
    n_orientations = length(unique(md$orientation)),
    # cos(bearing) contrasts the North-South axis, sin(bearing) East-West
    NS_axis_R2    = round(fit$R2[1], 4),
    NS_axis_p     = fit$`Pr(>F)`[1],
    EW_axis_R2    = round(fit$R2[2], 4),
    EW_axis_p     = fit$`Pr(>F)`[2],
    joint_R2      = if (!is.na(joint_row)) round(fit_joint$R2[joint_row], 4) else NA_real_,
    joint_p       = if (!is.na(joint_row)) fit_joint$`Pr(>F)`[joint_row] else NA_real_,
    stringsAsFactors = FALSE)
}

circ_rows <- list(circular_test(comm_orient, meta_orient, "All substrates"))
for (sub in unique(meta_orient$substrate)) {
  sel_sub <- meta_orient$substrate == sub
  circ_rows[[length(circ_rows) + 1]] <- circular_test(
    comm_orient[sel_sub, , drop = FALSE], meta_orient[sel_sub, , drop = FALSE], sub)
}
circ_table <- bind_rows(circ_rows)
if (nrow(circ_table) > 0) {
  write.csv(circ_table, "tables/orientation_circular_gradient.csv", row.names = FALSE)
  cat("Saved: tables/orientation_circular_gradient.csv\n")
  print(circ_table)
}

# ---- G7. Alpha diversity per orientation ----
cat("\n--- G7. Alpha diversity per orientation ---\n")

alpha_by <- function(data, grp_col, extra = NULL) {
  keys <- c(extra, grp_col)
  cm <- data %>% count(across(all_of(keys)), its_taxon) %>%
    pivot_wider(names_from = its_taxon, values_from = n, values_fill = 0) %>%
    as.data.frame()
  mat <- as.matrix(cm[, -seq_along(keys), drop = FALSE])
  S <- specnumber(mat); N <- rowSums(mat)
  H <- diversity(mat, "shannon")
  out <- cm[, keys, drop = FALSE]
  out$richness_S  <- S
  out$abundance_N <- N
  out$shannon_H   <- round(H, 4)
  out$simpson_1mD <- round(diversity(mat, "simpson"), 4)
  out$inv_simpson <- round(diversity(mat, "invsimpson"), 4)
  out$pielou_J    <- round(ifelse(S > 1, H / log(S), NA_real_), 4)
  # Richness rarefied to a common isolate count, so that unequal
  # sampling per bearing cannot drive the comparison. Groups holding
  # fewer than 10 isolates are too small to rarefy meaningfully: they
  # are excluded from the target and left as NA.
  min_n <- suppressWarnings(min(N[N >= 10]))
  if (is.finite(min_n)) {
    ok <- N >= min_n
    out$rarefied_S <- NA_real_
    out$rarefied_S[ok] <- round(
      as.numeric(rarefy(mat[ok, , drop = FALSE], sample = min_n)), 2)
    out$rarefied_to <- min_n
  } else {
    out$rarefied_S  <- NA_real_
    out$rarefied_to <- NA_integer_
  }
  out
}

alpha_orient <- alpha_by(samples_orient, "orientation")
alpha_orient_sub <- alpha_by(samples_orient, "orientation", extra = "substrate")
alpha_aspect <- alpha_by(samples_orient, "sun_aspect")

write.csv(alpha_orient, "tables/orientation_alpha_diversity.csv", row.names = FALSE)
write.csv(alpha_orient_sub, "tables/orientation_alpha_diversity_by_substrate.csv",
          row.names = FALSE)
write.csv(alpha_aspect, "tables/orientation_sun_aspect_alpha_diversity.csv",
          row.names = FALSE)
cat("\nAlpha diversity per orientation (all substrates):\n")
print(alpha_orient)
cat("\nAlpha diversity per sun aspect:\n")
print(alpha_aspect)

for (idx_name in c("richness_S", "rarefied_S", "shannon_H", "pielou_J")) {
  p <- ggplot(alpha_orient_sub,
              aes(x = orientation, y = .data[[idx_name]], fill = substrate)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = substrate_colours_3) +
    scale_x_discrete(limits = intersect(names(orientation_colours),
                                        alpha_orient_sub$orientation)) +
    labs(title = paste("Orientation vs", idx_name),
         subtitle = if (idx_name == "rarefied_S")
                      paste("Richness rarefied to", alpha_orient_sub$rarefied_to[1],
                            "isolates per group") else NULL,
         x = "Sampling orientation (compass bearing)",
         y = idx_name, fill = "Substrate") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), legend.position = "top")
  ggsave(file.path("plots/community", paste0("orientation_alpha_", idx_name, ".pdf")),
         plot = p, width = 9, height = 6)
}
cat("Saved: orientation alpha diversity plots\n")

# ---- G8. Taxa shared between orientations ----
cat("\n--- G8. Taxa shared between orientations ---\n")

orient_taxa <- split(samples_orient$its_taxon, samples_orient$orientation)
orient_taxa <- lapply(orient_taxa, unique)
if (length(orient_taxa) %in% 2:4) {
  p_venn <- ggVennDiagram(orient_taxa, label_alpha = 0) +
    scale_fill_gradient(low = "#F4FAFE", high = "#2E86AB") +
    labs(title = "ITS genotypes shared between sampling orientations") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.margin = margin(15, 25, 15, 25))
  ggsave("plots/community/orientation_venn_its_taxon.pdf",
         plot = p_venn, width = 10, height = 7)
  cat("Saved: orientation Venn diagram\n")
}

# Ficus leaves is the only substrate sampled on all four cardinal
# bearings, so it gives the cleanest four-way overlap picture.
fl_orient <- samples_orient %>% filter(substrate == "Ficus leaves")
fl_taxa <- lapply(split(fl_orient$its_taxon, fl_orient$orientation), unique)
if (length(fl_taxa) %in% 2:4) {
  p_venn_fl <- ggVennDiagram(fl_taxa, label_alpha = 0) +
    scale_fill_gradient(low = "#FDF2F7", high = "#A23B72") +
    labs(title = "Ficus leaves - ITS genotypes shared between orientations") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.margin = margin(15, 25, 15, 25))
  ggsave("plots/community/orientation_venn_ficus_leaves.pdf",
         plot = p_venn_fl, width = 10, height = 7)
}

orient_overlap <- samples_orient %>%
  distinct(orientation, its_taxon) %>%
  count(its_taxon, name = "n_orientations") %>%
  count(n_orientations, name = "n_taxa") %>%
  arrange(n_orientations)
write.csv(orient_overlap, "tables/orientation_taxon_sharing.csv", row.names = FALSE)
cat("\nHow many orientations each ITS genotype occurs in:\n")
print(as.data.frame(orient_overlap))

# ---- G9. Orientation across taxonomic ranks ----
cat("\n--- G9. Orientation across taxonomic ranks ---\n")

samples_lineage_orient <- samples_lineage %>% filter(!is.na(orientation))

sweep_orient <- run_level_sweep(
  rep(TRUE, nrow(meta_orient)), "orientation", orientation_colours,
  plot_prefix = "orientation_bylevel",
  meta = meta_orient, id_col = "sample_id_orient",
  data = samples_lineage_orient)
if (nrow(sweep_orient) > 0) {
  write.csv(sweep_orient, "tables/orientation_multilevel_summary.csv", row.names = FALSE)
  print(sweep_orient)
}

sweep_aspect <- run_level_sweep(
  rep(TRUE, nrow(meta_orient)), "sun_aspect", aspect_colours,
  plot_prefix = NULL,
  meta = meta_orient, id_col = "sample_id_orient",
  data = samples_lineage_orient)
if (nrow(sweep_aspect) > 0) {
  write.csv(sweep_aspect, "tables/orientation_sun_aspect_multilevel_summary.csv",
            row.names = FALSE)
  print(sweep_aspect)
}

sel_fl <- meta_orient$substrate == "Ficus leaves"
sweep_orient_fl <- run_level_sweep(
  sel_fl, "orientation", orientation_colours,
  plot_prefix = "orientation_ficus_leaves_bylevel",
  meta = meta_orient, id_col = "sample_id_orient",
  data = samples_lineage_orient)
if (nrow(sweep_orient_fl) > 0)
  write.csv(sweep_orient_fl, "tables/orientation_ficus_leaves_multilevel_summary.csv",
            row.names = FALSE)

# ---- G10. Headline summary of every orientation test ----
grab <- function(res, name) {
  if (is.null(res) || is.null(res$permanova))
    return(data.frame(analysis = name, n_units = NA_integer_, n_groups = NA_integer_,
                      permanova_F = NA_real_, permanova_R2 = NA_real_,
                      permanova_p = NA_real_, anosim_R = NA_real_,
                      anosim_p = NA_real_, betadisper_p = NA_real_,
                      n_indicators = NA_integer_, stringsAsFactors = FALSE))
  data.frame(analysis     = name,
             n_units      = res$n_samples,
             n_groups     = res$n_groups,
             permanova_F  = round(res$permanova$F[1], 4),
             permanova_R2 = round(res$permanova$R2[1], 4),
             permanova_p  = res$permanova$`Pr(>F)`[1],
             anosim_R     = if (!is.null(res$anosim_R)) res$anosim_R else NA_real_,
             anosim_p     = if (!is.null(res$anosim_p)) res$anosim_p else NA_real_,
             betadisper_p = if (!is.null(res$betadisper_p)) res$betadisper_p else NA_real_,
             n_indicators = if (!is.null(res$n_indicators)) res$n_indicators else NA_integer_,
             stringsAsFactors = FALSE)
}

orient_summary <- bind_rows(
  grab(results_orient_all, "Orientation - all substrates"),
  bind_rows(lapply(names(results_orient_sub),
                   function(s) grab(results_orient_sub[[s]],
                                    paste0("Orientation - ", s)))),
  grab(results_sub_ori, "Substrate x orientation"),
  grab(results_aspect, "Sun aspect - all substrates"),
  bind_rows(lapply(names(results_aspect_sub),
                   function(s) grab(results_aspect_sub[[s]],
                                    paste0("Sun aspect - ", s))))
)
orient_summary <- orient_summary[!is.na(orient_summary$n_units), , drop = FALSE]
write.csv(orient_summary, "tables/orientation_summary.csv", row.names = FALSE)
cat("\nOrientation analyses - headline summary:\n")
print(orient_summary)

}  # end of orientation section

# ============================================================
# Copy key files for backward compatibility
# ============================================================
for (src_dest in list(
  c("tables/substrate_all_pairwise_permanova.csv", "tables/pairwise_permanova.csv"),
  c("tables/substrate_all_indicator_species.txt",  "tables/indicator_species.txt"),
  c("tables/substrate_all_indicator_species_significant.csv",
    "tables/indicator_species_significant.csv"),
  c("tables/substrate_all_anosim.txt",     "tables/anosim_results.txt"),
  c("tables/substrate_all_betadisper.txt",  "tables/betadisper_results.txt")
)) {
  if (file.exists(src_dest[1]))
    file.copy(src_dest[1], src_dest[2], overwrite = TRUE)
}

# Move plot PDFs from tables/ to plots/community/
for (f in list.files("tables", pattern = "\\.(pdf|png)$", full.names = TRUE)) {
  file.copy(f, file.path("plots/community", basename(f)), overwrite = TRUE)
  file.remove(f)
}

cat("\n====== All community analyses complete ======\n")
