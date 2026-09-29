outdir <- "plots/abundance_pooled_incertae_sedis"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

tax_hierarchy <- c("phylum", "subphylum", "superclass", "class",
                   "subclass", "order", "family", "genus", "species")

uncertain_values <- c("?", "NO", "incertae sedis")
incertae_label   <- "Incertae sedis"

dat <- pooled %>%
  mutate(across(all_of(tax_hierarchy),
                ~ ifelse(is.na(.) | . %in% uncertain_values,
                         incertae_label, .)),
         its_taxon = ifelse(is.na(its_taxon) | its_taxon %in% uncertain_values,
                            incertae_label, its_taxon),
         across(c(n_fungi_laur_leaf, n_fungi_fic_leaf, n_fungi_fic_wood),
                ~ replace_na(., 0)))

make_unique_labels <- function(df, target_col, parent_cols) {
  labels <- df[[target_col]]
  for (pcol in rev(parent_cols)) {
    dups <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
    if (!any(dups)) break
    labels[dups] <- paste(df[[pcol]][dups], labels[dups], sep = " | ")
  }
  labels
}

substrate_colours <- c("Lauraceae leaves" = "#2E86AB",
                       "Ficus leaves"     = "#A23B72",
                       "Ficus wood"       = "#F18F01")

levels_to_plot <- c("culture_code", "its_taxon", tax_hierarchy)

for (level in levels_to_plot) {

  if (level %in% c("culture_code", "its_taxon")) {
    agg <- dat %>%
      group_by(.data[[level]]) %>%
      summarise(Lauraceae_leaves = sum(n_fungi_laur_leaf),
                Ficus_leaves     = sum(n_fungi_fic_leaf),
                Ficus_wood       = sum(n_fungi_fic_wood),
                .groups = "drop") %>%
      mutate(label = as.character(.data[[level]]))
  } else {
    idx <- which(tax_hierarchy == level)
    group_cols  <- tax_hierarchy[1:idx]
    parent_cols <- if (idx >= 2) tax_hierarchy[1:(idx - 1)] else character(0)

    agg <- dat %>%
      group_by(across(all_of(group_cols))) %>%
      summarise(Lauraceae_leaves = sum(n_fungi_laur_leaf),
                Ficus_leaves     = sum(n_fungi_fic_leaf),
                Ficus_wood       = sum(n_fungi_fic_wood),
                .groups = "drop")

    is_uncertain <- agg[[level]] == incertae_label
    if (any(is_uncertain)) {
      certain   <- agg[!is_uncertain, ]
      uncertain <- agg[is_uncertain, ] %>%
        summarise(across(c(Lauraceae_leaves, Ficus_leaves, Ficus_wood), sum))
      for (gc in group_cols) uncertain[[gc]] <- incertae_label
      agg <- bind_rows(certain, uncertain)
    }

    agg <- agg %>%
      mutate(label = make_unique_labels(., level, parent_cols))
  }

  plot_data <- agg %>%
    select(label, Lauraceae_leaves, Ficus_leaves, Ficus_wood) %>%
    pivot_longer(-label, names_to = "substrate", values_to = "count") %>%
    mutate(substrate = factor(substrate,
                              levels = c("Lauraceae_leaves", "Ficus_leaves", "Ficus_wood"),
                              labels = c("Lauraceae leaves", "Ficus leaves", "Ficus wood")))

  # Identified taxa ranked by abundance; Incertae sedis pinned to the top
  # row (last level under coord_flip) so it is never lost in the ranking
  label_order <- plot_data %>%
    group_by(label) %>%
    summarise(total = sum(count), .groups = "drop") %>%
    arrange(total) %>%
    pull(label)
  has_incertae <- incertae_label %in% label_order
  label_order <- c(setdiff(label_order, incertae_label),
                   if (has_incertae) incertae_label)
  plot_data$label <- factor(plot_data$label, levels = label_order)

  # Unresolved share of each substrate's isolates
  plot_data <- plot_data %>%
    group_by(substrate) %>%
    mutate(pct = 100 * count / sum(count)) %>%
    ungroup()
  # The dodged bars are too thin to carry readable value labels, so the
  # per-substrate shares go in the subtitle, which is always legible
  incertae_rows <- plot_data %>% filter(label == incertae_label)
  subtitle <- if (has_incertae)
    paste0(incertae_label, " (grey band, top row) = unknown or unresolved at this rank.\n",
           "Share of isolates: ",
           paste0(incertae_rows$substrate, " ", sprintf("%.1f", incertae_rows$pct), "%",
                  " (", incertae_rows$count, ")", collapse = ", "))
  else paste0("No isolates are ", incertae_label, " at this rank")

  n_labels <- length(label_order)
  plot_h <- max(6, n_labels * 0.25 + 2)
  x_size <- if (n_labels > 80) 4 else if (n_labels > 40) 6 else 8
  y_max  <- max(plot_data$count)

  p <- ggplot(plot_data, aes(x = label, y = count, fill = substrate))
  if (has_incertae)
    p <- p + geom_col(data = data.frame(label = factor(incertae_label, levels = label_order),
                                        y = y_max * 1.05),
                      aes(x = label, y = y), inherit.aes = FALSE,
                      width = 1, fill = "#EEEEEE")
  p <- p +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = substrate_colours) +
    # explicit limits: the band layer would otherwise reorder the axis
    scale_x_discrete(limits = label_order) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    coord_flip() +
    labs(title = paste("Fungal isolate abundance by", level),
         subtitle = subtitle,
         x = NULL, y = "Number of isolates", fill = "Substrate") +
    theme_minimal(base_size = 12) +
    theme(axis.text.y  = element_text(size = x_size),
          legend.position = "top",
          plot.title.position = "plot",
          plot.title  = element_text(face = "bold"))

  ggsave(file.path(outdir, paste0("abundance_by_", level, ".pdf")),
         plot = p, width = 12, height = plot_h, limitsize = FALSE)
  cat("Saved: ", file.path(outdir, paste0("abundance_by_", level, ".pdf")),
      "  (", n_labels, " groups)\n", sep = "")
}
