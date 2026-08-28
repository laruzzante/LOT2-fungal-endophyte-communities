dir.create("plots/abundance", recursive = TRUE, showWarnings = FALSE)

tax_hierarchy <- c("phylum", "subphylum", "superclass", "class",
                   "subclass", "order", "family", "genus", "species")

dat <- all_isolates %>%
  mutate(across(all_of(tax_hierarchy), ~ replace_na(., "Unknown")),
         its_taxon = replace_na(its_taxon, "Unknown"),
         across(c(n_fungi_laur_leaf, n_fungi_fic_leaf, n_fungi_fic_wood),
                ~ replace_na(., 0)))

# Progressively prepend parent ranks until all labels are unique
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
                .groups = "drop") %>%
      mutate(label = make_unique_labels(., level, parent_cols))
  }

  plot_data <- agg %>%
    select(label, Lauraceae_leaves, Ficus_leaves, Ficus_wood) %>%
    pivot_longer(-label, names_to = "substrate", values_to = "count") %>%
    mutate(substrate = factor(substrate,
                              levels = c("Lauraceae_leaves", "Ficus_leaves", "Ficus_wood"),
                              labels = c("Lauraceae leaves", "Ficus leaves", "Ficus wood")))

  label_order <- plot_data %>%
    group_by(label) %>%
    summarise(total = sum(count), .groups = "drop") %>%
    arrange(total) %>%
    pull(label)
  plot_data$label <- factor(plot_data$label, levels = label_order)

  n_labels <- length(label_order)
  plot_h <- max(6, n_labels * 0.25 + 2)
  x_size <- if (n_labels > 80) 4 else if (n_labels > 40) 6 else 8

  p <- ggplot(plot_data, aes(x = label, y = count, fill = substrate)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = substrate_colours) +
    coord_flip() +
    labs(title = paste("Fungal isolate abundance by", level),
         x = NULL, y = "Number of isolates", fill = "Substrate") +
    theme_minimal(base_size = 12) +
    theme(axis.text.y  = element_text(size = x_size),
          legend.position = "top",
          plot.title  = element_text(face = "bold"))

  ggsave(file.path("plots/abundance", paste0("abundance_by_", level, ".pdf")),
         plot = p, width = 12, height = plot_h, limitsize = FALSE)
  cat("Saved: abundance_by_", level, ".pdf  (", n_labels, " groups)\n", sep = "")
}
