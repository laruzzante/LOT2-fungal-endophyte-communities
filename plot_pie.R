piedir <- "plots/pie_charts"
dir.create(piedir, recursive = TRUE, showWarnings = FALSE)

tax_hierarchy <- c("phylum", "subphylum", "superclass", "class",
                   "subclass", "order", "family", "genus", "species")

uncertain_values <- c("?", "NO", "incertae sedis")
incertae_label   <- "Incertae sedis"

dat_pie <- pooled %>%
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

build_palette <- function(n) {
  if (n <= 8) {
    pal <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
             "#FF7F00", "#A65628", "#F781BF", "#999999")
  } else {
    base <- c(
      "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628",
      "#F781BF","#999999","#66C2A5","#FC8D62","#8DA0CB","#E78AC3",
      "#A6D854","#FFD92F","#E5C494","#B3B3B3","#1B9E77","#D95F02",
      "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666",
      "#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462",
      "#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F"
    )
    if (n <= length(base)) {
      pal <- base[1:n]
    } else {
      pal <- colorRampPalette(base)(n)
    }
  }
  pal[1:n]
}

substrate_levels <- c("Lauraceae_leaves", "Ficus_leaves", "Ficus_wood")
substrate_labels <- c("Lauraceae leaves", "Ficus leaves", "Ficus wood")

for (level in c("its_taxon", tax_hierarchy)) {

  if (level == "its_taxon") {
    agg <- dat_pie %>%
      group_by(its_taxon) %>%
      summarise(Lauraceae_leaves = sum(n_fungi_laur_leaf),
                Ficus_leaves     = sum(n_fungi_fic_leaf),
                Ficus_wood       = sum(n_fungi_fic_wood),
                .groups = "drop") %>%
      mutate(label = as.character(its_taxon))
  } else {

  idx <- which(tax_hierarchy == level)
  group_cols  <- tax_hierarchy[1:idx]
  parent_cols <- if (idx >= 2) tax_hierarchy[1:(idx - 1)] else character(0)

  agg <- dat_pie %>%
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
    filter(count > 0) %>%
    mutate(substrate = factor(substrate,
                              levels = substrate_levels,
                              labels = substrate_labels))

  label_order <- plot_data %>%
    group_by(label) %>%
    summarise(total = sum(count), .groups = "drop") %>%
    arrange(desc(total)) %>%
    pull(label)
  plot_data$label <- factor(plot_data$label, levels = label_order)

  n_labels <- length(label_order)
  taxon_colours <- setNames(build_palette(n_labels), label_order)

  plot_data <- plot_data %>%
    group_by(substrate) %>%
    mutate(pct = count / sum(count) * 100) %>%
    ungroup()

  plot_data <- plot_data %>%
    mutate(pie_label = ifelse(pct >= 3, paste0(round(pct, 1), "%"), ""))

  legend_cols <- min(4, ceiling(n_labels / 15))
  legend_size <- if (n_labels > 50) 5 else if (n_labels > 25) 6 else 8
  plot_w <- if (n_labels > 30) 18 else 14
  plot_h <- max(7, 5 + ceiling(n_labels / legend_cols) * 0.22)

  p <- ggplot(plot_data, aes(x = "", y = count, fill = label)) +
    geom_col(width = 1, colour = "white", linewidth = 0.3) +
    geom_text(aes(label = pie_label),
              position = position_stack(vjust = 0.5), size = 2.5) +
    coord_polar(theta = "y") +
    facet_wrap(~ substrate, nrow = 1, scales = "free_y") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = taxon_colours, name = level) +
    labs(title = paste("Fungal isolate composition by", level),
         subtitle = "Unknown/missing taxa pooled as Incertae sedis") +
    theme_void(base_size = 12) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          strip.text    = element_text(face = "bold", size = 12),
          legend.position = "bottom",
          legend.text   = element_text(size = legend_size),
          legend.title  = element_text(face = "bold"),
          legend.key.size = unit(0.35, "cm"),
          plot.margin   = margin(10, 10, 10, 10))

  if (n_labels > 60) {
    p <- p + guides(fill = guide_legend(ncol = legend_cols, override.aes = list(size = 1)))
  } else {
    p <- p + guides(fill = guide_legend(ncol = legend_cols))
  }

  ggsave(file.path(piedir, paste0("pie_by_", level, ".pdf")),
         plot = p, width = plot_w, height = plot_h, limitsize = FALSE)
  cat("Saved: ", file.path(piedir, paste0("pie_by_", level, ".pdf")),
      "  (", n_labels, " taxa)\n", sep = "")
}
