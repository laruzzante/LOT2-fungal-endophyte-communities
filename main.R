library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

# Every PDF figure is also written as PNG, for the HTML report
ggsave <- function(filename, plot = ggplot2::last_plot(), ...) {
  ggplot2::ggsave(filename, plot = plot, ...)
  if (grepl("\\.pdf$", filename))
    ggplot2::ggsave(sub("\\.pdf$", ".png", filename), plot = plot, ..., dpi = 90)
}

source("parse_LOT2.R")

# The whole analysis is run once per taxonomic unit set, each into its
# own results/<set>/{tables,plots} folder
analysis_scripts <- c("plot_abundance.R", "plot_abundance_pooled.R", "plot_pie.R",
                      "community_analysis.R", "diversity_analysis.R", "lifestyle_analysis.R")
project_root <- getwd()
for (unit_set in names(unit_sets)) {
  cat("\n################ Unit set:", unit_set, "-", unit_sets[[unit_set]]$label, "################\n")
  samples_raw <- unit_sets[[unit_set]]$samples
  pooled      <- unit_sets[[unit_set]]$pooled
  unit_label  <- unit_sets[[unit_set]]$label
  out_dir <- file.path(project_root, "results", unit_set)
  unlink(out_dir, recursive = TRUE)
  dir.create(out_dir, recursive = TRUE)
  setwd(out_dir)
  tryCatch(for (f in analysis_scripts) source(file.path(project_root, f)),
           finally = setwd(project_root))
}

source("compare_units.R")
source("generate_report.R")
