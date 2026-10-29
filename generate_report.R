# =============================================================
# Full results report: EVERY table and plot produced by the pipeline,
# for both taxonomic unit sets plus their side-by-side comparison.
# -> REPORT_all_results.md and REPORT_all_results.html (needs pandoc).
# Images are linked, not embedded: keep the HTML next to results/.
# =============================================================
html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

csv_block <- function(path) {
  # A table with no rows may have been written with no header either
  d <- tryCatch(read.csv(path, check.names = FALSE, stringsAsFactors = FALSE),
                error = function(e) data.frame())
  if (!nrow(d)) return(c("*(empty table)*", ""))
  tab <- paste(knitr::kable(d, format = "html", escape = TRUE,
                            table.attr = 'class="tbl"'), collapse = "\n")
  tab <- paste0('<div class="tblwrap">', tab, "</div>")
  if (nrow(d) > 25)
    tab <- paste0("<details><summary>Show table (", nrow(d), " rows)</summary>\n", tab, "\n</details>")
  c(tab, "")
}
txt_block <- function(path) c(paste0("<pre>", html_escape(paste(readLines(path, warn = FALSE), collapse = "\n")), "</pre>"), "")

files_section <- function(dir, level = "###") {
  out <- character(0)
  tabs <- sort(list.files(file.path(dir, "tables"), pattern = "\\.(csv|txt)$", full.names = TRUE))
  if (length(tabs)) {
    out <- c(out, paste(level, "Tables"), "")
    for (f in tabs)
      out <- c(out, paste0("**`", f, "`**"), "",
               if (grepl("\\.csv$", f)) csv_block(f) else txt_block(f))
  }
  pngs <- sort(list.files(file.path(dir, "plots"), pattern = "\\.png$", recursive = TRUE, full.names = TRUE))
  if (length(pngs)) {
    out <- c(out, paste(level, "Plots"), "")
    for (sub in unique(dirname(pngs))) {
      out <- c(out, paste0(level, "# ", sub), "")
      for (f in pngs[dirname(pngs) == sub])
        out <- c(out, paste0("**", basename(f), "** ([PDF](", sub("\\.png$", ".pdf", f), "))"), "",
                 paste0("![", basename(f), "](", f, ")"), "")
    }
  }
  out
}

ov <- read.csv("results/comparison/tables/unit_set_overview.csv", stringsAsFactors = FALSE)
ov_txt <- paste0(ov$unit_set, ": ", ov$taxa, " taxa", collapse = "; ")

report <- c(
  "---",
  "title: 'LOT2 fungal endophytes (Peru) - all results'",
  paste0("date: '", Sys.Date(), "'"),
  "---",
  "",
  "# About this report",
  "",
  "This report contains **every table and every plot** produced by `main.R`, for each of the two definitions of the taxonomic unit, followed by their side-by-side comparison. Interpretation is left to the reader; see `README.md` for the pipeline and the sampling-design assumptions.",
  "",
  "## Input data",
  "",
  paste0("- `LOT2_samples.xlsx` - one row per isolate (", nrow(samples_raw), " isolates): substrate, tree zone, sampling unit, reconciled orientation, ITS name at **100%** and at **99%** BLAST similarity, GenBank hit. Sheet `Feuil2` holds the ITS sequences."),
  "- `LOT2_pooled_counts_ITS100.xlsx` - per-taxon counts and taxonomy for the 100% names.",
  "- `LOT2_pooled_counts_ITS99.xlsx` - per-taxon counts and taxonomy for the 99% names (**in production**, V. Hofstetter).",
  "",
  "## The two taxonomic unit sets (run in parallel)",
  "",
  "| Unit set | Definition | Source |",
  "|---|---|---|",
  "| **ITS100** (reference) | Taxon named after the GenBank BLAST top score at **100%** ITS similarity, with current names (MycoBank / Index Fungorum) | V. Hofstetter, `LOT2_samples.xlsx` column 9 + `LOT2_pooled_counts_ITS100.xlsx` |",
  "| **ITS99** | Same, at **99%** ITS similarity | V. Hofstetter (in production), column 10 + `LOT2_pooled_counts_ITS99.xlsx` |",
  "",
  paste0("Number of taxa: ", ov_txt, "."),
  "",
  "**Name harmonisation (ITS100 / ITS99).** Isolate names are matched to the workbook names ignoring case, punctuation and the word \"or\", and - failing that - ignoring a trailing \"(...)\" comment. Names with no workbook match are kept as written and take their lineage from the genus in the name. All mismatches are listed under *Data checks*.",
  "",
  "## Changes in this version",
  "",
  "- Order name corrected: *Togniniale* -> **Togniniales**. This and the other trivial corrections made directly in the workbooks are logged in `LOT2_corrections_log.csv`.",
  "- **Correction list for the curator:** `results/data_checks/TO_VERIFY.xlsx` lists every remaining name/lineage inconsistency and flagged value with its cell, plus the applied corrections to confirm.",
  "- New data (re-issued isolate workbook with 100% and 99% names; new 100% and 99% pooled workbooks).",
  "- **97% OTUs (OTU97) are no longer used**: the OTU clustering, its isolate-to-OTU mapping and every analysis based on it were removed.",
  "- **Superclass column removed from the ITS100 and ITS99 workbooks; subclass is used.** Every rank is aggregated on its own value, so a class (e.g. Leotiomycetes, Sordariomycetes, Dothideomycetes, Mortierellomycetes) is never split by inconsistent higher-rank entries; phylum/subphylum are aligned with the class by majority vote (`results/data_checks/lineage_harmonisation.csv`). Species = genus + epithet.",
  "- **All analyses are run for ITS100 and ITS99** (`results/<set>/`), with side-by-side tables and plots in `results/comparison/`.",
  "- **Ficus trunk and branch wood are one substrate (\"Ficus wood\").** The trunk-vs-branch and substrate x position analyses were removed. Kept: the height gradient within Ficus wood (zones 1-6) and the zone-6 (top canopy) substrate comparison.",
  paste0("- **Lifestyle analyses** (pathogenic share, lifestyle composition, Lauraceae = strangled host vs Ficus = strangler): ",
         if (have_lifestyle) "run - see `lifestyle_*` tables and `plots/lifestyle/` in each unit set."
         else "**not run - no `Lifestyle` column was found in any workbook.** Add one (per isolate in `LOT2_samples.xlsx`, or per taxon in the ITS100/ITS99 workbooks) and re-run `main.R`; the analyses then run automatically."),
  "- In every unit set the column/level called `its_taxon` is that set's unit (ITS100 name or ITS99 name).",
  "",
  "# Data checks",
  "",
  unlist(lapply(sort(list.files("results/data_checks", pattern = "\\.csv$", full.names = TRUE)),
                function(f) c(paste0("**`", f, "`**"), "", csv_block(f)))),
  "# Comparison of the two unit sets",
  "",
  files_section("results/comparison", "##"),
  unlist(lapply(names(unit_sets), function(u)
    c(paste0("# Unit set ", u, " - ", unit_sets[[u]]$label), "", files_section(file.path("results", u), "##"))))
)

writeLines(report, "REPORT_all_results.md")
cat("REPORT_all_results.md generated\n")

css <- tempfile(fileext = ".html")
writeLines(c("<style>",
             "body{max-width:1400px;font-family:sans-serif;font-size:14px}",
             "img{max-width:100%;border:1px solid #ddd;margin-bottom:8px}",
             ".tblwrap{overflow-x:auto} .tbl{border-collapse:collapse;font-size:12px}",
             ".tbl td,.tbl th{border:1px solid #ccc;padding:2px 6px}",
             "pre{font-size:11px;background:#f6f6f6;padding:6px;overflow-x:auto}",
             "</style>"), css)
pandoc_bin <- Sys.which("pandoc")
if (!nzchar(pandoc_bin) && requireNamespace("rmarkdown", quietly = TRUE))
  pandoc_bin <- tryCatch(rmarkdown::pandoc_exec(), error = function(e) "")
if (length(pandoc_bin) && nzchar(pandoc_bin)) {
  ok <- system2(pandoc_bin, c("REPORT_all_results.md", "-s", "--toc", "--toc-depth=3",
                              "-H", shQuote(css), "-o", "REPORT_all_results.html"))
  if (ok == 0) cat("REPORT_all_results.html generated\n")
} else {
  cat("WARNING: pandoc not found - only the Markdown report was written\n")
}
