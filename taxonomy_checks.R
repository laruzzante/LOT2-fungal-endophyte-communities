# =============================================================
# Correction list for the data curator
#
# Every taxonomic-name inconsistency within and between the input
# workbooks, and every value flagged as needing verification, one row
# per issue with the cell to edit:
#   results/data_checks/TO_VERIFY.csv
#   results/data_checks/TO_VERIFY.xlsx  (filterable; needs openxlsx)
# The pipeline does not change any of these values; it only lists them.
# Trivial corrections already made in the workbooks are recorded in
# LOT2_corrections_log.csv and listed too, as "Applied - please confirm".
# Sourced by main.R after parse_LOT2.R.
# =============================================================

rank_cols <- c(phylum = "Phylum", subphylum = "Subphylum", class = "Class",
               subclass = "Subclass", order = "Order", family = "Family",
               genus = "Genus", species = "Species")
# Standard ending of a name at each rank (ICN); a value without it is
# either misspelled or a name of another rank entered in this column
rank_suffix <- c(Phylum = "mycota$", Subphylum = "mycotina$", Class = "mycetes$",
                 Subclass = "mycetidae$", Order = "ales$", Family = "aceae$")
is_unresolved <- function(x) is.na(x) | squish(x) %in% uncertain_values |
  tolower(squish(x)) == "incertae sedis"
col_letter <- function(i) {
  out <- character(length(i))
  for (k in seq_along(i)) {
    n <- i[k]; s <- ""
    while (n > 0) { r <- (n - 1) %% 26; s <- paste0(LETTERS[r + 1], s); n <- (n - 1) %/% 26 }
    out[k] <- s
  }
  out
}

issues <- list()
add_issue <- function(priority, check, workbook, sheet, row, column, culture_code,
                      taxon, field, current, issue, suggestion = NA_character_) {
  if (!length(row)) return(invisible())
  column <- rep_len(column, length(row))
  issues[[length(issues) + 1]] <<- data.frame(
    priority = priority, check = check, workbook = workbook, sheet = sheet,
    cell = ifelse(is.na(column), NA_character_, paste0(column, row)),
    row = row, culture_code = culture_code, taxon_name = taxon, field = field,
    current_value = current, issue = issue, suggested_value = suggestion,
    stringsAsFactors = FALSE)
}

# --- Raw taxonomy of the pooled workbooks, with cell positions ---
raw_wb <- list()
for (s in c("ITS100", "ITS99")) {
  f <- paste0("LOT2_pooled_counts_", s, ".xlsx")
  p <- read_excel(f, sheet = 1)
  sheet <- excel_sheets(f)[1]
  d <- data.frame(set = s, workbook = f, sheet = sheet, row = seq_len(nrow(p)) + 1,
                  code = squish(p[[1]]), name = squish(p[[2]]), stringsAsFactors = FALSE)
  for (r in rank_cols) d[[r]] <- squish(p[[r]])
  d <- d[!(is.na(d$code) & is.na(d$name)), ]
  attr(d, "col") <- setNames(col_letter(match(rank_cols, names(p))), rank_cols)
  attr(d, "name_col") <- col_letter(2)
  raw_wb[[s]] <- d
}
wb_all <- bind_rows(raw_wb)
col_of <- function(set, field) attr(raw_wb[[set]], "col")[[field]]

# 1. Value does not carry the ending of its rank
for (r in names(rank_suffix)) {
  v <- wb_all[[r]]
  bad <- !is_unresolved(v) & !grepl(rank_suffix[[r]], v)
  for (i in which(bad)) {
    # Suggest a correctly-ended value at this rank that is one or two
    # letters away (e.g. Mucoromycete -> Mucoromycetes)
    ok <- unique(v[!is_unresolved(v) & grepl(rank_suffix[[r]], v)])
    near <- ok[adist(v[i], ok, ignore.case = TRUE)[1, ] <= 2]
    add_issue("High", "Wrong ending for the rank (typo, or a name of another rank)",
              wb_all$workbook[i], wb_all$sheet[i], wb_all$row[i], col_of(wb_all$set[i], r),
              wb_all$code[i], wb_all$name[i], r, v[i],
              paste0(r, " names end in -", sub("\\$$", "", rank_suffix[[r]])),
              if (length(near)) near[1] else NA_character_)
  }
}

# 2. Near-identical spellings of one value at the same rank (both workbooks
#    pooled, so a spelling that differs between them is caught too)
for (r in rank_cols[1:7]) {
  tab <- table(wb_all[[r]][!is_unresolved(wb_all[[r]])])
  v <- names(tab)
  if (length(v) < 2) next
  m <- adist(v, ignore.case = TRUE)
  for (a in seq_along(v)) for (b in seq_along(v)) {
    if (b <= a || m[a, b] > 2 || substr(v[a], 1, 2) != substr(v[b], 1, 2)) next
    minor <- if (tab[[v[a]]] < tab[[v[b]]]) v[a] else v[b]
    major <- setdiff(c(v[a], v[b]), minor)
    rows <- which(wb_all[[r]] == minor)
    add_issue("High", "Two spellings of the same name", wb_all$workbook[rows], wb_all$sheet[rows],
              wb_all$row[rows], sapply(wb_all$set[rows], col_of, field = r),
              wb_all$code[rows], wb_all$name[rows], r, minor,
              paste0("also spelled '", major, "' (", tab[[major]], " rows)"), major)
  }
}

# 3. A name placed under different parents (e.g. one class under two
#    phyla). The minority placement is listed against the majority one.
parent_pairs <- list(c("Genus", "Family"), c("Family", "Order"), c("Order", "Subclass"),
                     c("Order", "Class"), c("Subclass", "Class"),
                     c("Class", "Subphylum"), c("Class", "Phylum"), c("Subphylum", "Phylum"))
for (pp in parent_pairs) {
  ch <- pp[1]; pa <- pp[2]
  d <- wb_all[!is_unresolved(wb_all[[ch]]) & !is_unresolved(wb_all[[pa]]), ]
  for (g in unique(d[[ch]])) {
    x <- d[d[[ch]] == g, ]
    if (n_distinct(x[[pa]]) < 2) next
    major <- names(sort(table(x[[pa]]), decreasing = TRUE))[1]
    rows <- which(x[[pa]] != major)
    add_issue("High", paste(ch, "placed under two different", pa, "names"),
              x$workbook[rows], x$sheet[rows], x$row[rows], sapply(x$set[rows], col_of, field = pa),
              x$code[rows], x$name[rows], pa, x[[pa]][rows],
              paste0(ch, " ", g, " is under ", pa, " '", major, "' in ",
                     sum(x[[pa]] == major), " other row(s)"), major)
  }
}

# 4. Unresolved rank ("?", "NO", "incertae sedis", blank) that other rows
#    resolve for the same lower-rank name - fillable from the workbooks
fill_from <- list(Phylum = "Class", Subphylum = "Class", Class = "Order", Subclass = "Order",
                  Order = "Family", Family = "Genus")
for (r in names(fill_from)) {
  ch <- fill_from[[r]]
  known <- wb_all[!is_unresolved(wb_all[[ch]]) & !is_unresolved(wb_all[[r]]), ] %>%
    group_by(child = .data[[ch]]) %>%
    summarise(value = names(sort(table(.data[[r]]), decreasing = TRUE))[1],
              n = n(), .groups = "drop")
  gap <- which(is_unresolved(wb_all[[r]]) & !is_unresolved(wb_all[[ch]]) &
                 wb_all[[ch]] %in% known$child)
  k <- match(wb_all[[ch]][gap], known$child)
  add_issue("Medium", "Unresolved rank that other rows resolve",
            wb_all$workbook[gap], wb_all$sheet[gap], wb_all$row[gap],
            sapply(wb_all$set[gap], col_of, field = r), wb_all$code[gap], wb_all$name[gap], r,
            ifelse(is.na(wb_all[[r]][gap]), "(blank)", wb_all[[r]][gap]),
            paste0(ch, " ", wb_all[[ch]][gap], " has ", r, " '", known$value[k], "' in ",
                   known$n[k], " other row(s)"), known$value[k])
}
# Genus unresolved although the name starts with a genus used elsewhere
genera <- unique(wb_all$Genus[!is_unresolved(wb_all$Genus)])
tok <- sub(" .*", "", wb_all$name)
gap <- which(is_unresolved(wb_all$Genus) & tok %in% genera)
add_issue("Medium", "Unresolved rank that other rows resolve",
          wb_all$workbook[gap], wb_all$sheet[gap], wb_all$row[gap],
          sapply(wb_all$set[gap], col_of, field = "Genus"), wb_all$code[gap], wb_all$name[gap],
          "Genus", wb_all$Genus[gap], "the name starts with a genus used in other rows", tok[gap])

# 5. Genus / Species columns that disagree with the name itself
binomial <- grepl("^[A-Z][a-z]+ [a-z-]+$", wb_all$name) &
  !grepl(" (sp|spp|cf|aff|or|gen)\\.?$", wb_all$name)
epithet  <- sub("^\\S+ ", "", wb_all$name)
bad <- which(binomial & (is_unresolved(wb_all$Species) | wb_all$Species != epithet))
# A missing epithet is filled from the name; a different epithet means the
# name or the Species column is wrong, which only the curator can tell
add_issue("Medium", "Species column disagrees with the name",
          wb_all$workbook[bad], wb_all$sheet[bad], wb_all$row[bad],
          sapply(wb_all$set[bad], col_of, field = "Species"), wb_all$code[bad], wb_all$name[bad],
          "Species", wb_all$Species[bad],
          ifelse(is_unresolved(wb_all$Species[bad]), "the name is a binomial",
                 paste0("the name has epithet '", epithet[bad], "'; correct the name or this column")),
          ifelse(is_unresolved(wb_all$Species[bad]), epithet[bad], NA_character_))
bad <- which(!is_unresolved(wb_all$Genus) & tok != wb_all$Genus)
add_issue("High", "Genus column disagrees with the name",
          wb_all$workbook[bad], wb_all$sheet[bad], wb_all$row[bad],
          sapply(wb_all$set[bad], col_of, field = "Genus"), wb_all$code[bad], wb_all$name[bad],
          "Genus", wb_all$Genus[bad], "the name starts with another genus", tok[bad])

# 6. The same taxon listed twice in one workbook
for (s in names(raw_wb)) {
  d <- raw_wb[[s]]; k <- name_key(d$name)
  dup <- which(!is.na(d$name) & (duplicated(k) | duplicated(k, fromLast = TRUE)))
  add_issue("High", "Taxon listed more than once in the workbook", d$workbook[dup], d$sheet[dup],
            d$row[dup], attr(d, "name_col"), d$code[dup], d$name[dup], "ITS name", d$name[dup],
            "same name (ignoring case and punctuation) on another row; merge the counts")
}

# 7. The same taxon name with a different lineage in the two workbooks
k100 <- name_key(raw_wb$ITS100$name); k99 <- name_key(raw_wb$ITS99$name)
for (i in which(k100 %in% k99)) {
  j <- match(k100[i], k99)
  for (r in rank_cols) {
    a <- raw_wb$ITS100[[r]][i]; b <- raw_wb$ITS99[[r]][j]
    if (identical(coalesce(a, ""), coalesce(b, ""))) next
    add_issue("High", "Same name, different lineage in the ITS100 and ITS99 workbooks",
              raw_wb$ITS99$workbook[j], raw_wb$ITS99$sheet[j], raw_wb$ITS99$row[j],
              col_of("ITS99", r), raw_wb$ITS99$code[j], raw_wb$ITS99$name[j], r,
              coalesce(b, "(blank)"),
              paste0("ITS100 workbook row ", raw_wb$ITS100$row[i], " has '", coalesce(a, "(blank)"), "'"),
              if (is_unresolved(a)) NA_character_ else a)
  }
}

# 8. Names carrying a query or a note to resolve
for (s in names(raw_wb)) {
  d <- raw_wb[[s]]
  q <- which(grepl("\\?", d$name))
  add_issue("Medium", "Name contains a query", d$workbook[q], d$sheet[q], d$row[q],
            attr(d, "name_col"), d$code[q], d$name[q], "ITS name", d$name[q],
            "the name contains '?' - a note to resolve before publication")
}

# --- Isolate workbook (LOT2_samples.xlsx, first sheet) --------------
samp_file  <- "LOT2_samples.xlsx"
samp_sheet <- excel_sheets(samp_file)[1]
iso <- unit_sets$ITS100$samples
iso$xl_row <- seq_len(nrow(iso)) + 1     # the sheet has no blank rows

# 9. Isolate names absent from the pooled workbook of their set
for (s in names(raw_wb)) {
  raw_col <- if (s == "ITS100") "its100_raw" else "its99_raw"
  xl_col  <- if (s == "ITS100") "I" else "J"
  idx <- map_to_workbook(squish(iso[[raw_col]]), pooled_wb[[s]])
  m <- which(is.na(idx) & !is.na(iso[[raw_col]]))
  add_issue("High", paste("Isolate name not found in the", s, "workbook"),
            samp_file, samp_sheet, iso$xl_row[m], xl_col, iso$culture_code1[m],
            squish(iso[[raw_col]][m]), paste("ITS name", s), squish(iso[[raw_col]][m]),
            paste0("no row of LOT2_pooled_counts_", s, ".xlsx carries this name; ",
                   "spelling differs or the taxon is missing there"))
}

# 10. Per-taxon counts that differ between a pooled workbook and the isolates
for (s in names(raw_wb)) {
  f <- file.path("results/data_checks", paste0(s, "_workbook_vs_isolate_counts.csv"))
  cmp <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  if (!nrow(cmp)) next
  d <- raw_wb[[s]]
  r <- d$row[match(name_key(cmp$its_taxon), name_key(d$name))]
  txt <- sprintf("workbook %d / %d / %d vs isolate records %d / %d / %d (Lauraceae leaves / Ficus leaves / Ficus wood)",
                 cmp$wb_lauraceae, cmp$wb_ficus_leaves, cmp$wb_ficus_wood,
                 cmp$iso_lauraceae, cmp$iso_ficus_leaves, cmp$iso_ficus_wood)
  add_issue("High", paste("Counts differ between the", s, "workbook and the isolate records"),
            ifelse(is.na(r), samp_file, d$workbook[1]), ifelse(is.na(r), samp_sheet, d$sheet[1]),
            r, ifelse(is.na(r), NA, "C"), NA, cmp$its_taxon, "counts (columns C-E)",
            NA, ifelse(is.na(r), paste(txt, "- taxon only in the isolate records"), txt))
}

# 11. One isolate, two lineages: its 100% and 99% names fall in different
#     taxa at genus rank or above
s99 <- unit_sets$ITS99$samples
for (r in c("class", "order", "family", "genus")) {
  a <- iso[[r]]; b <- s99[[r]]
  bad <- which(!is_unresolved(a) & !is_unresolved(b) & a != b)
  add_issue("High", "ITS100 and ITS99 names of one isolate disagree on the lineage",
            samp_file, samp_sheet, iso$xl_row[bad], "J", iso$culture_code1[bad],
            squish(iso$its100_raw[bad]), r, b[bad],
            paste0(rank_cols[[r]], ": '", a[bad], "' from the 100% name (column I) vs '", b[bad],
                   "' from the 99% name (column J) - ", squish(s99$its99_raw[bad])))
}

# 12. Corrections the pipeline applies on reading (to make in the workbook)
raw_unit <- read_excel(samp_file, sheet = 1, .name_repair = "minimal")[[3]]
fix <- which(iso$substrate == "Ficus wood" & iso$zone %in% 1:5 & squish(raw_unit) != iso$unit)
add_issue("Medium", "Sampling unit corrected by the pipeline on reading",
          samp_file, samp_sheet, iso$xl_row[fix], "C", iso$culture_code1[fix],
          squish(iso$its100_raw[fix]), "Sampling unit", squish(raw_unit[fix]),
          "Ficus trunk units are S11-S15, one per zone (5 to 1)", iso$unit[fix])

# 13. Orientation not reconciled
unres <- which(is.na(iso$orientation))
add_issue("Low", "Orientation unresolved", samp_file, samp_sheet, iso$xl_row[unres], "G",
          iso$culture_code1[unres], squish(iso$its100_raw[unres]), "Orientation (the one to use)",
          ifelse(is.na(iso$orientation_raw[unres]), "(blank)", squish(iso$orientation_raw[unres])),
          "excluded from the orientation analyses until resolved")

# --- Write ------------------------------------------------------
to_verify <- bind_rows(issues) %>% mutate(status = "To fix")

log_file <- "LOT2_corrections_log.csv"
if (file.exists(log_file)) {
  cl <- read.csv(log_file, stringsAsFactors = FALSE, check.names = FALSE)
  to_verify <- bind_rows(to_verify, data.frame(
    priority = "Applied", check = paste("Corrected in the workbook on", cl$date),
    workbook = cl$workbook, sheet = cl$sheet, cell = cl$cell,
    row = suppressWarnings(as.integer(sub("^[A-Z]+", "", cl$cell))),
    field = cl$field, current_value = cl$new_value, issue = cl$reason,
    suggested_value = paste("was:", cl$old_value),
    status = "Applied - please confirm", stringsAsFactors = FALSE))
}
to_verify <- to_verify %>%
  mutate(priority = factor(priority, levels = c("High", "Medium", "Low", "Applied"))) %>%
  arrange(priority, check, workbook, row) %>%
  mutate(priority = as.character(priority), done = "") %>%
  relocate(status, .after = priority)
write.csv(to_verify, "results/data_checks/TO_VERIFY.csv", row.names = FALSE, na = "")

check_summary <- to_verify %>% count(priority, check, name = "rows") %>%
  mutate(priority = factor(priority, levels = c("High", "Medium", "Low", "Applied"))) %>%
  arrange(priority, desc(rows)) %>% mutate(priority = as.character(priority))
write.csv(check_summary, "results/data_checks/TO_VERIFY_summary.csv", row.names = FALSE)

if (requireNamespace("openxlsx", quietly = TRUE)) {
  xw <- openxlsx::createWorkbook()
  hs <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#DDDDDD", wrapText = TRUE)
  openxlsx::addWorksheet(xw, "Summary")
  openxlsx::writeData(xw, "Summary", check_summary, headerStyle = hs)
  openxlsx::setColWidths(xw, "Summary", 1:3, c(10, 70, 8))
  openxlsx::addWorksheet(xw, "To verify")
  openxlsx::writeData(xw, "To verify", to_verify, headerStyle = hs, withFilter = TRUE, keepNA = FALSE)
  openxlsx::freezePane(xw, "To verify", firstRow = TRUE)
  openxlsx::setColWidths(xw, "To verify", seq_along(to_verify),
                         c(9, 22, 40, 30, 14, 7, 6, 12, 45, 14, 25, 60, 25, 7))
  openxlsx::saveWorkbook(xw, "results/data_checks/TO_VERIFY.xlsx", overwrite = TRUE)
}
cat("Correction list:", nrow(to_verify), "item(s) -> results/data_checks/TO_VERIFY.xlsx\n")
print(as.data.frame(check_summary), right = FALSE)
