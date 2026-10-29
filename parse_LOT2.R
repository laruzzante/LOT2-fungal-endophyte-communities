# =============================================================
# Parse LOT2 input data
#
# Inputs (first sheet of each workbook is authoritative):
#   LOT2_samples.xlsx              one row per isolate, with the ITS name
#                                  at 100% AND at 99% BLAST similarity
#   LOT2_pooled_counts_ITS100.xlsx per-taxon counts + taxonomy, 100% names
#   LOT2_pooled_counts_ITS99.xlsx  per-taxon counts + taxonomy, 99% names
#                                  (in production - V. Hofstetter)
# =============================================================

# The superclass column was removed from the workbooks (it was
# inconsistently filled); subclass is consistent and is used.
tax_hierarchy    <- c("phylum", "subphylum", "class", "subclass", "order",
                      "family", "genus", "species")
uncertain_values <- c("?", "NO", "incertae sedis", "")
incertae_label   <- "Incertae sedis"

squish    <- function(x) trimws(gsub("\\s+", " ", as.character(x)))
clean_tax <- function(x) { x <- squish(x); ifelse(is.na(x) | x %in% uncertain_values, NA_character_, x) }
majority  <- function(x) {
  x <- clean_tax(x); x <- x[!is.na(x)]
  if (!length(x)) NA_character_ else names(sort(table(x), decreasing = TRUE))[1]
}
# Spelling-insensitive key: "A, or B" == "A or B" == "A B"
name_key <- function(x) gsub("[^a-z0-9]", "", gsub("\\bor\\b", "", tolower(squish(x))))
# Same, ignoring a trailing "(...)" comment, e.g. "Diaporthe sp. 7 (>99% similarity with sp. 30)"
name_key_stripped <- function(x) name_key(sub("\\s*\\([^()]*\\)\\s*$", "", squish(x)))

dir.create("results/data_checks", recursive = TRUE, showWarnings = FALSE)

# --- Pooled per-taxon workbooks -------------------------------
read_pooled <- function(path) {
  p <- read_excel(path, sheet = 1)
  ls_col <- grep("^\\s*life\\s*style", names(p), ignore.case = TRUE)
  out <- tibble(
    culture_code      = squish(p[[1]]),
    its_taxon         = squish(p[[2]]),
    n_fungi_laur_leaf = as.numeric(p[[3]]),
    n_fungi_fic_leaf  = as.numeric(p[[4]]),
    n_fungi_fic_wood  = as.numeric(p[[5]]),
    phylum = p$Phylum, subphylum = p$Subphylum, class = p$Class,
    subclass = p$Subclass, order = p$Order, family = p$Family, genus = p$Genus,
    species_epithet = p$Species,
    lifestyle = if (length(ls_col)) squish(p[[ls_col[1]]]) else NA_character_)
  # The sheet ends with a TOTALS row (no code, no name); read as a taxon
  # it would become a phantom slice holding a whole substrate's isolates.
  n0 <- nrow(out)
  out <- out %>% filter(!(is.na(culture_code) & is.na(its_taxon)))
  if (n0 > nrow(out)) cat(basename(path), ": dropped", n0 - nrow(out), "total row(s)\n")
  out %>% mutate(
    across(c(phylum, subphylum, class, subclass, order, family, genus), clean_tax),
    species = ifelse(!is.na(genus) & !is.na(clean_tax(species_epithet)),
                     paste(genus, clean_tax(species_epithet)), NA_character_),
    lifestyle = clean_tax(lifestyle)) %>%
    select(-species_epithet)
}

# Higher ranks are made consistent with the class (majority vote), so a
# class is never split by a mis-typed phylum/subphylum (e.g. one
# Mortierellomycetes row with phylum "incertae sedis").
harmonise_lineage <- function(p, set_name) {
  changes <- list()
  for (parent in c("phylum", "subphylum")) {
    ref <- p %>% filter(!is.na(class)) %>% group_by(class) %>%
      summarise(ref = majority(.data[[parent]]), .groups = "drop")
    new <- ref$ref[match(p$class, ref$class)]
    chg <- !is.na(new) & (is.na(p[[parent]]) | p[[parent]] != new)
    if (any(chg))
      changes[[parent]] <- data.frame(unit_set = set_name, its_taxon = p$its_taxon[chg],
                                      class = p$class[chg], rank = parent,
                                      original = p[[parent]][chg], harmonised = new[chg])
    p[[parent]][chg] <- new[chg]
  }
  list(data = p, changes = bind_rows(changes))
}

pooled_wb <- list()
lineage_changes <- list()
for (s in c("ITS100", "ITS99")) {
  h <- harmonise_lineage(read_pooled(paste0("LOT2_pooled_counts_", s, ".xlsx")), s)
  pooled_wb[[s]] <- h$data
  lineage_changes[[s]] <- h$changes
}
lineage_changes <- bind_rows(lineage_changes)
write.csv(lineage_changes, "results/data_checks/lineage_harmonisation.csv", row.names = FALSE)
if (nrow(lineage_changes))
  cat("Lineage harmonisation:", nrow(lineage_changes), "phylum/subphylum value(s) aligned with their class\n")

# --- Sample-level isolate records (first sheet only) ---
# Layout (12 columns):
#   1 substrate              5 branch colour (project database)   9 ITS name, 100% similarity
#   2 tree zone              6 orientation (pre-reconciliation)  10 ITS name, 99% similarity
#   3 sampling unit          7 orientation (reconciled - USE)    11 GenBank hit
#   4 branch colour (field)  8 culture code                      12 culture code (dup)
# Columns 4-6 are reconciliation working columns kept for traceability;
# only column 7 is authoritative for orientation.
samples_raw <- read_excel("LOT2_samples.xlsx", sheet = 1, .name_repair = "minimal")
iso_ls_col  <- grep("^\\s*life\\s*style", names(samples_raw), ignore.case = TRUE)
lifestyle_isolate <- if (length(iso_ls_col)) clean_tax(samples_raw[[iso_ls_col[1]]]) else NA_character_
samples_raw <- samples_raw[, 1:12]
names(samples_raw) <- c("substrate", "zone", "unit",
                        "branch_colour_field", "branch_colour_db",
                        "orientation_old", "orientation_raw",
                        "culture_code1", "its100_raw", "its99_raw", "genbank",
                        "culture_code2")
samples_raw$lifestyle_isolate <- lifestyle_isolate

# Drop any fully empty trailing rows
samples_raw <- samples_raw[!is.na(samples_raw$substrate), ]

# --- Sampling-unit corrections from the data curator ----------
# Two labelling errors in the workbook, both confirmed by the data
# curator with the field team (LoT). They are corrected here, by rule,
# so the delivered workbook stays untouched and each correction becomes
# a no-op once a fixed workbook is issued.
#
# 1. Every Lauraceae branch was collected in ZONE 6. Some rows carry
#    zone 5 (all of S7, most of S8, one S4 isolate), copied from a field
#    file that later proved wrong. A sampling unit is substrate x zone x
#    unit, so this split S4 and S8 in two and turned the 8 Lauraceae
#    branches (S1-S8) into 10 units.
# 2. Ficus trunk wood is numbered S11-S15 (zone 5 down to zone 1). It was
#    renamed from S9-S13 so that it cannot be mistaken for the wood of
#    canopy branches S9/S10, whose leaves only were sampled. The rename
#    missed 3 zone-4 isolates still labelled S10. One trunk unit was
#    taken per zone, so the unit is re-derived from the zone.
n_laur_rezoned <- sum(samples_raw$substrate == "Lauraceae" &
                        samples_raw$zone != 6, na.rm = TRUE)
samples_raw$zone[samples_raw$substrate == "Lauraceae"] <- 6

is_trunk_wood   <- samples_raw$substrate == "Ficus wood" & samples_raw$zone %in% 1:5
trunk_unit      <- paste0("S", 16 - samples_raw$zone)
n_trunk_renamed <- sum(is_trunk_wood & trimws(samples_raw$unit) != trunk_unit)
samples_raw$unit[is_trunk_wood] <- trunk_unit[is_trunk_wood]

cat("Curator corrections:", n_laur_rezoned, "Lauraceae isolates moved to zone 6;",
    n_trunk_renamed, "Ficus trunk isolates relabelled to S11-S15\n")

# --- Orientation cleaning -------------------------------------
# Field sheets carry a replicate index on some bearings (N1..N4,
# NW1..NW4): the digit identifies the branch, not the direction, and
# the field team collapsed it for most - but not all - records.
# "?" marks material whose orientation could not be reconciled.
clean_orientation <- function(x) {
  o <- toupper(trimws(as.character(x)))
  o <- sub("[0-9]+$", "", o)          # N3 -> N, NW1 -> NW
  o <- trimws(o)
  o[o %in% c("", "?", "NA", "X")] <- NA_character_
  o
}

# Compass bearing (degrees clockwise from North) for circular tests
orientation_bearing <- c(N = 0, NE = 45, E = 90, SE = 135,
                         S = 180, SW = 225, W = 270, NW = 315)

# Standardise substrate labels
samples_raw <- samples_raw %>%
  mutate(
    substrate = case_when(
      substrate == "Ficus leaves" ~ "Ficus leaves",
      substrate == "Ficus wood"   ~ "Ficus wood",
      substrate == "Lauraceae"    ~ "Lauraceae leaves",
      TRUE                        ~ substrate
    ),
    zone = as.integer(zone),
    # Two individual trees were sampled, and they are physically
    # interlocked: a strangler Ficus growing on a single Lauraceae host.
    # This is the inferential unit, so "host effect" here means
    # "these two individuals", not "these two species" (see README).
    tree_id = ifelse(substrate == "Lauraceae leaves",
                     "Lauraceae (host tree)", "Ficus (strangler)"),
    # The Lauraceae trunk was inaccessible - entirely encased by the
    # strangling Ficus - so ALL Lauraceae material came from exposed
    # branches. Stated explicitly rather than left to the Ficus rule
    # (zone <= 5 = trunk), which only happens to agree now that every
    # Lauraceae unit is recorded in zone 6.
    position = case_when(
      substrate == "Lauraceae leaves" ~ "Branch",
      zone <= 5                       ~ "Trunk",
      TRUE                            ~ "Branch"
    ),
    unit = trimws(unit),
    orientation = clean_orientation(orientation_raw),
    orientation_deg = unname(orientation_bearing[orientation]),
    # LOT2 was collected in Peru (southern hemisphere): the sun-facing
    # side of a tree is the NORTHERN one, the shaded side the southern.
    sun_aspect = case_when(
      orientation %in% c("N", "NE", "NW") ~ "Sun-facing (N)",
      orientation %in% c("S", "SE", "SW") ~ "Shaded (S)",
      orientation %in% c("E", "W")        ~ "Lateral (E/W)",
      TRUE                                ~ NA_character_
    ),
    sample_id = paste(substrate, paste0("Z", zone), unit, sep = "__"),
    # Orientation-aware sampling unit: a single unit can span more than
    # one bearing, so orientation analyses need this finer grouping.
    sample_id_orient = ifelse(
      is.na(orientation), NA_character_,
      paste(substrate, paste0("Z", zone), unit, orientation, sep = "__"))
  )

samples_raw <- samples_raw %>% mutate(culture_code1 = squish(culture_code1))

# --- Sampling depth ------------------------------------------
# Isolates recovered per sampling unit varies roughly two-fold between
# substrates, and depth is itself a predictor of community composition,
# so it is carried through as a covariate rather than ignored.
depth_by_unit <- samples_raw %>% count(sample_id, name = "unit_depth")
samples_raw <- samples_raw %>% left_join(depth_by_unit, by = "sample_id")

# --- Sampling design -----------------------------------------
# The units that were actually collected, as specified by the data
# curator. Units are compared against this list because a mislabelled
# zone or unit does not fail - it silently becomes an extra "replicate"
# (which is how the Lauraceae ended up with 10 units instead of 8).
sampling_design <- bind_rows(
  data.frame(substrate = "Ficus leaves",     zone = 6L,  unit = paste0("S", 1:10)),
  data.frame(substrate = "Ficus wood",       zone = 6L,  unit = paste0("S", 1:8)),
  data.frame(substrate = "Ficus wood",       zone = 5:1, unit = paste0("S", 11:15)),
  data.frame(substrate = "Lauraceae leaves", zone = 6L,  unit = paste0("S", 1:8))
) %>%
  mutate(sample_id = paste(substrate, paste0("Z", zone), unit, sep = "__"))

units_off_design   <- setdiff(unique(samples_raw$sample_id), sampling_design$sample_id)
units_no_isolates  <- setdiff(sampling_design$sample_id, unique(samples_raw$sample_id))
design_units <- sampling_design %>%
  group_by(substrate) %>%
  summarise(n_design  = n(),
            n_with_isolates = sum(sample_id %in% samples_raw$sample_id),
            .groups = "drop")

# --- Cross-file integrity checks -----------------------------
# The workbooks are maintained separately and have drifted before
# (a totals row read as a taxon, genotype names spelled two ways).
# Fail loudly rather than silently analysing inconsistent inputs.
check_inputs <- function() {
  for (s in names(pooled_wb)) {
    pw <- pooled_wb[[s]]
    pooled_tot <- c(
      `Lauraceae leaves` = sum(pw$n_fungi_laur_leaf, na.rm = TRUE),
      `Ficus leaves`     = sum(pw$n_fungi_fic_leaf,  na.rm = TRUE),
      `Ficus wood`       = sum(pw$n_fungi_fic_wood,  na.rm = TRUE))
    sample_tot <- table(samples_raw$substrate)[names(pooled_tot)]
    diffs <- pooled_tot - as.integer(sample_tot)
    if (any(diffs != 0)) {
      cat("WARNING:", s, "workbook and sample-level totals disagree\n")
      print(data.frame(substrate = names(pooled_tot), pooled = as.integer(pooled_tot),
                       samples = as.integer(sample_tot), difference = as.integer(diffs)))
    } else {
      cat("Integrity check:", s, "workbook and sample-level totals agree",
          paste0("(", sum(pooled_tot), " isolates)\n"))
    }
  }
  # Same culture code = subsamples of one colony, so repeated codes would
  # have to be collapsed before abundances mean anything.
  n_codes <- n_distinct(samples_raw$culture_code1)
  if (n_codes < nrow(samples_raw)) {
    cat("WARNING:", nrow(samples_raw) - n_codes, "repeated culture codes -",
        "these are subsamples of one colony and must be de-replicated.\n")
  } else {
    cat("Integrity check:", n_codes, "distinct culture codes for",
        nrow(samples_raw), "isolates - no colony subsampling to collapse\n")
  }

  if (length(units_off_design) > 0) {
    warning(length(units_off_design), " sampling unit(s) are not in the sampling design ",
            "(mislabelled zone or unit?): ", paste(units_off_design, collapse = ", "))
  } else {
    cat("Integrity check: all", n_distinct(samples_raw$sample_id),
        "sampling units match the sampling design\n")
  }
  if (length(units_no_isolates) > 0)
    cat("Note: no isolates recorded for designed unit(s):",
        paste(units_no_isolates, collapse = ", "), "\n")
  invisible(NULL)
}
check_inputs()

# =============================================================
# Taxonomic unit sets
# Every analysis is run in parallel on two definitions of the unit:
#   ITS100  ITS name from the GenBank BLAST top score at 100% similarity
#           (V. Hofstetter) - the reference set
#   ITS99   same at 99% similarity (V. Hofstetter, in production)
# =============================================================
pool_isolates <- function(s) s %>%
  group_by(its_taxon) %>%
  summarise(culture_code      = sort(culture_code1)[1],
            n_fungi_laur_leaf = sum(substrate == "Lauraceae leaves"),
            n_fungi_fic_leaf  = sum(substrate == "Ficus leaves"),
            n_fungi_fic_wood  = sum(substrate == "Ficus wood"),
            across(all_of(tax_hierarchy), majority),
            lifestyle = majority(lifestyle),
            .groups = "drop")

# Isolate names are mapped onto the curated workbook names, tolerating
# spelling variants and trailing "(...)" comments.
map_to_workbook <- function(iso_names, wb) {
  pk   <- name_key(wb$its_taxon)
  idx  <- match(name_key(iso_names), pk)
  miss <- is.na(idx)
  idx[miss] <- match(name_key_stripped(iso_names[miss]), pk)
  idx
}

build_name_set <- function(set, raw_col) {
  wb  <- pooled_wb[[set]]
  s   <- samples_raw
  raw <- squish(s[[raw_col]])
  idx <- map_to_workbook(raw, wb)
  s$its_taxon <- ifelse(is.na(idx), raw, wb$its_taxon[idx])
  for (r in tax_hierarchy) s[[r]] <- wb[[r]][idx]
  s$lifestyle <- wb$lifestyle[idx]
  s$lineage_source <- ifelse(is.na(idx), NA_character_, "workbook")

  # Unmatched names: lineage from the leading genus token
  um <- is.na(idx)
  left <- is.na(s$lineage_source)
  if (any(left)) {
    gen_lin <- wb %>% filter(!is.na(genus)) %>% group_by(genus) %>%
      summarise(across(c(phylum, subphylum, class, subclass, order, family), majority), .groups = "drop")
    tok <- sub(" .*", "", raw)
    j <- match(tok, gen_lin$genus)
    take <- left & !is.na(j)
    for (r in c("phylum", "subphylum", "class", "subclass", "order", "family"))
      s[[r]][take] <- gen_lin[[r]][j[take]]
    s$genus[take] <- tok[take]
    s$lineage_source[take] <- "genus name"
  }
  s$lineage_source[is.na(s$lineage_source)] <- "none"

  write.csv(s[um, ] %>% count(isolate_name = raw[um], unit = its_taxon, lineage_source, name = "n_isolates"),
            file.path("results/data_checks", paste0(set, "_names_not_in_workbook.csv")), row.names = FALSE)

  # Per-taxon counts: workbook vs isolate records
  cmp <- full_join(
    wb %>% group_by(its_taxon) %>%
      summarise(wb_lauraceae = sum(n_fungi_laur_leaf), wb_ficus_leaves = sum(n_fungi_fic_leaf),
                wb_ficus_wood = sum(n_fungi_fic_wood), .groups = "drop"),
    pool_isolates(s) %>% select(its_taxon, iso_lauraceae = n_fungi_laur_leaf,
                                iso_ficus_leaves = n_fungi_fic_leaf, iso_ficus_wood = n_fungi_fic_wood),
    by = "its_taxon") %>%
    mutate(across(-its_taxon, ~ replace_na(., 0))) %>%
    filter(wb_lauraceae != iso_lauraceae | wb_ficus_leaves != iso_ficus_leaves |
             wb_ficus_wood != iso_ficus_wood)
  write.csv(cmp, file.path("results/data_checks", paste0(set, "_workbook_vs_isolate_counts.csv")),
            row.names = FALSE)
  cat(set, ":", n_distinct(s$its_taxon), "taxa;", sum(um), "isolates with a name not in the workbook;",
      nrow(cmp), "taxa whose workbook counts differ from the isolate records\n")
  s
}

unit_sets <- list(
  ITS100 = list(label = "ITS 100% BLAST name",
                samples = build_name_set("ITS100", "its100_raw")),
  ITS99  = list(label = "ITS 99% BLAST name (in production)",
                samples = build_name_set("ITS99", "its99_raw")))

for (u in names(unit_sets)) {
  s <- unit_sets[[u]]$samples
  # A Lifestyle column in the samples workbook overrides the per-taxon one
  s$lifestyle <- dplyr::coalesce(s$lifestyle_isolate, s$lifestyle)
  unit_sets[[u]]$samples <- s
  unit_sets[[u]]$pooled  <- pool_isolates(s)
}
have_lifestyle <- any(!is.na(unit_sets$ITS100$samples$lifestyle))
cat("Lifestyle column:", if (have_lifestyle) "found" else
      "not found in any workbook - lifestyle analyses will be skipped", "\n")

cat("Sample records:", nrow(samples_raw), "isolates\n")
cat("Substrates:", paste(sort(unique(samples_raw$substrate)), collapse = ", "), "\n")
cat("Zones:", paste(sort(unique(samples_raw$zone)), collapse = ", "), "\n")
cat("Orientations:", paste(sort(unique(na.omit(samples_raw$orientation))), collapse = ", "),
    "(", sum(is.na(samples_raw$orientation)), "isolates unresolved )\n")
