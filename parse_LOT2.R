# =============================================================
# Parse LOT2 input data
# =============================================================

# --- Pooled genotype counts (first sheet only) ---
pooled <- read_excel("LOT2_pooled_counts.xlsx", sheet = 1) %>%
  rename(
    culture_code      = `Hofstetter-culture code (CTAB 1x)`,
    its_taxon         = `ITS sequences named after GenBank BLAST top score(s) result(s) and taxon current names (Mycobank/Index Fungorum)`,
    n_fungi_laur_leaf = `Number of fungi isolated from the Lauraceae leaves`,
    n_fungi_fic_leaf  = `Number of fungi isolated from the Ficus leaves`,
    n_fungi_fic_wood  = `Number of fungi isolated from from theFicus wood`,
    phylum            = Phylum,
    subphylum         = Subphylum,
    superclass        = `Superclass (Fungi are sometimes reported under both Leotiomycetes and Sordariomycetes due to historical classification errors, overlapping physical traits, and updates from modern DNA sequencing.`,
    class             = Class,
    subclass          = Subclass,
    order             = Order,
    family            = Family,
    genus             = Genus,
    species           = Species
  )

# The workbook ends with a spreadsheet TOTALS row: no culture code, no
# taxon name, just the column sums. Read as a genotype it becomes a
# phantom taxon carrying an entire substrate's isolate count, which the
# pie and abundance charts then draw as a giant "Incertae sedis" slice.
# Drop every row that carries neither an identifier nor a name.
pooled_raw_rows <- nrow(pooled)
pooled <- pooled %>% filter(!(is.na(culture_code) & is.na(its_taxon)))
if (pooled_raw_rows > nrow(pooled))
  cat("Dropped", pooled_raw_rows - nrow(pooled),
      "unidentified/total row(s) from the pooled counts\n")

# --- Sample-level isolate records (first sheet only) ---
# The samples workbook was re-issued by the field team with the
# reconciled branch numbering / colour codes, which finally makes the
# branch & trunk sampling ORIENTATION usable. Layout (12 columns):
#   1 substrate            5 branch colour (project database)   9  ITS taxon
#   2 tree zone            6 orientation (pre-reconciliation)  10  GenBank hit
#   3 sampling unit        7 orientation (reconciled - USE)    11  culture code (dup)
#   4 branch colour (field) 8 culture code                     12  field notes
# Columns 4-6 and 12 are the reconciliation working columns kept for
# traceability; only column 7 is authoritative for orientation.
samples_raw <- read_excel("LOT2_samples.xlsx", sheet = 1)
names(samples_raw) <- c("substrate", "zone", "unit",
                        "branch_colour_field", "branch_colour_db",
                        "orientation_old", "orientation_raw",
                        "culture_code1", "its_taxon", "genbank",
                        "culture_code2", "field_note")

# Drop any fully empty trailing rows
samples_raw <- samples_raw[!is.na(samples_raw$substrate), ]

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
    # Zone-to-position is NOT the same rule for both trees. The
    # Lauraceae trunk was inaccessible - entirely encased by the
    # strangling Ficus - so ALL Lauraceae material came from exposed
    # upper branches, including the zone-5 units. Applying the Ficus
    # rule (zone <= 5 = trunk) to it mislabels 3 units / 53 isolates.
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

# --- Sequence-based OTUs -------------------------------------
# `its_taxon` is a BLAST-derived name string, not a clustered unit: it
# over-splits (near-identical sequences under two spellings) and
# over-lumps (single labels covering "26 spp."). `cluster_otus.sh`
# clusters the ITS sequences themselves; the mapping is committed so the
# pipeline runs without vsearch.
otu_map_file <- "LOT2_otu_map.csv"
if (file.exists(otu_map_file)) {
  otu_map <- read.csv(otu_map_file, colClasses = "character")
  samples_raw <- samples_raw %>%
    mutate(culture_code1 = trimws(as.character(culture_code1))) %>%
    left_join(otu_map, by = c("culture_code1" = "culture_code"))
  n_missing_otu <- sum(is.na(samples_raw$otu_97))
  if (n_missing_otu > 0)
    warning(n_missing_otu, " isolates have no OTU assignment; ",
            "re-run cluster_otus.sh after changing the input sequences.")
  cat("OTUs (97% ITS identity):", n_distinct(na.omit(samples_raw$otu_97)),
      "from", nrow(samples_raw), "isolates",
      sprintf("(%d name-based labels)\n", n_distinct(samples_raw$its_taxon)))
  have_otus <- TRUE
} else {
  samples_raw$otu_97 <- NA_character_
  samples_raw$otu_985 <- NA_character_
  have_otus <- FALSE
  cat("NOTE:", otu_map_file, "not found - falling back to name-based genotypes.\n",
      "     Run `bash cluster_otus.sh` to enable sequence-based OTUs.\n")
}

# --- Sampling depth ------------------------------------------
# Isolates recovered per sampling unit varies roughly two-fold between
# substrates, and depth is itself a predictor of community composition,
# so it is carried through as a covariate rather than ignored.
depth_by_unit <- samples_raw %>% count(sample_id, name = "unit_depth")
samples_raw <- samples_raw %>% left_join(depth_by_unit, by = "sample_id")

# --- Cross-file integrity checks -----------------------------
# The two workbooks are maintained separately and have drifted before
# (a totals row read as a taxon, genotype names spelled two ways).
# Fail loudly rather than silently analysing inconsistent inputs.
check_inputs <- function() {
  pooled_tot <- c(
    `Lauraceae leaves` = sum(pooled$n_fungi_laur_leaf, na.rm = TRUE),
    `Ficus leaves`     = sum(pooled$n_fungi_fic_leaf,  na.rm = TRUE),
    `Ficus wood`       = sum(pooled$n_fungi_fic_wood,  na.rm = TRUE))
  sample_tot <- table(samples_raw$substrate)[names(pooled_tot)]
  diffs <- pooled_tot - as.integer(sample_tot)
  if (any(diffs != 0)) {
    cat("WARNING: pooled and sample-level totals disagree\n")
    print(data.frame(substrate = names(pooled_tot),
                     pooled = as.integer(pooled_tot),
                     samples = as.integer(sample_tot),
                     difference = as.integer(diffs)))
  } else {
    cat("Integrity check: pooled and sample-level totals agree",
        paste0("(", sum(pooled_tot), " isolates)\n"))
  }
  # Same culture code = subsamples of one colony, so repeated codes would
  # have to be collapsed before abundances mean anything. They are all
  # distinct here, so every row is an independent isolate.
  n_codes <- n_distinct(samples_raw$culture_code1)
  if (n_codes < nrow(samples_raw)) {
    cat("WARNING:", nrow(samples_raw) - n_codes, "repeated culture codes -",
        "these are subsamples of one colony and must be de-replicated.\n")
  } else {
    cat("Integrity check:", n_codes, "distinct culture codes for",
        nrow(samples_raw), "isolates - no colony subsampling to collapse\n")
  }

  unmatched <- setdiff(unique(samples_raw$its_taxon), unique(pooled$its_taxon))
  if (length(unmatched) > 0)
    cat("Note:", length(unmatched), "sample genotype names are spelled",
        "differently in the pooled file\n     (lineage is resolved per OTU,",
        "so these are rescued by their OTU-mates).\n")
  invisible(NULL)
}
check_inputs()

cat("Pooled counts:", nrow(pooled), "genotypes\n")
cat("Sample records:", nrow(samples_raw), "isolates\n")
cat("Substrates:", paste(sort(unique(samples_raw$substrate)), collapse = ", "), "\n")
cat("Zones:", paste(sort(unique(samples_raw$zone)), collapse = ", "), "\n")
cat("Orientations:", paste(sort(unique(na.omit(samples_raw$orientation))), collapse = ", "),
    "(", sum(is.na(samples_raw$orientation)), "isolates unresolved )\n")
