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
    position = ifelse(zone <= 5, "Trunk", "Branch"),
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

cat("Pooled counts:", nrow(pooled), "genotypes\n")
cat("Sample records:", nrow(samples_raw), "isolates\n")
cat("Substrates:", paste(sort(unique(samples_raw$substrate)), collapse = ", "), "\n")
cat("Zones:", paste(sort(unique(samples_raw$zone)), collapse = ", "), "\n")
cat("Orientations:", paste(sort(unique(na.omit(samples_raw$orientation))), collapse = ", "),
    "(", sum(is.na(samples_raw$orientation)), "isolates unresolved )\n")
