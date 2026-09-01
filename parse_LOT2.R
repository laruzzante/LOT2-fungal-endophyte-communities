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
samples_raw <- read_excel("LOT2_samples.xlsx", sheet = 1)
names(samples_raw) <- c("substrate", "zone", "unit", "orientation",
                         "culture_code1", "its_taxon", "genbank",
                         "culture_code2")

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
    sample_id = paste(substrate, paste0("Z", zone), unit, sep = "__")
  )

cat("Pooled counts:", nrow(pooled), "genotypes\n")
cat("Sample records:", nrow(samples_raw), "isolates\n")
cat("Substrates:", paste(sort(unique(samples_raw$substrate)), collapse = ", "), "\n")
cat("Zones:", paste(sort(unique(samples_raw$zone)), collapse = ", "), "\n")


