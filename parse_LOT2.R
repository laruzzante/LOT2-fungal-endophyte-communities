xlsx <- "LOT2_for_Livio.xlsx"

# --- Sheet: ALL isolates ---
all_isolates <- read_excel(xlsx, sheet = "ALL isolates") %>%
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

# --- Sheet: 67. Fungi-Endo leaf (Ficus) ---
endo_leaf_ficus <- read_excel(xlsx, sheet = "67. Fungi-Endo leaf (Ficus)") %>%
  rename(
    sample_code      = `LOT-SampleCode`,
    culture_code     = `Hofstetter-culture code (CTAB 1x)`,
    tree             = Tree,
    tree_orient      = `Tree orientation`,
    its_taxon        = `ITS sequences named after GenBank BLAST top score(s) result(s) and taxon current names (Mycobank/Index Fungorum)`,
    n_species        = `Number of species`,
    n_genotypes      = `Number of genotypes`,
    genbank_blast_acc = `GenBank BLAST top score sequence(s) accession(s); sequence similarity/sequence query coverage both expressed in %`,
    seq_code         = `Hofstetter-sequence code`,
    its_seq          = `LOT2  ITS sequences`
  )

# --- Sheet: 65. Fungi-Endo wood (Ficus) ---
endo_wood_ficus <- read_excel(xlsx, sheet = "65. Fungi-Endo wood (Ficus)") %>%
  rename(
    sample_code      = `LOT-SampleCode`,
    culture_code     = `Hofstetter-culture code (CTAB 1x)`,
    tree_zone        = `Tree zone`,
    tree_branch      = `Tree branch/trunk`,
    its_taxon        = `ITS sequences named after GenBank BLAST top score(s) result(s) and taxon current names (Mycobank/Index Fungorum)`,
    n_species        = `Number of species`,
    n_genotypes      = `Number of genotypes`,
    genbank_blast_acc = `GenBank BLAST top score sequence(s) accession(s); sequence similarity/sequence query coverage both expressed in %`,
    its_seq          = `ITS sequences`
  )

# --- Sheet: 66. Fungi-Endo wood (Host) ---
endo_wood_host <- read_excel(xlsx, sheet = "66. Fungi-Endo wood (Host)") %>%
  rename(
    sample_code      = `LOT-SampleCode`,
    culture_code     = `Hofstetter-Code`,
    tree_zone        = `Tree Zone`,
    tree_branch      = `Tree branch`,
    its_taxon        = `ITS sequences named after GenBank BLAST top score(s) result(s) and taxon current names (Mycobank/Index Fungorum)`,
    n_species        = `Number of species`,
    n_genotypes      = `Number of genotypes`,
    genbank_acc      = `GenBank accessions for ITS sequences`,
    genbank_blast_acc = `GenBank BLAST top score sequence(s) accession(s); sequence similarity/sequence query coverage both expressed in %`,
    classif_genbank  = `Classification after GenBank`,
    seq_code         = `Hofstetter-sequence code`,
    its_seq          = `ITS sequences`
  )

# --- Sheet: Feuil2 (sequence data) ---
seq_data <- read_excel(xlsx, sheet = "Feuil2") %>%
  setNames(c("lot_prefix", "culture_code", "its_taxon", "its_seq"))


