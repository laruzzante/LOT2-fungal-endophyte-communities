# =============================================================
# Sampling completeness and replicated alpha diversity
#
# The pooled-count diversity table (community_analysis.R, section A)
# gives ONE value per substrate, so nothing there can be tested. The
# sample-level file has 10-13 independent sampling units per substrate,
# which supports both a proper test of alpha diversity and an estimate
# of how much of the community was actually recovered.
#
# Reported here:
#   1. Richness estimators (Chao1 / ACE) and Good's coverage, which say
#      how far the sampling is from complete.
#   2. Alpha diversity per SAMPLING UNIT, with a Kruskal-Wallis test
#      across substrates and pairwise Wilcoxon follow-ups.
#   3. Hill numbers (q = 0, 1, 2), the modern replacement for reporting
#      richness / Shannon / Simpson on different scales.
# =============================================================
library(vegan)

dir.create("plots/diversity", recursive = TRUE, showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)

cat("\n====== Sampling completeness & replicated alpha diversity ======\n")
cat("Analytical unit:", unit_label, "\n")

div_unit_col <- primary_unit

# ------------------------------------------------------------
# 1. How complete is the sampling?
# ------------------------------------------------------------
comm_by_substrate <- samples_raw %>%
  filter(!is.na(.data[[div_unit_col]])) %>%
  count(substrate, .data[[div_unit_col]]) %>%
  pivot_wider(names_from = all_of(div_unit_col), values_from = n, values_fill = 0) %>%
  as.data.frame()
rownames(comm_by_substrate) <- comm_by_substrate$substrate
mat_sub <- as.matrix(comm_by_substrate[, -1, drop = FALSE])

est <- as.data.frame(t(estimateR(mat_sub)))
# Good's coverage: the share of individuals belonging to taxa seen more
# than once - i.e. the probability that the next isolate collected would
# already be represented in the sample.
goods <- apply(mat_sub, 1, function(x) 1 - sum(x == 1) / sum(x))

completeness <- data.frame(
  substrate       = rownames(mat_sub),
  n_units         = as.integer(table(samples_raw$substrate)[rownames(mat_sub)] > 0) * NA,
  n_isolates      = as.integer(rowSums(mat_sub)),
  observed_S      = as.integer(est$S.obs),
  chao1           = round(est$S.chao1, 1),
  chao1_se        = round(est$se.chao1, 1),
  ace             = round(est$S.ACE, 1),
  pct_of_chao1    = round(100 * est$S.obs / est$S.chao1),
  singletons      = as.integer(apply(mat_sub, 1, function(x) sum(x == 1))),
  doubletons      = as.integer(apply(mat_sub, 1, function(x) sum(x == 2))),
  goods_coverage  = round(100 * goods, 1),
  row.names = NULL
)
completeness$n_units <- as.integer(
  samples_raw %>% distinct(substrate, sample_id) %>% count(substrate) %>%
    arrange(match(substrate, completeness$substrate)) %>% pull(n))

write.csv(completeness, "tables/sampling_completeness.csv", row.names = FALSE)
cat("\nSampling completeness per substrate:\n")
print(completeness)
cat("\nObserved richness recovers ",
    min(completeness$pct_of_chao1), "-", max(completeness$pct_of_chao1),
    "% of the Chao1 estimate; Good's coverage ",
    min(completeness$goods_coverage), "-", max(completeness$goods_coverage),
    "%.\n", sep = "")

p_comp <- completeness %>%
  select(substrate, Observed = observed_S, Chao1 = chao1, ACE = ace) %>%
  pivot_longer(-substrate, names_to = "estimator", values_to = "richness") %>%
  mutate(estimator = factor(estimator, levels = c("Observed", "Chao1", "ACE")))
se_lookup <- setNames(completeness$chao1_se, completeness$substrate)
p1 <- ggplot(p_comp, aes(x = substrate, y = richness, fill = estimator)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    data = transform(subset(p_comp, estimator == "Chao1"),
                     se = se_lookup[substrate]),
    aes(ymin = richness - se, ymax = richness + se),
    width = 0.15, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c(Observed = "#2E86AB", Chao1 = "#F18F01", ACE = "#A23B72")) +
  labs(title = "Observed vs estimated richness",
       subtitle = paste0("Unit: ", unit_label,
                         " | error bar = 1 SE of Chao1 | less than ",
                         max(completeness$pct_of_chao1),
                         "% of estimated richness was recovered"),
       x = NULL, y = "Number of taxa", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")
ggsave("plots/diversity/richness_estimators.pdf", p1, width = 9, height = 6)

# ------------------------------------------------------------
# 2. Alpha diversity per sampling unit (this is what can be tested)
# ------------------------------------------------------------
comm_unit <- samples_raw %>%
  filter(!is.na(.data[[div_unit_col]])) %>%
  count(sample_id, .data[[div_unit_col]]) %>%
  pivot_wider(names_from = all_of(div_unit_col), values_from = n, values_fill = 0) %>%
  as.data.frame()
rownames(comm_unit) <- comm_unit$sample_id
mat_unit <- as.matrix(comm_unit[, -1, drop = FALSE])

unit_meta <- samples_raw %>%
  distinct(sample_id, substrate, zone, position) %>%
  as.data.frame()
unit_meta <- unit_meta[match(rownames(mat_unit), unit_meta$sample_id), ]

H <- diversity(mat_unit, "shannon")
S <- specnumber(mat_unit)
N <- rowSums(mat_unit)

alpha_units <- data.frame(
  sample_id   = rownames(mat_unit),
  substrate   = unit_meta$substrate,
  zone        = unit_meta$zone,
  position    = unit_meta$position,
  n_isolates  = as.integer(N),
  richness_S  = as.integer(S),
  shannon_H   = round(H, 4),
  simpson_1mD = round(diversity(mat_unit, "simpson"), 4),
  pielou_J    = round(ifelse(S > 1, H / log(S), NA_real_), 4),
  # Hill numbers: all three on the same scale (effective taxon counts)
  hill_q0     = as.integer(S),
  hill_q1     = round(exp(H), 3),
  hill_q2     = round(diversity(mat_unit, "invsimpson"), 3),
  row.names   = NULL
)
write.csv(alpha_units, "tables/alpha_diversity_per_unit.csv", row.names = FALSE)

alpha_summary <- alpha_units %>%
  group_by(substrate) %>%
  summarise(n_units = n(),
            isolates = sum(n_isolates),
            mean_S = round(mean(richness_S), 2), sd_S = round(sd(richness_S), 2),
            mean_H = round(mean(shannon_H), 3),  sd_H = round(sd(shannon_H), 3),
            mean_q1 = round(mean(hill_q1), 2),
            mean_J = round(mean(pielou_J, na.rm = TRUE), 3),
            .groups = "drop") %>%
  as.data.frame()
write.csv(alpha_summary, "tables/alpha_diversity_per_unit_summary.csv", row.names = FALSE)
cat("\nAlpha diversity per sampling unit (mean +/- SD across units):\n")
print(alpha_summary)

# ------------------------------------------------------------
# 3. Test it. Non-parametric: 10-13 units per group, and the indices
#    are bounded and skewed, so Kruskal-Wallis rather than ANOVA.
# ------------------------------------------------------------
alpha_tests <- list()
for (idx in c("richness_S", "shannon_H", "hill_q1", "pielou_J")) {
  v <- alpha_units[[idx]]
  keep <- !is.na(v)
  kw <- kruskal.test(v[keep] ~ factor(alpha_units$substrate[keep]))
  alpha_tests[[idx]] <- data.frame(
    index = idx, test = "Kruskal-Wallis",
    statistic = round(unname(kw$statistic), 4),
    df = unname(kw$parameter),
    p_value = signif(kw$p.value, 4),
    stringsAsFactors = FALSE)
}
alpha_test_table <- bind_rows(alpha_tests)
# Four indices on the same data are not independent hypotheses, but a
# holm adjustment keeps the family-wise claim honest.
alpha_test_table$p_adj <- signif(p.adjust(alpha_test_table$p_value, "holm"), 4)
write.csv(alpha_test_table, "tables/alpha_diversity_tests.csv", row.names = FALSE)
cat("\nAlpha diversity across substrates:\n")
print(alpha_test_table)

pw_list <- list()
for (idx in c("richness_S", "shannon_H")) {
  v <- alpha_units[[idx]]
  pw <- pairwise.wilcox.test(v, alpha_units$substrate,
                             p.adjust.method = "holm", exact = FALSE)
  m <- as.data.frame(as.table(pw$p.value))
  names(m) <- c("group1", "group2", "p_adj")
  m <- m[!is.na(m$p_adj), , drop = FALSE]
  m$index <- idx
  m$p_adj <- signif(m$p_adj, 4)
  pw_list[[idx]] <- m[, c("index", "group1", "group2", "p_adj")]
}
alpha_pairwise <- bind_rows(pw_list)
write.csv(alpha_pairwise, "tables/alpha_diversity_pairwise.csv", row.names = FALSE)
cat("\nPairwise (Wilcoxon, Holm-adjusted):\n")
print(alpha_pairwise)

substrate_colours_div <- c("Ficus leaves" = "#A23B72", "Ficus wood" = "#F18F01",
                           "Lauraceae leaves" = "#2E86AB")
for (idx in c("richness_S", "shannon_H", "hill_q1")) {
  pv <- alpha_test_table$p_adj[alpha_test_table$index == idx]
  p <- ggplot(alpha_units, aes(x = substrate, y = .data[[idx]], fill = substrate)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.55) +
    geom_jitter(width = 0.12, height = 0, size = 2.2, alpha = 0.9) +
    scale_fill_manual(values = substrate_colours_div, guide = "none") +
    labs(title = paste("Per-sampling-unit", idx),
         subtitle = paste0("Each point is one sampling unit | Kruskal-Wallis Holm-adjusted p = ",
                           pv, " | unit: ", unit_label),
         x = NULL, y = idx) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path("plots/diversity", paste0("alpha_per_unit_", idx, ".pdf")),
         p, width = 8, height = 6)
}

# ------------------------------------------------------------
# 4. Coverage-standardised comparison.
#    Comparing richness at equal COUNT still favours whichever group
#    happens to have the flatter abundance distribution; comparing at
#    equal COVERAGE is the fairer contrast (Chao & Jost 2012).
# ------------------------------------------------------------
target_cov <- floor(min(completeness$goods_coverage))
cat("\nRarefying each substrate to a common Good's coverage of ~",
    target_cov, "%\n", sep = "")

# Coverage of a size-n subsample has no simple closed form, so it is
# estimated by repeated subsampling: draw n isolates, count how many
# fall in taxa seen only once in that draw.
cov_at_n <- function(counts, n, reps = 200) {
  if (n >= sum(counts)) return(1 - sum(counts == 1) / sum(counts))
  m <- matrix(counts, nrow = 1)
  mean(replicate(reps, {
    sub <- as.numeric(rrarefy(m, n))
    1 - sum(sub == 1) / sum(sub)
  }))
}
set.seed(42)
cov_rows <- list()
for (sb in rownames(mat_sub)) {
  cnt <- mat_sub[sb, ]; cnt <- cnt[cnt > 0]
  tot <- sum(cnt)
  # coverage rises monotonically with n, so scan a grid and take the
  # smallest n that reaches the target
  ns <- unique(round(seq(10, tot, length.out = 40)))
  covs <- vapply(ns, function(n) cov_at_n(cnt, n), 1)
  hit <- which(covs >= target_cov / 100)
  n_at <- if (length(hit)) ns[min(hit)] else tot
  cov_rows[[sb]] <- data.frame(
    substrate = sb, isolates_needed = n_at,
    coverage_reached = round(100 * cov_at_n(cnt, n_at), 1),
    richness_at_coverage = round(as.numeric(rarefy(matrix(cnt, nrow = 1), n_at)), 1))
}
coverage_std <- bind_rows(cov_rows)
write.csv(coverage_std, "tables/richness_at_equal_coverage.csv", row.names = FALSE)
cat("\nRichness standardised to equal coverage:\n")
print(coverage_std)

cat("\n====== Diversity analyses complete ======\n")
