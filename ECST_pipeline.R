#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--abundance", type = "character"),
  make_option("--metadata",  type = "character"),
  make_option("--outdir",    type = "character", default = "output"),
  make_option("--exclude",   type = "character", default = ""),
  make_option("--mean_abd_thr", type = "double", default = 0.0054),
  make_option("--min_cohorts",  type = "integer", default = 2L),
  make_option("--prev",      type = "double",    default = 0),
  make_option("--manual_k",  type = "integer",   default = 4L),
  make_option("--kmax",      type = "integer",   default = 8L),
  make_option("--seed",      type = "integer",   default = 2025L)
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$abundance) || is.null(opt$metadata))
  stop("Must specify --abundance and --metadata")

ABUNDANCE_FILE   <- opt$abundance
METADATA_FILE    <- opt$metadata
BASE_OUT         <- opt$outdir
EXCLUDED_COHORTS <- if (nchar(opt$exclude) > 0) strsplit(opt$exclude, ",")[[1]] else character(0)
MEAN_ABD_THR     <- opt$mean_abd_thr
MIN_COHORTS      <- opt$min_cohorts
PREV_THR         <- opt$prev
MANUAL_K         <- if (opt$manual_k == 0L) NULL else opt$manual_k
K_MAX            <- opt$kmax
SEED             <- opt$seed
ABD_THR          <- 1e-5
PSEUDO           <- 1e-5

if (!is.null(MANUAL_K)) {
  if (MANUAL_K < 2 || MANUAL_K > K_MAX)
    stop(sprintf("--manual_k must be between 2 and %d, got %d", K_MAX, MANUAL_K))
}

cat("════════════════════════════════════════════\n")
cat("ECST Pipeline (Kraken2 species-level)\n")
cat("  abundance:", ABUNDANCE_FILE, "\n")
cat("  metadata: ", METADATA_FILE, "\n")
cat("  outdir:   ", BASE_OUT, "\n")
cat("  mean_abd_thr:", MEAN_ABD_THR*100, "% (per-cohort mean)\n")
cat("  min_cohorts: ", MIN_COHORTS, "\n")
if (!is.null(MANUAL_K)) cat("  k:        ", MANUAL_K, "\n") else cat("  k:        auto\n")
cat("  seed:     ", SEED, "\n")
cat("════════════════════════════════════════════\n\n")

suppressPackageStartupMessages({
  library(MMUPHin); library(DirichletMultinomial)
  library(vegan); library(cluster); library(ape)
  library(tidyverse); library(patchwork); library(ggpubr); library(ggrepel)
  library(ComplexHeatmap); library(circlize); library(gridExtra)
})

if ("package:plyr" %in% search()) {
  detach("package:plyr", unload = FALSE, force = TRUE)
}
library(dplyr, warn.conflicts = FALSE)

make_outdir <- function(m) { d <- file.path(BASE_OUT, m); dir.create(d, recursive = TRUE, showWarnings = FALSE); d }

CNS_THEME <- theme_classic(base_size = 7) +
  theme(axis.line = element_line(linewidth = 0.3, color = "black"),
        axis.ticks = element_line(linewidth = 0.2, color = "black"),
        axis.text = element_text(size = 6, color = "black", family = "Times"),
        axis.title = element_text(size = 7, color = "black", family = "Times"),
        legend.text = element_text(size = 5.5, family = "Times"),
        legend.title = element_text(size = 6.5, face = "bold", family = "Times"),
        legend.key.size = unit(2.5, "mm"),
        legend.background = element_rect(fill = "white", color = NA),
        strip.text = element_text(size = 6.5, face = "bold", family = "Times"),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        plot.tag = element_text(size = 8, face = "bold", family = "Times"),
        plot.title = element_text(size = 8, face = "bold", hjust = 0, family = "Times"),
        plot.subtitle = element_text(size = 6, face = "italic", family = "Times"),
        plot.margin = margin(3, 3, 3, 3, "pt"),
        text = element_text(family = "Times"))

AXIS_ROTATE_45 <- theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 5))

COHORT_COLORS <- c("PRJCA002217"="#3C5488","PRJCA016672"="#E64B35","PRJEB38989"="#4DBBD5",
                   "PRJEB55149"="#00A087","PRJNA1261000"="#F39B7F","PRJNA774492"="#8491B4",
                   "PRJNA886972"="#91D1C2")

SPECIES_PALETTE <- c(
  "#D94F4F","#3674B5","#1A7A5E","#E08A28","#8555B3",
  "#2BB57A","#C74E8A","#5DA4D9","#A67B2E","#6B8E23",
  "#E06060","#4B7FA8","#9C5E2E","#7CB87C","#B074C8",
  "#D48C6A","#408080","#C25D5D","#6A6AAD","#88B04B"
)

assign_sp_colors <- function(sp_ranked) {
  n <- length(sp_ranked)
  pal <- rep(SPECIES_PALETTE, length.out = n)
  cols <- setNames(pal[1:n], sp_ranked)
  cols["Others"] <- "#E5E5E5"
  cols["Unclassified"] <- "#CCCCCC"
  cols
}

ECST_PAL <- c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#91D1C2","#7E6148")

to_proportion <- function(mat) { mat[mat < 0] <- 0; cs <- colSums(mat); cs[cs == 0] <- 1; sweep(mat, 2, cs, "/") }
clr_transform <- function(mat) { mat[mat < 0] <- 0; m <- mat + PSEUDO; apply(m, 2, function(x) log(x) - mean(log(x))) }

save_fig <- function(p, path, w = 180, h = 120) {
  ggsave(paste0(path, ".pdf"), p, width = w, height = h, units = "mm", dpi = 300, device = cairo_pdf)
  ggsave(paste0(path, ".jpg"), p, width = w, height = h, units = "mm", dpi = 300)
  ggsave(paste0(path, ".svg"), p, width = w, height = h, units = "mm", dpi = 300)
  cat("  ✓ fig:", basename(path), "\n")
}
save_tsv <- function(df, path) {
  write.table(df, paste0(path, ".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
}

# ====================================================================
# M01 — Cross-cohort consistency filter + MMUPHin batch correction
# ====================================================================
cat("\n╔══ M01: Cross-cohort filter + MMUPHin ══╗\n")
OUT_M01 <- make_outdir("M01_MMUPHin")

species_raw <- read.table(ABUNDANCE_FILE, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
meta_raw    <- read.table(METADATA_FILE, sep = "\t", header = TRUE, check.names = FALSE)
if (length(EXCLUDED_COHORTS) > 0) {
  meta_raw <- meta_raw[!meta_raw$Dataset %in% EXCLUDED_COHORTS, ]
}

original_names <- rownames(species_raw)
parsed_names <- sapply(original_names, function(nm) {
  m <- regmatches(nm, regexpr("s__[^|]+", nm))
  if (length(m) > 0) sub("^s__", "", m) else nm
})
if (any(duplicated(parsed_names))) {
  species_raw <- as.data.frame(species_raw)
  species_raw$sp_parsed <- parsed_names
  species_raw <- species_raw %>%
    group_by(sp_parsed) %>%
    summarise(across(everything(), sum)) %>%
    column_to_rownames("sp_parsed")
} else {
  rownames(species_raw) <- parsed_names
}

common  <- intersect(colnames(species_raw), meta_raw$Sample)
species <- species_raw[, common]
meta    <- meta_raw[match(common, meta_raw$Sample), ]
rownames(meta) <- meta$Sample
meta$Cohort <- meta$Dataset
cat("  Samples:", length(common), "| Species:", nrow(species), "\n")

species_prop <- to_proportion(species)
cohort_list <- unique(meta$Cohort)
cohort_mean_mat <- matrix(0, nrow=nrow(species_prop), ncol=length(cohort_list),
                          dimnames=list(rownames(species_prop), cohort_list))
for (coh in cohort_list) {
  coh_samples <- meta$Sample[meta$Cohort == coh]
  coh_samples <- intersect(coh_samples, colnames(species_prop))
  if (length(coh_samples) > 0) {
    cohort_mean_mat[, coh] <- rowMeans(species_prop[, coh_samples, drop=FALSE])
  }
}

cohort_pass <- apply(cohort_mean_mat, 1, function(x) sum(x >= MEAN_ABD_THR))
species_keep <- names(cohort_pass[cohort_pass >= MIN_COHORTS])

filter_summary <- data.frame(
  species = rownames(cohort_mean_mat),
  n_cohorts_pass = cohort_pass,
  stringsAsFactors = FALSE
)
for (coh in cohort_list) {
  filter_summary[[paste0("mean_", coh)]] <- cohort_mean_mat[, coh]
}
filter_summary <- filter_summary[order(-filter_summary$n_cohorts_pass), ]

if (length(species_keep) < 3) {
  prev <- apply(species_prop, 1, function(x) mean(x > ABD_THR))
  species_filt <- to_proportion(species_prop[prev >= PREV_THR, , drop = FALSE])
} else {
  species_filt <- to_proportion(species_prop[species_keep, , drop = FALSE])
  cat("  ✓ Cross-cohort filter retained:", length(species_keep), "species\n")
}

save_tsv(filter_summary, file.path(OUT_M01, "species_cohort_filter_summary"))

gpc <- tapply(meta$Group, meta$Cohort, function(x) length(unique(x)))
if (any(gpc <= 1)) {
  fit <- adjust_batch(feature_abd = species_filt, batch = "Cohort",
                      data = meta[colnames(species_filt),], control = list(verbose = FALSE))
} else {
  fit <- adjust_batch(feature_abd = species_filt, batch = "Cohort", covariates = "Group",
                      data = meta[colnames(species_filt),], control = list(verbose = FALSE))
}
corrected <- fit$feature_abd_adj; corrected[corrected < 0] <- 0
corrected <- to_proportion(corrected)

write.table(corrected, file.path(OUT_M01, "abundance_corrected.tsv"), sep="\t", quote=FALSE)

# ====================================================================
# M02 — Cross-cohort filter visualizations (Figure 2)
# ====================================================================
cat("\n╔══ M02: Filter visualizations ══╗\n")
OUT_M02 <- make_outdir("M02_Filter_Plots")

cohort_means_long <- as.data.frame(cohort_mean_mat) %>%
  rownames_to_column("species") %>%
  pivot_longer(-species, names_to = "Cohort", values_to = "mean_abd") %>%
  filter(mean_abd > 0)

cohort_means_long$log10_abd <- log10(cohort_means_long$mean_abd * 100 + 1e-6)

n_above_1 <- cohort_means_long %>% group_by(Cohort) %>%
  summarise(n = sum(mean_abd >= 0.01), .groups = "drop")

p_kde <- ggplot(cohort_means_long, aes(x = log10_abd, fill = Cohort)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, alpha = 0.4, color = NA) +
  geom_density(aes(color = Cohort), linewidth = 0.4, fill = NA) +
  geom_vline(xintercept = log10(1), linetype = "dashed", linewidth = 0.3, color = "grey40") +
  geom_text(data = n_above_1, aes(x = Inf, y = Inf, label = paste0("n≥1%: ", n)),
            hjust = 1.1, vjust = 1.5, size = 2, inherit.aes = FALSE) +
  facet_wrap(~ Cohort, scales = "free_y") +
  scale_color_manual(values = COHORT_COLORS, guide = "none") +
  scale_fill_manual(values = COHORT_COLORS, guide = "none") +
  labs(x = "Cohort-level mean abundance (log10 %)", y = "Density",
       title = "Per-cohort species abundance distribution") +
  CNS_THEME

save_fig(p_kde, file.path(OUT_M02, "Fig2A_KDE_perCohort"), w = 180, h = 110)

threshold_seq <- seq(0.0001, 0.05, by = 0.0001)
retention_data <- list()
for (k_min in c(2, 3, 4)) {
  retention <- sapply(threshold_seq, function(t) {
    pass <- apply(cohort_mean_mat, 1, function(x) sum(x >= t))
    sum(pass >= k_min)
  })
  retention_data[[paste0("k", k_min)]] <- data.frame(
    threshold = threshold_seq * 100,
    retained = retention,
    criterion = paste0("≥ ", k_min, " cohorts")
  )
}
retention_df <- bind_rows(retention_data)

find_kneedle <- function(x, y, S = 1.0) {
  if (length(x) < 5) return(x[1])
  norm_x <- (x - min(x)) / (max(x) - min(x))
  norm_y <- (y - min(y)) / (max(y) - min(y))
  diff_y <- norm_y - norm_x
  norm_d <- max(diff_y) - diff_y
  threshold_vals <- norm_d - S * mean(diff(norm_x))
  idx <- which(norm_d > threshold_vals)[1]
  if (is.na(idx) || idx > length(x)) idx <- which.max(abs(diff(diff_y))) + 1
  x[idx]
}

k2_data <- retention_data[["k2"]]
optimal_thr <- find_kneedle(k2_data$threshold, k2_data$retained, S = 1.0)
n_at_optimal <- k2_data$retained[which.min(abs(k2_data$threshold - optimal_thr))]
cat(sprintf("  Kneedle optimal threshold: %.3f%%, retained: %d species\n",
            optimal_thr, n_at_optimal))

p_kneedle <- ggplot(retention_df, aes(threshold, retained, color = criterion)) +
  geom_line(linewidth = 0.5) +
  geom_vline(xintercept = optimal_thr, linetype = "dashed", color = "#E64B35", linewidth = 0.4) +
  annotate("text", x = optimal_thr, y = max(retention_df$retained) * 0.9,
           label = sprintf("Elbow: %.3f%%\n(%d species)", optimal_thr, n_at_optimal),
           hjust = -0.1, size = 2.2, color = "#E64B35", fontface = "bold") +
  scale_color_manual(values = c("≥ 2 cohorts" = "#3C5488",
                                "≥ 3 cohorts" = "#00A087",
                                "≥ 4 cohorts" = "#E64B35"), name = "Criterion") +
  scale_x_continuous(breaks = seq(0, 5, 0.5)) +
  labs(x = "Abundance threshold (%)", y = "Number of retained species",
       title = "Kneedle elbow detection") + CNS_THEME

save_fig(p_kneedle, file.path(OUT_M02, "Fig2B_Kneedle"), w = 110, h = 80)

retained_species <- filter_summary %>% filter(n_cohorts_pass >= MIN_COHORTS)
retained_species$global_mean <- rowMeans(cohort_mean_mat[retained_species$species, , drop = FALSE]) * 100
retained_species$max_cohort_mean <- apply(cohort_mean_mat[retained_species$species, , drop = FALSE], 1, max) * 100
retained_species <- retained_species %>% arrange(desc(max_cohort_mean))

filtered_species <- filter_summary %>% filter(n_cohorts_pass < MIN_COHORTS)
filtered_species$global_mean <- rowMeans(cohort_mean_mat[filtered_species$species, , drop = FALSE]) * 100
filtered_species$max_cohort_mean <- apply(cohort_mean_mat[filtered_species$species, , drop = FALSE], 1, max) * 100
filtered_species <- filtered_species %>% arrange(desc(max_cohort_mean))

save_tsv(retained_species, file.path(OUT_M02, "retained_species"))
save_tsv(filtered_species, file.path(OUT_M02, "filtered_species"))

top10_retained <- head(retained_species$species, 10)
top10_filtered <- head(filtered_species$species, 10)

build_stack_data <- function(sp_list) {
  abd_long <- as.data.frame(species_prop) %>%
    rownames_to_column("species") %>%
    filter(species %in% sp_list) %>%
    pivot_longer(-species, names_to = "Sample", values_to = "abundance") %>%
    left_join(meta[, c("Sample", "Cohort")], by = "Sample")
  abd_long$species <- factor(abd_long$species, levels = sp_list)
  abd_long
}

stack_retained <- build_stack_data(top10_retained)
stack_filtered <- build_stack_data(top10_filtered)

sp_cols_retained <- assign_sp_colors(top10_retained)
sp_cols_filtered <- assign_sp_colors(top10_filtered)

p_stack_retained <- ggplot(stack_retained, aes(Sample, abundance * 100, fill = species)) +
  geom_bar(stat = "identity", width = 1, linewidth = 0) +
  facet_grid(~ Cohort, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = sp_cols_retained, name = "Species",
                    labels = function(x) gsub("_", " ", x)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Relative abundance (%)",
       title = "Top 10 retained core species") +
  CNS_THEME +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.text = element_text(face = "italic", size = 5))

p_stack_filtered <- ggplot(stack_filtered, aes(Sample, abundance * 100, fill = species)) +
  geom_bar(stat = "identity", width = 1, linewidth = 0) +
  facet_grid(~ Cohort, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = sp_cols_filtered, name = "Species",
                    labels = function(x) gsub("_", " ", x)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Relative abundance (%)",
       title = "Top 10 filtered cohort-specific species") +
  CNS_THEME +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.text = element_text(face = "italic", size = 5))

save_fig(p_stack_retained, file.path(OUT_M02, "Fig2C_Retained_Stacked"), w = 180, h = 70)
save_fig(p_stack_filtered, file.path(OUT_M02, "Fig2D_Filtered_Stacked"), w = 180, h = 70)

profile_df <- retained_species %>%
  select(species, global_mean, max_cohort_mean, n_cohorts_pass) %>%
  arrange(desc(max_cohort_mean))
profile_df$species <- factor(profile_df$species, levels = profile_df$species)

profile_long <- profile_df %>%
  select(species, global_mean, max_cohort_mean) %>%
  pivot_longer(-species, names_to = "metric", values_to = "value")

scale_factor <- max(profile_df$max_cohort_mean) / max(profile_df$n_cohorts_pass)

p_profile <- ggplot() +
  geom_bar(data = profile_long, aes(species, value, fill = metric),
           stat = "identity", position = position_dodge(0.7), width = 0.6) +
  geom_line(data = profile_df, aes(species, n_cohorts_pass * scale_factor, group = 1),
            color = "#00A087", linewidth = 0.5) +
  geom_point(data = profile_df, aes(species, n_cohorts_pass * scale_factor),
             color = "#00A087", size = 1.5) +
  scale_fill_manual(values = c("global_mean" = "#3C5488", "max_cohort_mean" = "#E64B35"),
                    labels = c("Global mean", "Max cohort mean"), name = NULL) +
  scale_y_continuous(name = "Mean abundance (%)",
                     sec.axis = sec_axis(~./scale_factor, name = "Cohorts passing threshold")) +
  labs(x = NULL, title = "Retained core species abundance profile") +
  CNS_THEME +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 5),
        axis.title.y.right = element_text(color = "#00A087"),
        axis.text.y.right = element_text(color = "#00A087"))

save_fig(p_profile, file.path(OUT_M02, "Fig2E_AbundanceProfile"), w = 160, h = 90)

# ====================================================================
# M03 — DMM modeling (Figure 3A-B)
# ====================================================================
cat("\n╔══ M03: DMM clustering ══╗\n")
OUT_M03 <- make_outdir("M03_DMM")

abd  <- corrected
SCALE <- 10000
count_mat <- round(t(abd) * SCALE)
count_mat <- count_mat[, colSums(count_mat) > 0, drop = FALSE]
count_mat <- count_mat[rowSums(count_mat) > 0, , drop = FALSE]

set.seed(SEED)
fit_list <- lapply(1:K_MAX, function(k) {
  cat("  k=", k, "..."); f <- dmn(count_mat, k=k, verbose=FALSE)
  cat(" Lap=", round(laplace(f)), "\n"); f
})

laps <- sapply(fit_list, laplace)
bics <- sapply(fit_list, BIC)

bc_dmm <- vegdist(count_mat, method = "bray")
bc_mat <- as.matrix(bc_dmm); bc_mat[is.nan(bc_mat)] <- 0; bc_dmm <- as.dist(bc_mat)

sil_scan <- sapply(2:K_MAX, function(k) {
  a <- apply(mixture(fit_list[[k]]), 1, which.max)
  mean(silhouette(a, dist = bc_dmm)[, "sil_width"])
})
names(sil_scan) <- 2:K_MAX
sil_scan_df <- data.frame(k = 2:K_MAX, avg_sil = sil_scan)

if (!is.null(MANUAL_K)) {
  optimal_k <- MANUAL_K
} else {
  sil_top3 <- sort(sil_scan, decreasing = TRUE)[1:min(3, length(sil_scan))]
  candidates <- as.integer(names(sil_top3))
  optimal_k <- candidates[which.min(laps[candidates])]
}
cat(sprintf("  Selected k = %d\n", optimal_k))

model_sel <- data.frame(k = 1:K_MAX, Laplace = laps, BIC = bics,
                        Silhouette = c(NA, sil_scan))

best_fit <- fit_list[[optimal_k]]
post_probs <- mixture(best_fit)
ecst_assign <- apply(post_probs, 1, which.max)
comp_means <- fitted(best_fit)
colnames(comp_means) <- paste0("ECST_", 1:optimal_k)

ecst_df <- data.frame(sample_id = rownames(count_mat),
                      ECST = paste0("ECST_", ecst_assign),
                      posterior_prob = round(apply(post_probs, 1, max), 4),
                      stringsAsFactors = FALSE) %>%
  left_join(meta, by = c("sample_id" = "Sample"))
cat("  ECST distribution:\n"); print(table(ecst_df$ECST))

sil <- silhouette(ecst_assign, dist = bc_dmm)
avg_sil <- mean(sil[, "sil_width"])
sil_df <- data.frame(ECST = paste0("ECST_", sil[,"cluster"]), sil_width = sil[,"sil_width"])

ecst_colors <- setNames(ECST_PAL[1:optimal_k], paste0("ECST_", 1:optimal_k))

p3a <- ggplot(model_sel, aes(k, Laplace)) +
  geom_line(color="#3C5488", linewidth=0.5) + geom_point(color="#3C5488", size=1.5) +
  geom_point(data=model_sel[optimal_k,], color="#E64B35", size=2.5, shape=18) +
  scale_x_continuous(breaks=1:K_MAX) +
  labs(x="Number of ECSTs (k)", y="Laplace approximation",
       title="DMM model selection") + CNS_THEME

p3b <- ggplot(sil_scan_df, aes(k, avg_sil)) +
  geom_line(color="#00A087", linewidth=0.5) + geom_point(color="#00A087", size=1.5) +
  geom_point(data=sil_scan_df[sil_scan_df$k==optimal_k,], color="#E64B35", size=2.5, shape=18) +
  geom_hline(yintercept=0.5, linetype="dashed", linewidth=0.2) +
  scale_x_continuous(breaks=2:K_MAX) +
  labs(x="k", y="Avg Silhouette width",
       title=sprintf("Optimal k = %d (Silhouette = %.3f)", optimal_k, avg_sil)) +
  CNS_THEME

save_fig(p3a, file.path(OUT_M03, "Fig3A_Laplace"), w=85, h=60)
save_fig(p3b, file.path(OUT_M03, "Fig3B_Silhouette"), w=85, h=60)

save_tsv(model_sel, file.path(OUT_M03, "dmm_model_selection"))
save_tsv(ecst_df,   file.path(OUT_M03, "ecst_assignment"))
write.table(as.data.frame(comp_means), file.path(OUT_M03, "dmm_component_means.tsv"),
            sep="\t", quote=FALSE)
save_tsv(sil_scan_df, file.path(OUT_M03, "silhouette_scan"))
save_tsv(sil_df, file.path(OUT_M03, "silhouette_values"))

abd <- abd[, colnames(abd) %in% ecst_df$sample_id, drop=FALSE]
meta <- meta[rownames(meta) %in% ecst_df$sample_id, ]

# ====================================================================
# M04 — ECST composition: stacked bar with dendrogram (Figure 3C-E)
# ====================================================================
cat("\n╔══ M04: ECST composition ══╗\n")
OUT_M04 <- make_outdir("M04_Composition")

ecst_naming <- data.frame(
  ECST_id = colnames(comp_means),
  dominant = apply(comp_means, 2, function(x) rownames(comp_means)[which.max(x)]),
  dom_prop = apply(comp_means, 2, function(x) max(x)/sum(x)),
  stringsAsFactors = FALSE
)
save_tsv(ecst_naming, file.path(OUT_M04, "ecst_naming"))

core_list <- list()
for (ec in names(ecst_colors)) {
  s <- ecst_df$sample_id[ecst_df$ECST == ec]
  s <- s[s %in% colnames(abd)]
  if (length(s) < 2) next
  sub <- abd[, s, drop=FALSE]
  ma <- rowMeans(sub)*100; pv <- rowSums(sub>0)/ncol(sub)*100
  core <- data.frame(species=names(ma), mean_abd=ma, prevalence=pv, ECST=ec, stringsAsFactors=FALSE) %>%
    filter(mean_abd > 1 & prevalence > 50) %>% arrange(desc(mean_abd))
  core_list[[ec]] <- core
}
core_all <- bind_rows(core_list)
save_tsv(core_all, file.path(OUT_M04, "core_species_by_ecst"))

top10 <- names(sort(rowMeans(abd), decreasing=TRUE))[1:10]
sp_order <- c(top10, "Others")
sp_colors <- assign_sp_colors(top10)
sp_labels <- setNames(gsub("_"," ", sp_order), sp_order)
ecst_levels <- sort(unique(ecst_df$ECST))

suppressPackageStartupMessages(library(ggdendro))

tree_samples <- ecst_df$sample_id[ecst_df$sample_id %in% colnames(abd) & !is.na(ecst_df$ECST)]
sample_ecst_map <- setNames(ecst_df$ECST[match(tree_samples, ecst_df$sample_id)], tree_samples)

dend_order <- c()
dend_segments_all <- list()
dend_leaves_all <- list()
x_offset <- 0

for (ec in ecst_levels) {
  ec_samps <- names(sample_ecst_map[sample_ecst_map == ec])
  ec_samps <- ec_samps[ec_samps %in% tree_samples]
  if (length(ec_samps) < 2) {
    dend_order <- c(dend_order, ec_samps)
    if (length(ec_samps) == 1) {
      dend_leaves_all[[ec]] <- data.frame(x=x_offset+1, y=0, label=ec_samps, ECST=ec)
      x_offset <- x_offset + 1
    }
    next
  }
  bc_ec <- vegdist(t(abd[, ec_samps, drop=FALSE]), method="bray")
  bc_ec_m <- as.matrix(bc_ec); bc_ec_m[is.nan(bc_ec_m)] <- 0
  hc_ec <- hclust(as.dist(bc_ec_m), method="average")
  dend_ec <- as.dendrogram(hc_ec)
  ec_order <- hc_ec$labels[hc_ec$order]
  dd <- dendro_data(dend_ec, type="rectangle")
  seg <- segment(dd)
  seg$x    <- seg$x + x_offset
  seg$xend <- seg$xend + x_offset
  dend_segments_all[[ec]] <- seg
  lf <- dd$labels
  lf$x <- lf$x + x_offset
  lf$ECST <- ec
  dend_leaves_all[[ec]] <- lf
  dend_order <- c(dend_order, ec_order)
  x_offset <- x_offset + length(ec_samps)
}

seg_df <- bind_rows(dend_segments_all)
leaf_df <- bind_rows(dend_leaves_all)

p_dend <- ggplot(seg_df) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend), linewidth=0.15, color="grey40") +
  geom_point(data=leaf_df, aes(x=x, y=0, color=ECST), size=0.3, alpha=0.8) +
  scale_color_manual(values=ecst_colors, guide="none") +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0.002,0.002)) +
  coord_cartesian(clip="off") +
  labs(y="BC distance") +
  theme_void() +
  theme(axis.text.y=element_text(size=4, color="grey50"),
        axis.title.y=element_text(size=5, angle=90, color="grey50"),
        plot.margin=margin(2,3,0,3,"pt"))

abd_long_tree <- as.data.frame(abd[, dend_order, drop=FALSE]) %>%
  rownames_to_column("species") %>%
  pivot_longer(-species, names_to="sample_id", values_to="abundance") %>%
  left_join(ecst_df[,c("sample_id","ECST","Dataset")], by="sample_id") %>%
  filter(!is.na(ECST)) %>%
  mutate(sp_grp = ifelse(species %in% top10, species, "Others")) %>%
  group_by(sample_id, sp_grp, ECST, Dataset) %>%
  summarise(abundance=sum(abundance), .groups="drop")
abd_long_tree$sample_id <- factor(abd_long_tree$sample_id, levels=dend_order)
abd_long_tree$sp_grp    <- factor(abd_long_tree$sp_grp, levels=rev(sp_order))

p3c_tree <- ggplot(abd_long_tree, aes(sample_id, abundance*100, fill=sp_grp)) +
  geom_bar(stat="identity", width=1, linewidth=0) +
  scale_fill_manual(values=sp_colors, labels=sp_labels, name="Species",
                    guide=guide_legend(ncol=2)) +
  scale_y_continuous(expand=c(0,0), labels=function(x) paste0(x,"%")) +
  scale_x_discrete(expand=c(0.002,0.002)) +
  labs(x=NULL, y="Relative abundance") +
  CNS_THEME + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                    plot.margin=margin(0,3,0,3,"pt"),
                    legend.text=element_text(size=5, face="italic"))

actual_cohorts <- sort(unique(ecst_df$Dataset[ecst_df$sample_id %in% dend_order]))
actual_cohort_colors <- COHORT_COLORS[names(COHORT_COLORS) %in% actual_cohorts]

strip_theme_base <- theme_void() +
  theme(legend.key.size=unit(2,"mm"), legend.text=element_text(size=4),
        legend.position="right", plot.margin=margin(0,3,0,25,"pt"),
        plot.title=element_text(size=5, hjust=0, face="bold", color="grey40",
                                margin=margin(0,0,0,0)))

ecst_strip_df <- data.frame(
  sample_id=factor(dend_order, levels=dend_order),
  ECST=sample_ecst_map[dend_order], stringsAsFactors=FALSE)
p3_ecst_strip <- ggplot(ecst_strip_df, aes(sample_id, 1, fill=ECST)) +
  geom_tile(width=1) +
  scale_fill_manual(values=ecst_colors, name="ECST") +
  scale_x_discrete(expand=c(0.002,0.002)) +
  labs(y=NULL, title="ECST") + strip_theme_base

cohort_strip_df <- data.frame(
  sample_id=factor(dend_order, levels=dend_order),
  Dataset=ecst_df$Dataset[match(dend_order, ecst_df$sample_id)], stringsAsFactors=FALSE)
cohort_strip_df$Dataset <- factor(cohort_strip_df$Dataset, levels=names(actual_cohort_colors))
p3_cohort_strip <- ggplot(cohort_strip_df, aes(sample_id, 1, fill=Dataset)) +
  geom_tile(width=1) +
  scale_fill_manual(values=actual_cohort_colors, name="Cohort", drop=TRUE) +
  scale_x_discrete(expand=c(0.002,0.002)) +
  labs(y=NULL, title="Cohort") + strip_theme_base

site_strip_df <- data.frame(
  sample_id=factor(dend_order, levels=dend_order),
  SamplingSite=ecst_df$SamplingSite[match(dend_order, ecst_df$sample_id)], stringsAsFactors=FALSE)
site_strip_df$SamplingSite[is.na(site_strip_df$SamplingSite)] <- "Unknown"
actual_sites <- sort(unique(site_strip_df$SamplingSite))
site_pal <- c("#F39B7F","#8491B4","#91D1C2","#7E6148","#B09C85","#FFD700","#DAA520","#76D7C4")
SITE_COLORS <- setNames(site_pal[seq_along(actual_sites)], actual_sites)
site_strip_df$SamplingSite <- factor(site_strip_df$SamplingSite, levels=actual_sites)

p3_site_strip <- ggplot(site_strip_df, aes(sample_id, 1, fill=SamplingSite)) +
  geom_tile(width=1) +
  scale_fill_manual(values=SITE_COLORS, name="Sampling Site", drop=TRUE) +
  scale_x_discrete(expand=c(0.002,0.002)) +
  labs(y=NULL, title="Site") + strip_theme_base

group_strip_df <- data.frame(
  sample_id=factor(dend_order, levels=dend_order),
  Group=ecst_df$Group[match(dend_order, ecst_df$sample_id)], stringsAsFactors=FALSE)
group_strip_df$Group[is.na(group_strip_df$Group)] <- "Unknown"
actual_groups <- sort(unique(group_strip_df$Group))
GROUP_COLORS_STRIP <- c("Control"="#4DBBD5","MGD"="#E64B35","Pterygium"="#00A087",
                         "DED"="#F39B7F","Dry_Eye"="#F39B7F","Unknown"="#D9D9D9",
                         "Blepharitis"="#8491B4","Keratitis"="#91D1C2")
for (g in actual_groups) {
  if (!g %in% names(GROUP_COLORS_STRIP)) {
    extra <- setdiff(c("#7E6148","#B09C85","#D4AC0D","#943126","#6C3483"), GROUP_COLORS_STRIP)
    GROUP_COLORS_STRIP[g] <- extra[1]
  }
}
group_strip_df$Group <- factor(group_strip_df$Group, levels=actual_groups)

p3_group_strip <- ggplot(group_strip_df, aes(sample_id, 1, fill=Group)) +
  geom_tile(width=1) +
  scale_fill_manual(values=GROUP_COLORS_STRIP, name="Group", drop=TRUE) +
  scale_x_discrete(expand=c(0.002,0.002)) +
  labs(y=NULL, title="Group") + strip_theme_base

country_strip_df <- data.frame(
  sample_id=factor(dend_order, levels=dend_order),
  Country=ecst_df$Country[match(dend_order, ecst_df$sample_id)], stringsAsFactors=FALSE)
country_strip_df$Country[is.na(country_strip_df$Country)] <- "Unknown"
actual_countries <- sort(unique(country_strip_df$Country))
country_pal <- c("#E64B35","#3C5488","#00A087","#F39B7F","#8491B4","#D9D9D9")
COUNTRY_COLORS <- setNames(country_pal[seq_along(actual_countries)], actual_countries)
country_strip_df$Country <- factor(country_strip_df$Country, levels=actual_countries)

p3_country_strip <- ggplot(country_strip_df, aes(sample_id, 1, fill=Country)) +
  geom_tile(width=1) +
  scale_fill_manual(values=COUNTRY_COLORS, name="Country", drop=TRUE) +
  scale_x_discrete(expand=c(0.002,0.002)) +
  labs(y=NULL, title="Country") + strip_theme_base

p_panel_c_dendro <- p_dend / p3c_tree / p3_ecst_strip / p3_cohort_strip /
  p3_site_strip / p3_group_strip / p3_country_strip +
  plot_layout(heights=c(3, 12, 0.5, 0.5, 0.5, 0.5, 0.5), guides="collect") &
  theme(legend.position="right")

save_fig(p_panel_c_dendro, file.path(OUT_M04, "Fig3C_Composition_Dendro"), w=220, h=120)

ecst_valid <- ecst_df %>% filter(!is.na(ECST))
prop_df <- ecst_valid %>% dplyr::count(Dataset, ECST) %>%
  group_by(Dataset) %>% mutate(prop=n/sum(n)) %>% ungroup()
co <- ecst_valid %>% dplyr::count(Dataset) %>% arrange(desc(n)) %>% pull(Dataset)
prop_df$Dataset <- factor(prop_df$Dataset, levels=co)

cramers_v <- function(tab) {
  chi2 <- suppressWarnings(chisq.test(tab, simulate.p.value=TRUE, B=10000))
  n <- sum(tab); k <- min(nrow(tab), ncol(tab)) - 1
  if (k == 0 || n == 0) return(list(V=NA, chi2=chi2$statistic, p=chi2$p.value))
  V <- sqrt(as.numeric(chi2$statistic) / (n * k))
  list(V=V, chi2=as.numeric(chi2$statistic), p=chi2$p.value)
}

ft <- fisher.test(table(ecst_valid$Dataset, ecst_valid$ECST), simulate.p.value=TRUE, B=10000)
cv_cohort <- cramers_v(table(ecst_valid$Dataset, ecst_valid$ECST))

lbl_full <- ecst_valid %>% dplyr::count(Dataset) %>%
  mutate(l = paste0(Dataset, "\n(n=", n, ")")) %>%
  { setNames(.$l, .$Dataset) }

p3d <- ggplot(prop_df, aes(Dataset, prop, fill=ECST)) +
  geom_bar(stat="identity", width=0.7, color="white", linewidth=0.15) +
  scale_fill_manual(values=ecst_colors) +
  scale_y_continuous(expand=c(0,0), labels=scales::percent) +
  scale_x_discrete(labels=lbl_full) +
  annotate("text", x=0.5, y=0.98,
           label=sprintf("Fisher P = %.1e\nCramér's V = %.3f", ft$p.value, cv_cohort$V),
           hjust=0, vjust=1, size=2, fontface="italic") +
  labs(x=NULL, y="Proportion", title="ECST distribution per cohort") +
  CNS_THEME + AXIS_ROTATE_45

save_fig(p3d, file.path(OUT_M04, "Fig3D_ECST_perCohort"), w=95, h=70)
save_tsv(prop_df, file.path(OUT_M04, "ecst_cohort_distribution"))

core_consistency <- list()
for (ec in names(ecst_colors)) {
  dom_sp <- ecst_naming$dominant[ecst_naming$ECST_id == ec]
  if (length(dom_sp)==0 || is.na(dom_sp)) next
  if (!dom_sp %in% rownames(abd)) next
  for (coh in unique(ecst_df$Dataset)) {
    s <- ecst_df$sample_id[ecst_df$ECST == ec & ecst_df$Dataset == coh]
    s <- s[s %in% colnames(abd)]
    if (length(s) < 1) next
    core_consistency[[paste(ec, coh, dom_sp)]] <- data.frame(
      ECST=ec, cohort=coh, species=dom_sp,
      mean_abd=mean(abd[dom_sp, s])*100,
      prevalence=sum(abd[dom_sp, s]>0)/length(s)*100,
      n=length(s), stringsAsFactors=FALSE)
  }
}
core_con_df <- bind_rows(core_consistency)
save_tsv(core_con_df, file.path(OUT_M04, "core_species_consistency"))

if (nrow(core_con_df) > 0) {
  p3e <- ggplot(core_con_df, aes(x=cohort, y=ECST, size=mean_abd, color=prevalence)) +
    geom_point() +
    scale_size_continuous(range=c(1,6), name="Mean abd. (%)") +
    scale_color_gradient(low="#FDDBC7", high="#B2182B", name="Prevalence (%)") +
    labs(x=NULL, y=NULL, title="Dominant species consistency") +
    CNS_THEME + AXIS_ROTATE_45
  save_fig(p3e, file.path(OUT_M04, "Fig3E_DomSpecies_Consistency"), w=110, h=75)
}

# ====================================================================
# M05 — Robustness validation: LOO + variance decomposition (Figure 3F-G)
# ====================================================================
cat("\n╔══ M05: Robustness validation ══╗\n")
OUT_M05 <- make_outdir("M05_Validation")

ref_labels <- setNames(ecst_df$ECST, ecst_df$sample_id)

ari_calc <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  ai <- rowSums(tab); bj <- colSums(tab)
  sum_cij2 <- sum(choose(tab, 2))
  sum_ai2 <- sum(choose(ai, 2)); sum_bj2 <- sum(choose(bj, 2))
  en <- sum_ai2 * sum_bj2 / choose(n, 2)
  num <- sum_cij2 - en
  den <- (sum_ai2 + sum_bj2) / 2 - en
  if (den == 0) return(1)
  num / den
}

leave_one_level_out <- function(ecst_df_in, abd_in, group_col, optimal_k_in, SEED_in, SCALE_in) {
  results <- list()
  levels_all <- unique(ecst_df_in[[group_col]])
  levels_all <- levels_all[!is.na(levels_all) & levels_all != ""]
  if (length(levels_all) < 2) return(NULL)
  for (lv in levels_all) {
    cat(sprintf("    LOO(%s): removing '%s' ...", group_col, lv))
    keep_samps <- ecst_df_in$sample_id[ecst_df_in[[group_col]] != lv | is.na(ecst_df_in[[group_col]])]
    keep_samps <- keep_samps[keep_samps %in% colnames(abd_in)]
    if (length(keep_samps) < 20) { cat(" too few, skip\n"); next }
    abd_loo <- abd_in[, keep_samps, drop=FALSE]
    abd_loo <- abd_loo[rowSums(abd_loo) > 0, , drop=FALSE]
    count_loo <- round(t(abd_loo) * SCALE_in)
    count_loo <- count_loo[rowSums(count_loo) > 0, colSums(count_loo) > 0, drop=FALSE]
    set.seed(SEED_in)
    fit_loo <- tryCatch(dmn(count_loo, k=optimal_k_in, verbose=FALSE), error=function(e) NULL)
    if (is.null(fit_loo)) { cat(" DMM failed\n"); next }
    loo_assign <- paste0("ECST_", apply(mixture(fit_loo), 1, which.max))
    names(loo_assign) <- rownames(count_loo)
    common_s <- intersect(names(loo_assign), names(ref_labels))
    if (length(common_s) < 10) { cat(" too few common, skip\n"); next }
    ari_val <- ari_calc(ref_labels[common_s], loo_assign[common_s])
    cat(sprintf(" ARI=%.3f (n=%d)\n", ari_val, length(common_s)))
    results[[lv]] <- data.frame(
      dimension = group_col, left_out = lv,
      n_remaining = length(keep_samps), ARI = ari_val,
      stringsAsFactors = FALSE)
  }
  bind_rows(results)
}

loo_dim_results <- list()
for (gcol in c("Dataset", "SamplingSite", "Group", "Country")) {
  if (gcol %in% colnames(ecst_df) && length(unique(ecst_df[[gcol]])) >= 2) {
    cat("  ── LOO by", gcol, "──\n")
    loo_dim_results[[gcol]] <- leave_one_level_out(
      ecst_df, abd, gcol, optimal_k, SEED, SCALE)
  }
}

loo_all_dim <- bind_rows(loo_dim_results)
save_tsv(loo_all_dim, file.path(OUT_M05, "loo_results"))

if (nrow(loo_all_dim) > 0) {
  loo_all_dim$dimension <- factor(loo_all_dim$dimension,
    levels=c("Dataset","SamplingSite","Group","Country"))

  p_loo_box <- ggplot(loo_all_dim, aes(dimension, ARI, fill=dimension)) +
    geom_boxplot(alpha=0.6, width=0.5, linewidth=0.3, outlier.shape=NA) +
    geom_jitter(width=0.15, size=1.5, alpha=0.7) +
    geom_hline(yintercept=c(0.5, 0.8), linetype=c("dashed","dotted"),
               color="grey50", linewidth=0.2) +
    scale_fill_manual(values=c("Dataset"="#3C5488","SamplingSite"="#F39B7F",
                               "Group"="#4DBBD5","Country"="#E64B35")) +
    scale_y_continuous(limits=c(0,1)) +
    annotate("text", x=Inf, y=0.82, label="Good (0.8)", hjust=1.1, size=2, color="grey50") +
    annotate("text", x=Inf, y=0.52, label="Moderate (0.5)", hjust=1.1, size=2, color="grey50") +
    labs(x=NULL, y="ARI",
         title="ECST robustness across grouping dimensions") +
    CNS_THEME + theme(legend.position="none")
  save_fig(p_loo_box, file.path(OUT_M05, "Fig3F_LOO_Comparison"), w=100, h=80)
}

bc_ord_full <- vegdist(t(abd[, colnames(abd) %in% ecst_df$sample_id]), method="bray")
bc_m_full <- as.matrix(bc_ord_full); bc_m_full[is.nan(bc_m_full)] <- 0
bc_ord_full <- as.dist(bc_m_full)

meta_vd <- ecst_df[!duplicated(ecst_df$sample_id), ]
meta_vd <- meta_vd[meta_vd$sample_id %in% colnames(abd), ]
rownames(meta_vd) <- meta_vd$sample_id

vd_terms <- c("ECST", "Dataset")
if ("Country" %in% colnames(meta_vd)) vd_terms <- c(vd_terms, "Country")
if ("SamplingSite" %in% colnames(meta_vd)) vd_terms <- c(vd_terms, "SamplingSite")

bc_vd <- as.dist(as.matrix(bc_ord_full)[meta_vd$sample_id, meta_vd$sample_id])

set.seed(SEED)
fml_str <- paste("bc_vd ~", paste(vd_terms, collapse=" + "))
perm_vd <- adonis2(as.formula(fml_str), data=meta_vd, permutations=999, by="margin")
cat("  Variance decomposition (post-correction):\n"); print(perm_vd)

vd_df <- data.frame(
  term = rownames(perm_vd), R2 = perm_vd$R2, P = perm_vd$`Pr(>F)`,
  stage = "After MMUPHin", stringsAsFactors = FALSE
) %>% filter(!term %in% c("Total","Residual"))

pre_samps <- intersect(colnames(species_filt), meta_vd$sample_id)
if (length(pre_samps) > 20) {
  bc_pre_vd <- vegdist(t(species_filt[, pre_samps, drop=FALSE]), method="bray")
  bc_pre_m <- as.matrix(bc_pre_vd); bc_pre_m[is.nan(bc_pre_m)] <- 0
  bc_pre_vd <- as.dist(bc_pre_m)
  meta_vd_pre <- meta_vd[pre_samps, ]
  fml_pre <- paste("bc_pre_vd ~", paste(vd_terms, collapse=" + "))
  set.seed(SEED)
  perm_vd_pre <- adonis2(as.formula(fml_pre), data=meta_vd_pre, permutations=999, by="margin")
  cat("  Variance decomposition (pre-correction):\n"); print(perm_vd_pre)
  vd_df_pre <- data.frame(
    term = rownames(perm_vd_pre), R2 = perm_vd_pre$R2, P = perm_vd_pre$`Pr(>F)`,
    stage = "Before MMUPHin", stringsAsFactors = FALSE
  ) %>% filter(!term %in% c("Total","Residual"))
  vd_compare <- bind_rows(vd_df_pre, vd_df)
} else {
  vd_compare <- vd_df
}

vd_compare$stage <- factor(vd_compare$stage, levels=c("Before MMUPHin","After MMUPHin"))
vd_compare$sig <- ifelse(vd_compare$P < 0.001, "***",
                  ifelse(vd_compare$P < 0.01, "**",
                  ifelse(vd_compare$P < 0.05, "*", "ns")))

vd_colors <- c("ECST"="#E64B35", "Dataset"="#3C5488", "Country"="#00A087",
               "SamplingSite"="#F39B7F")

p_f4d_compare <- ggplot(vd_compare, aes(term, R2, fill=term)) +
  geom_col(width=0.6, color="black", linewidth=0.3) +
  geom_text(aes(label=sprintf("%.3f\n%s", R2, sig)),
            vjust=-0.2, size=2) +
  facet_wrap(~ stage) +
  scale_fill_manual(values=vd_colors, guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0, 0.25))) +
  labs(x=NULL, y="PERMANOVA R² (marginal)",
       title="Variance decomposition: before vs after batch correction") +
  CNS_THEME + AXIS_ROTATE_45

save_fig(p_f4d_compare, file.path(OUT_M05, "Fig3G_Variance_Decomposition"), w=160, h=85)
save_tsv(vd_compare, file.path(OUT_M05, "variance_decomposition"))

# ====================================================================
# M06 — Alpha diversity + PCoA + Slingshot (Figure 4A-B)
# ====================================================================
cat("\n╔══ M06: Alpha diversity + Ordination + Slingshot ══╗\n")
OUT_M06 <- make_outdir("M06_Diversity_Ordination")

div_df <- data.frame(sample_id=colnames(abd), Shannon=diversity(t(abd),"shannon"),
                     Richness=specnumber(t(abd)), stringsAsFactors=FALSE) %>%
  left_join(ecst_df[,c("sample_id","ECST","Dataset","Group")], by="sample_id") %>%
  filter(!is.na(ECST))

kw_shannon  <- kruskal.test(Shannon ~ ECST, data=div_df)
kw_richness <- kruskal.test(Richness ~ ECST, data=div_df)
cat(sprintf("  Kruskal-Wallis: Shannon P=%.2e | Richness P=%.2e\n",
            kw_shannon$p.value, kw_richness$p.value))

ecst_pairs <- combn(sort(unique(div_df$ECST)), 2, simplify=FALSE)
pw_results <- list()
for (pair in ecst_pairs) {
  d1 <- div_df$Shannon[div_df$ECST == pair[1]]
  d2 <- div_df$Shannon[div_df$ECST == pair[2]]
  wt <- wilcox.test(d1, d2, exact=FALSE)
  r1 <- div_df$Richness[div_df$ECST == pair[1]]
  r2 <- div_df$Richness[div_df$ECST == pair[2]]
  wt_r <- wilcox.test(r1, r2, exact=FALSE)
  pw_results[[paste(pair, collapse=" vs ")]] <- data.frame(
    group1=pair[1], group2=pair[2],
    Shannon_P=wt$p.value, Richness_P=wt_r$p.value,
    stringsAsFactors=FALSE)
}
pw_df <- bind_rows(pw_results)
pw_df$Shannon_Padj  <- p.adjust(pw_df$Shannon_P, method="BH")
pw_df$Richness_Padj <- p.adjust(pw_df$Richness_P, method="BH")
save_tsv(pw_df, file.path(OUT_M06, "ecst_pairwise_diversity"))

lm_div <- tryCatch({
  fit_s <- lm(Shannon ~ ECST + Dataset, data=div_df)
  anova_s <- anova(fit_s)
  fit_r <- lm(Richness ~ ECST + Dataset, data=div_df)
  anova_r <- anova(fit_r)
  data.frame(
    metric=c("Shannon","Richness"),
    ECST_F=c(anova_s["ECST","F value"], anova_r["ECST","F value"]),
    ECST_P=c(anova_s["ECST","Pr(>F)"], anova_r["ECST","Pr(>F)"]),
    Dataset_F=c(anova_s["Dataset","F value"], anova_r["Dataset","F value"]),
    Dataset_P=c(anova_s["Dataset","Pr(>F)"], anova_r["Dataset","Pr(>F)"]),
    stringsAsFactors=FALSE)
}, error=function(e) data.frame(metric="Shannon", ECST_F=NA, ECST_P=NA, Dataset_F=NA, Dataset_P=NA))
save_tsv(lm_div, file.path(OUT_M06, "ecst_diversity_cohort_adjusted"))

sig_pairs_shannon <- pw_df %>% filter(Shannon_Padj < 0.1) %>%
  mutate(pair = map2(group1, group2, ~c(.x, .y))) %>% pull(pair)
sig_pairs_richness <- pw_df %>% filter(Richness_Padj < 0.1) %>%
  mutate(pair = map2(group1, group2, ~c(.x, .y))) %>% pull(pair)

p_div_shannon <- ggplot(div_df, aes(ECST, Shannon, fill=ECST)) +
  geom_violin(alpha=0.4, width=0.8, linewidth=0.2, trim=FALSE, scale="width") +
  geom_boxplot(width=0.12, outlier.shape=NA, linewidth=0.2, fill="white", alpha=0.8) +
  geom_jitter(width=0.08, size=0.3, alpha=0.2, color="black") +
  scale_fill_manual(values=ecst_colors) +
  labs(x=NULL, y="Shannon diversity",
       subtitle=sprintf("Kruskal-Wallis P = %.2e", kw_shannon$p.value)) +
  CNS_THEME + theme(legend.position="none")
if (length(sig_pairs_shannon) > 0) {
  p_div_shannon <- p_div_shannon +
    stat_compare_means(comparisons=sig_pairs_shannon, method="wilcox.test",
                       label="p.signif", size=2.5, step.increase=0.08,
                       tip.length=0.01, bracket.size=0.3)
}

p_div_richness <- ggplot(div_df, aes(ECST, Richness, fill=ECST)) +
  geom_violin(alpha=0.4, width=0.8, linewidth=0.2, trim=FALSE, scale="width") +
  geom_boxplot(width=0.12, outlier.shape=NA, linewidth=0.2, fill="white", alpha=0.8) +
  geom_jitter(width=0.08, size=0.3, alpha=0.2, color="black") +
  scale_fill_manual(values=ecst_colors) +
  labs(x=NULL, y="Observed species (Richness)",
       subtitle=sprintf("Kruskal-Wallis P = %.2e", kw_richness$p.value)) +
  CNS_THEME + theme(legend.position="none")
if (length(sig_pairs_richness) > 0) {
  p_div_richness <- p_div_richness +
    stat_compare_means(comparisons=sig_pairs_richness, method="wilcox.test",
                       label="p.signif", size=2.5, step.increase=0.08,
                       tip.length=0.01, bracket.size=0.3)
}

fig_div_ecst <- (p_div_shannon | p_div_richness) +
  plot_annotation(tag_levels="a",
                  title="Alpha diversity comparison across ECST types",
                  subtitle=ifelse(!is.na(lm_div$ECST_P[1]),
                    sprintf("Cohort-adjusted ANOVA: Shannon P=%.2e, Richness P=%.2e",
                            lm_div$ECST_P[1], lm_div$ECST_P[2]), ""))

save_fig(fig_div_ecst, file.path(OUT_M06, "Fig4A_AlphaDiversity"), w=180, h=80)
save_tsv(div_df, file.path(OUT_M06, "diversity_by_ecst"))

abd_nz <- abd[, colSums(abd) > 0, drop=FALSE]
ecst_nz <- ecst_df[ecst_df$sample_id %in% colnames(abd_nz), ]
bc_ord <- vegdist(t(abd_nz), method="bray")
bc_m <- as.matrix(bc_ord); bc_m[is.nan(bc_m)] <- 0; bc_ord <- as.dist(bc_m)

pcoa_r <- cmdscale(bc_ord, k=2, eig=TRUE)
ve_ord <- round(pcoa_r$eig[1:2] / sum(abs(pcoa_r$eig)) * 100, 1)
pcoa_df <- data.frame(sample_id=colnames(abd_nz), PC1=pcoa_r$points[,1], PC2=pcoa_r$points[,2]) %>%
  left_join(ecst_nz[,c("sample_id","ECST","Dataset")], by="sample_id")

set.seed(SEED)
perm_meta <- ecst_nz[match(colnames(abd_nz), ecst_nz$sample_id), ]
perm_meta <- perm_meta[!is.na(perm_meta$ECST), ]
bc_perm <- as.dist(as.matrix(bc_ord)[perm_meta$sample_id, perm_meta$sample_id])
perm_ecst <- adonis2(bc_perm ~ ECST + Dataset, data=perm_meta, permutations=999, by="margin")
r2_ecst <- round(perm_ecst$R2[1], 3); pv_ecst <- perm_ecst$`Pr(>F)`[1]
r2_cohort <- round(perm_ecst$R2[2], 3); pv_cohort <- perm_ecst$`Pr(>F)`[2]
cat(sprintf("  PERMANOVA: ECST R²=%.3f P=%.4f; Cohort R²=%.3f P=%.4f\n",
            r2_ecst, pv_ecst, r2_cohort, pv_cohort))

hull_pcoa <- pcoa_df %>% group_by(ECST) %>% slice(chull(PC1, PC2)) %>% ungroup()
centroids_pcoa <- pcoa_df %>% group_by(ECST) %>%
  summarise(x=mean(PC1), y=mean(PC2), n=n(), .groups="drop") %>%
  mutate(label=paste0(ECST, "\n(n=", n, ")"))

p4b <- ggplot(pcoa_df, aes(PC1, PC2)) +
  geom_polygon(data=hull_pcoa, aes(fill=ECST), alpha=0.08, color=NA, show.legend=FALSE) +
  geom_polygon(data=hull_pcoa, aes(color=ECST), fill=NA, linewidth=0.3, linetype="dashed", show.legend=FALSE) +
  geom_point(aes(color=ECST), size=1.5, alpha=0.6, stroke=0.15) +
  geom_label(data=centroids_pcoa, aes(x, y, label=label, fill=ECST),
             color="white", fontface="bold", size=1.8, alpha=0.85,
             label.padding=unit(1.2,"pt"), label.size=0, show.legend=FALSE) +
  scale_color_manual(values=ecst_colors) + scale_fill_manual(values=ecst_colors) +
  labs(x=paste0("PCoA1 (", ve_ord[1], "%)"), y=paste0("PCoA2 (", ve_ord[2], "%)"),
       subtitle=sprintf("PERMANOVA: ECST R²=%.3f (P=%.3f) | Cohort R²=%.3f (P=%.3f)",
                        r2_ecst, pv_ecst, r2_cohort, pv_cohort),
       title="PCoA (Bray-Curtis)") +
  CNS_THEME

save_fig(p4b, file.path(OUT_M06, "Fig4B_PCoA"), w=110, h=85)
save_tsv(pcoa_df, file.path(OUT_M06, "pcoa_coords"))

tryCatch({
  suppressPackageStartupMessages({
    library(slingshot); library(SingleCellExperiment)
  })
  sling_samps <- pcoa_df$sample_id[!is.na(pcoa_df$ECST)]
  sling_coords <- as.matrix(pcoa_df[pcoa_df$sample_id %in% sling_samps, c("PC1","PC2")])
  rownames(sling_coords) <- sling_samps
  sling_ecst <- pcoa_df$ECST[pcoa_df$sample_id %in% sling_samps]

  abd_sling <- abd[, sling_samps, drop=FALSE]
  abd_sling <- abd_sling[rowSums(abd_sling) > 0, , drop=FALSE]

  sce <- SingleCellExperiment(
    assays = list(counts = as.matrix(round(abd_sling * SCALE))),
    reducedDims = list(PCoA = sling_coords)
  )
  colData(sce)$ECST <- factor(sling_ecst)

  ecst_sizes <- table(sling_ecst)
  start_cluster <- names(which.max(ecst_sizes))

  set.seed(SEED)
  sds <- slingshot(sce, clusterLabels = "ECST", reducedDim = "PCoA",
                   start.clus = start_cluster, stretch = 0, approx_points = 150)

  lineages <- slingLineages(sds)
  pseudotime_mat <- slingPseudotime(sds)
  curves <- slingCurves(sds)

  cat(sprintf("    Detected %d lineages:\n", length(lineages)))
  for (i in seq_along(lineages)) {
    cat(sprintf("      Lineage %d: %s\n", i, paste(lineages[[i]], collapse=" → ")))
  }

  pt_df <- data.frame(sample_id = sling_samps, ECST = sling_ecst)
  for (i in seq_len(ncol(pseudotime_mat))) {
    pt_df[[paste0("pseudotime_L", i)]] <- pseudotime_mat[, i]
  }
  pt_df <- pt_df %>% left_join(ecst_df[,c("sample_id","Dataset","Group")], by="sample_id")
  save_tsv(pt_df, file.path(OUT_M06, "slingshot_pseudotime"))

  for (lin_idx in seq_along(lineages)) {
    curve_coords <- curves[[lin_idx]]$s[curves[[lin_idx]]$ord, ]
    curve_df <- data.frame(PCoA1 = curve_coords[,1], PCoA2 = curve_coords[,2])

    pdata <- data.frame(
      PC1 = sling_coords[,1], PC2 = sling_coords[,2],
      ECST = sling_ecst, pseudotime = pseudotime_mat[, lin_idx]
    )

    p_traj <- ggplot(pdata, aes(PC1, PC2)) +
      geom_point(data = pdata[is.na(pdata$pseudotime),],
                 color="grey80", alpha=0.3, size=1) +
      geom_point(data = pdata[!is.na(pdata$pseudotime),],
                 aes(color=pseudotime), alpha=0.7, size=1.5) +
      scale_color_viridis_c(name="Pseudotime", option="C", na.value="grey80") +
      geom_path(data=curve_df, aes(x=PCoA1, y=PCoA2),
                linewidth=1.2, color="#E64B35", linetype="dashed") +
      geom_label(data=centroids_pcoa, aes(x, y, label=label),
                 fill="white", alpha=0.8, size=1.8, fontface="bold",
                 label.padding=unit(1.2,"pt"), label.size=0.2) +
      labs(x=paste0("PCoA1 (", ve_ord[1], "%)"), y=paste0("PCoA2 (", ve_ord[2], "%)"),
           title=paste0("Trajectory L", lin_idx, ": ",
                        paste(lineages[[lin_idx]], collapse=" → "))) +
      CNS_THEME

    save_fig(p_traj, file.path(OUT_M06, paste0("Fig4B_Trajectory_L", lin_idx)), w=110, h=85)
  }

  lineage_info <- data.frame(
    lineage = seq_along(lineages),
    path = sapply(lineages, paste, collapse=" → "),
    start = start_cluster, stringsAsFactors = FALSE)
  save_tsv(lineage_info, file.path(OUT_M06, "slingshot_lineages"))

}, error=function(e) {
  cat("  ⚠️ Slingshot skipped:", conditionMessage(e), "\n")
})

# ====================================================================
# M07 — Disease association (Figure 5A)
# ====================================================================
cat("\n╔══ M07: Disease association ══╗\n")
OUT_M07 <- make_outdir("M07_Disease")

disease_cohorts <- ecst_df %>%
  filter(!is.na(ECST) & !is.na(Group) & Group != "Unknown") %>%
  group_by(Dataset) %>%
  summarise(has_ctrl = "Control" %in% Group,
            has_dis  = any(!Group %in% c("Control","Unknown")),
            n_ctrl   = sum(Group == "Control"),
            n_dis    = sum(!Group %in% c("Control","Unknown")),
            diseases = paste(sort(unique(Group[Group != "Control"])), collapse="/"),
            .groups  = "drop") %>%
  filter(has_ctrl & has_dis & n_ctrl >= 3 & n_dis >= 3) %>%
  pull(Dataset)

cat("  Disease cohorts detected:\n")
for (dc in disease_cohorts) {
  dc_sub <- ecst_df %>% filter(Dataset == dc & !is.na(Group) & Group != "Unknown")
  cat(sprintf("    %s: %s\n", dc,
              paste(names(table(dc_sub$Group)), table(dc_sub$Group), sep="=", collapse=", ")))
}

if (length(disease_cohorts) > 0) {

ecst_dis <- ecst_df %>% filter(Dataset %in% disease_cohorts & Group != "Unknown" & !is.na(ECST))
ecst_dis_levels <- sort(unique(ecst_dis$ECST))
ecst_dis_colors <- ecst_colors[ecst_dis_levels]

or_list <- list()
fisher_per_cohort <- list()
for (coh in disease_cohorts) {
  sub <- ecst_dis %>% filter(Dataset==coh)
  dis <- setdiff(unique(sub$Group), "Control")
  if (length(dis)==0 || !"Control" %in% sub$Group) next
  cross <- table(sub$Group, sub$ECST); cross <- cross[,colSums(cross)>0, drop=FALSE]
  if (nrow(cross)<2) next
  ft <- fisher.test(cross, simulate.p.value=TRUE, B=10000)
  fisher_per_cohort[[coh]] <- data.frame(cohort=coh, fisher_p=ft$p.value)
  dis_label <- paste(dis, collapse="/")
  cat(sprintf("  %s: Fisher P = %.4e\n", coh, ft$p.value))
  for (ec in colnames(cross)) {
    a <- cross["Control",ec]; b <- sum(cross["Control",])-a
    cv <- sum(cross[dis, ec, drop=FALSE]); d <- sum(cross[dis, , drop=FALSE]) - cv
    ft2 <- tryCatch(fisher.test(matrix(c(cv,d,a,b),nrow=2)),
                    error=function(e) list(estimate=NA, conf.int=c(NA,NA), p.value=NA))
    or_list[[paste(coh,ec)]] <- data.frame(cohort=coh, disease=dis_label, ECST=ec,
      OR=unname(ft2$estimate), CI_lo=ft2$conf.int[1], CI_hi=ft2$conf.int[2],
      P=ft2$p.value, stringsAsFactors=FALSE)
  }
}
or_df <- bind_rows(or_list)
or_df$P_adj <- p.adjust(or_df$P, method="BH")
save_tsv(or_df, file.path(OUT_M07, "ecst_disease_OR"))
save_tsv(bind_rows(fisher_per_cohort), file.path(OUT_M07, "fisher_per_cohort"))

prop_dis <- ecst_dis %>% dplyr::count(Dataset, Group, ECST) %>%
  group_by(Dataset, Group) %>% mutate(prop=n/sum(n)) %>% ungroup()

fisher_labels <- bind_rows(fisher_per_cohort) %>%
  mutate(label = sprintf("Fisher P = %.2e", fisher_p))
facet_labels <- fisher_labels %>%
  mutate(facet_lab = paste0(cohort, "\n", label)) %>%
  { setNames(.$facet_lab, .$cohort) }

p5a <- ggplot(prop_dis, aes(Group, prop, fill=ECST)) +
  geom_bar(stat="identity", width=0.6, color="white", linewidth=0.15) +
  facet_wrap(~ Dataset, scales="free_x",
             labeller=labeller(Dataset=facet_labels)) +
  scale_fill_manual(values=ecst_dis_colors) +
  scale_y_continuous(labels=scales::percent, expand=expansion(mult=c(0, 0.05))) +
  labs(x=NULL, y="Proportion", title="ECST distribution: disease vs control") +
  CNS_THEME +
  theme(strip.text=element_text(size=6, face="bold", lineheight=1.2))

save_fig(p5a, file.path(OUT_M07, "Fig5A_Disease_ECST"), w=130, h=80)
save_tsv(prop_dis, file.path(OUT_M07, "ecst_disease_proportion"))

}

cat("\n════════════════════════════════════════════\n")
cat("✓ Pipeline finished\n")
cat("  Output:", BASE_OUT, "\n")
cat("  ECST k =", optimal_k, "| Avg Silhouette =", round(avg_sil, 3), "\n")
cat("  Species:", nrow(abd), "| Samples:", ncol(abd), "\n")
if (nrow(loo_all_dim) > 0) cat("  Mean LOO ARI =", round(mean(loo_all_dim$ARI), 3), "\n")
cat("════════════════════════════════════════════\n")
