# zoop_functions.R
# Shared helper functions for zooplankton phyloseq analysis

# ---- Data Prep ----

# Make sure Lake order always remains the same
set_lake_order <- function(ps, lake_order) {
  sd <- data.frame(sample_data(ps))
  sd$Lake <- factor(sd$Lake, levels = lake_order, ordered = TRUE)
  if ("Mesh" %in% colnames(sd)) sd$Mesh <- factor(sd$Mesh)
  sample_data(ps) <- sample_data(sd)
  ps
}

# Agglomerates to a specific tax rank
agg_rank <- function(ps, rank = "Species") {
  if (rank == "ASV") return(ps)
  tt <- access(ps, "tax_table", errorIfNULL = FALSE)
  if (is.null(tt)) {
    warning("No tax_table, skipping aggregation")
    return(ps)
  }
  if (!(rank %in% rank_names(ps))) {
    warning("Rank '", rank, "' not in tax_table, skipping")
    return(ps)
  }
  tax_glom(ps, taxrank = rank, NArm = FALSE)
}

# Converts ASV table to presence/absence data
to_pa <- function(ps) {
  ot <- as(otu_table(ps), "matrix")
  ot[ot > 0] <- 1
  otu_table(ps) <- otu_table(ot, taxa_are_rows = taxa_are_rows(ps))
  ps
}

# Specify taxa list to include/exclude, and subset
subset_taxa_custom <- function(ps, tsub = NULL) {
  if (is.null(tsub)) return(ps)
  tt <- data.frame(tax_table(ps), stringsAsFactors = FALSE)
  rc <- tsub$rank
  if (!is.null(tsub$include)) {
    keep <- which(tt[[rc]] %in% tsub$include)
    ps <- prune_taxa(taxa_names(ps)[keep], ps)
  }
  if (!is.null(tsub$exclude)) {
    keep <- which(!(tt[[rc]] %in% tsub$exclude))
    ps <- prune_taxa(taxa_names(ps)[keep], ps)
  }
  prune_samples(sample_sums(ps) > 0, ps)
}

# Filters a  list of ps objects by specified variables
filter_ps_list <- function(ps_list, markers = NULL, lakes = NULL,
                           mesh = NULL, tsub = NULL) {
  if (!is.null(markers))
    ps_list <- ps_list[intersect(names(ps_list), markers)]
  lapply(ps_list, function(ps) {
    if (!is.null(lakes)) {
      sd <- data.frame(sample_data(ps))
      keep <- rownames(sd)[sd$Lake %in% lakes]
      if (length(keep) > 0) ps <- prune_samples(keep, ps)
    }
    if (!is.null(mesh)) {
      sd <- data.frame(sample_data(ps))
      keep <- rownames(sd)[sd$Mesh %in% mesh]
      if (length(keep) > 0) ps <- prune_samples(keep, ps)
    }
    if (!is.null(tsub)) ps <- subset_taxa_custom(ps, tsub)
    ps
  })
}

# ---- Output Functions ----

# Significance stars from p-value
sig_stars <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

# Save plot
save_plot <- function(p, filepath, width = 12, height = 8) {
  dir.create(dirname(filepath), recursive = TRUE, showWarnings = FALSE)
  ggsave(filepath, p, width = width, height = height)
  cat("Saved:", filepath, "\n")
}

# Save a data frame as CSV and Word - formatted table (editable)
save_stats <- function(df, filepath_base, caption = NULL) {
  dir.create(dirname(filepath_base), recursive = TRUE, showWarnings = FALSE)

  csv_path <- paste0(filepath_base, ".csv")
  write.csv(df, csv_path, row.names = FALSE)
  cat("Saved:", csv_path, "\n")

  docx_path <- paste0(filepath_base, ".docx")
  ft <- flextable(as.data.frame(df))
  ft <- autofit(ft)
  ft <- theme_booktabs(ft)
  if (!is.null(caption)) ft <- set_caption(ft, caption)
  flextable::save_as_docx(ft, path = docx_path)
  cat("Saved:", docx_path, "\n")
}

# Save text summary to file
save_summary <- function(text, filepath) {
  dir.create(dirname(filepath), recursive = TRUE, showWarnings = FALSE)
  writeLines(text, filepath)
  cat("Saved:", filepath, "\n")
}

# ---- Exploratory ----

# Overall summary function for exploration.R
summarise_ps <- function(ps, name = "dataset",
                         tax_ranks = c("Phylum", "Class", "Order",
                                       "Family", "Genus", "Species")) {
  tt <- data.frame(tax_table(ps), stringsAsFactors = FALSE)
  otu <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu <- t(otu)
  taxon_totals <- rowSums(otu)
  grand_total <- sum(taxon_totals)

  rows <- list()
  for (r in tax_ranks) {
    if (r %in% colnames(tt)) {
      n_unique <- length(unique(na.omit(tt[[r]])))
      n_na <- sum(is.na(tt[[r]]))
      pct_un <- if (grand_total > 0) {
        round(100 * sum(taxon_totals[is.na(tt[[r]])]) / grand_total, 1)
      } else {
        0
      }
      rows[[r]] <- tibble(rank = r, unique_taxa = n_unique,
                          unassigned_asvs = n_na, pct_reads_unassigned = pct_un)
    }
  }
  summary_df <- bind_rows(rows)
  summary_df$dataset <- name
  summary_df$samples <- nsamples(ps)
  summary_df$total_asvs <- ntaxa(ps)
  summary_df
}

# Print summary to consol
summarise_ps_print <- function(ps, name, tax_ranks) {
  cat("===", name, "===\n")
  cat("  Samples:", nsamples(ps), "| ASVs:", ntaxa(ps), "\n")
  df <- summarise_ps(ps, name, tax_ranks)
  for (i in seq_len(nrow(df))) {
    cat("  ", df$rank[i], ": ", df$unique_taxa[i], " unique, ",
        df$pct_reads_unassigned[i], "% reads unassigned\n", sep = "")
  }
  cat("\n")
}

# ---- Alpha Diversity ----

# Computes alpha diversity, compares with metadata variables
compute_alpha <- function(ps, marker_name,
                          metrics = c("Observed", "InvSimpson"),
                          lake_order = NULL) {
  ad <- estimate_richness(ps, measures = metrics)
  ad$Sample_ID <- sample_names(ps)
  ad$Marker <- marker_name
  meta <- data.frame(sample_data(ps))
  meta$Sample_ID <- rownames(meta)
  ad <- left_join(ad, meta, by = "Sample_ID")
  if (!is.null(lake_order))
    ad$Lake <- factor(ad$Lake, levels = lake_order, ordered = TRUE)
  ad
}

# For many markers
compute_alpha_all <- function(ps_list,
                              metrics = c("Observed", "InvSimpson"),
                              lake_order = NULL) {
  bind_rows(imap(ps_list, ~ compute_alpha(.x, .y, metrics, lake_order)))
}

# Kruskal-Wallis test across groups
run_kruskal <- function(alpha_long, group_var, group_by_vars = NULL) {
  grp <- c(group_by_vars, "Metric")
  alpha_long %>%
    group_by(across(all_of(grp))) %>%
    summarise(
      H = tryCatch(kruskal.test(Value ~ .data[[group_var]])$statistic,
                   error = function(e) NA_real_),
      df = tryCatch(kruskal.test(Value ~ .data[[group_var]])$parameter,
                    error = function(e) NA_real_),
      p_value = tryCatch(kruskal.test(Value ~ .data[[group_var]])$p.value,
                         error = function(e) NA_real_),
      .groups = "drop") %>%
    mutate(sig = sig_stars(p_value))
}

# Pairwise Wilcoxon with BH adjustment
run_pairwise_wilcox <- function(alpha_long, group_var,
                                group_by_vars = NULL) {
  grp <- c(group_by_vars, "Metric")
  alpha_long %>%
    group_by(across(all_of(grp))) %>%
    reframe({
      pw <- tryCatch(
        pairwise.wilcox.test(Value, .data[[group_var]],
                             p.adjust.method = "BH"),
        error = function(e) NULL)
      if (is.null(pw)) return(tibble())
      as.data.frame(pw$p.value) %>%
        rownames_to_column("Group1") %>%
        pivot_longer(-Group1, names_to = "Group2", values_to = "p.adj") %>%
        filter(!is.na(p.adj))
    }) %>%
    mutate(sig = sig_stars(p.adj))
}

# ---- Beta Diversity ----

# Runs ordination and returns a list with the plot, ordination, and ps
run_ordination <- function(ps, method = "NMDS", distance = "jaccard",
                           binary = TRUE, color_var = "Lake",
                           shape_var = NULL, title = "",
                           color_palette = NULL,
                           point_size = 2.5, text_size = 10) {
  ps_use <- if (binary) to_pa(ps) else ps
  ord <- ordinate(ps_use, method = method, distance = distance)

  p <- plot_ordination(ps_use, ord, color = color_var, shape = shape_var) +
    geom_point(size = point_size, alpha = 0.8) +
    theme_minimal(base_size = text_size) +
    labs(title = title)

  if (!is.null(color_palette))
    p <- p + scale_color_manual(values = color_palette)

  # 95% confidence ellipses when enough groups have >= 3 points
  grps <- data.frame(sample_data(ps_use))[[color_var]]
  grp_n <- table(grps)
  if (sum(grp_n >= 3) >= 2)
    p <- p + stat_ellipse(aes(group = .data[[color_var]]),
                          type = "t", level = 0.95, linetype = 2)

  list(plot = p, ordination = ord, ps = ps_use)
}

# Runs PERMANOVA on a ps object
run_permanova <- function(ps, formula_str, distance = "jaccard",
                          binary = TRUE, nperm = 999) {
  ps_use <- if (binary) to_pa(ps) else ps
  otu <- as(otu_table(ps_use), "matrix")
  if (taxa_are_rows(ps_use)) otu <- t(otu)
  meta <- data.frame(sample_data(ps_use))
  dm <- vegdist(otu, method = distance, binary = binary)
  adonis2(as.formula(paste("dm", formula_str)), data = meta,
          permutations = nperm)
}

# Betadisper
run_betadisper <- function(ps, group_var = "Lake",
                           distance = "jaccard", binary = TRUE) {
  ps_use <- if (binary) to_pa(ps) else ps
  otu <- as(otu_table(ps_use), "matrix")
  if (taxa_are_rows(ps_use)) otu <- t(otu)
  meta <- data.frame(sample_data(ps_use))
  dm <- vegdist(otu, method = distance, binary = binary)
  bd <- betadisper(dm, meta[[group_var]])
  list(betadisper = bd, permtest = permutest(bd, permutations = 999))
}

# ---- Differential Abundance ----

# SIMPER analysis
run_simper_analysis <- function(ps, group_var = "Lake",
                                rank = "Genus", top_n = 15, tsub = NULL) {
  ps_agg <- agg_rank(ps, rank) %>% subset_taxa_custom(tsub)
  otu <- as(otu_table(ps_agg), "matrix")
  if (taxa_are_rows(ps_agg)) otu <- t(otu)

  tt <- data.frame(tax_table(ps_agg), stringsAsFactors = FALSE)
  if (rank %in% colnames(tt)) {
    tax_names_vec <- tt[[rank]][match(colnames(otu), rownames(tt))]
    tax_names_vec[is.na(tax_names_vec)] <- colnames(otu)[is.na(tax_names_vec)]
    colnames(otu) <- make.unique(tax_names_vec)
  }

  meta <- data.frame(sample_data(ps_agg))
  sim <- simper(otu, meta[[group_var]], permutations = 99)
  summ <- lapply(names(summary(sim)), function(comp) {
    s <- summary(sim)[[comp]]
    s$taxon <- rownames(s)
    s$comparison <- comp
    head(s, top_n)
  })
  list(simper = sim, top = bind_rows(summ))
}

# ---- Heatmap ----

# Top-N taxa heatmap with lake annotation
plot_top_heatmap <- function(ps, rank = "Genus", top_n = 30,
                             annotation_var = "Lake",
                             lake_colors = NULL, lake_order = NULL,
                             tsub = NULL) {
  ps_agg <- agg_rank(ps, rank) %>% subset_taxa_custom(tsub)
  ps_rel <- transform_sample_counts(ps_agg, function(x) x / sum(x))

  otu <- as(otu_table(ps_rel), "matrix")
  if (!taxa_are_rows(ps_rel)) otu <- t(otu)

  tt <- data.frame(tax_table(ps_rel), stringsAsFactors = FALSE)
  rn <- tt[[rank]]
  rn[is.na(rn)] <- paste0("Unassigned_", seq_along(rn[is.na(rn)]))
  rownames(otu) <- make.unique(rn)

  idx <- order(rowMeans(otu), decreasing = TRUE)[1:min(top_n, nrow(otu))]
  otu_top <- otu[idx, , drop = FALSE]

  ann <- data.frame(sample_data(ps_rel))[, annotation_var, drop = FALSE]
  ann_colors <- list()
  if ("Lake" %in% colnames(ann) && !is.null(lake_colors) && !is.null(lake_order))
    ann_colors[["Lake"]] <- lake_colors[levels(factor(ann$Lake, lake_order))]

  pheatmap(otu_top, annotation_col = ann, annotation_colors = ann_colors,
           cluster_cols = TRUE, cluster_rows = TRUE,
           fontsize_row = 7, fontsize_col = 5,
           color = viridis(100), border_color = NA)
}

# ---- Venn/Upset ----

# Extract unique taxa at a given rank
get_taxa_set <- function(ps, rank = "Species") {
  ps_agg <- agg_rank(ps, rank) %>% to_pa()
  tt <- data.frame(tax_table(ps_agg), stringsAsFactors = FALSE)
  if (rank == "ASV") return(taxa_names(ps_agg))
  ids <- tt[[rank]]
  unique(ids[!is.na(ids)])
}

# Binary matrix for UpSetR from a named list of taxa sets
make_upset_matrix <- function(taxa_sets) {
  all_taxa <- unique(unlist(taxa_sets))
  mat <- data.frame(row.names = all_taxa)
  for (nm in names(taxa_sets))
    mat[[nm]] <- as.integer(all_taxa %in% taxa_sets[[nm]])
  mat
}

# ---- Trebitz Compare ----

# Read Trebitz CSV(s)
load_trebitz <- function(filepath) {
  read.csv(filepath, stringsAsFactors = FALSE)
}

# Compare marker taxa against Trebitz lists, return unexpected detections
compare_to_trebitz <- function(ps, trebitz_df, rank = "Species", lake = NULL) {
  ps_agg <- agg_rank(ps, rank) %>% to_pa()
  if (!is.null(lake)) {
    sd <- data.frame(sample_data(ps_agg))
    keep <- rownames(sd)[sd$Lake == lake]
    if (length(keep) == 0) return(tibble(taxon = character(), rank = character()))
    ps_agg <- prune_samples(keep, ps_agg)
    ps_agg <- prune_taxa(taxa_sums(ps_agg) > 0, ps_agg)
  }
  tt <- data.frame(tax_table(ps_agg), stringsAsFactors = FALSE)
  detected <- unique(na.omit(tt[[rank]]))
  known <- unique(na.omit(trebitz_df[[rank]]))
  unexpected <- setdiff(detected, known)
  if (length(unexpected) == 0) return(tibble(taxon = character(), rank = character()))
  tibble(taxon = unexpected, rank = rank)
}
