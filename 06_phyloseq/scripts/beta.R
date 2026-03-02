# beta.R
# Beta diversity: ordination and PERMANOVA

source("setup.R")

outdir <- file.path(output_root, "beta")

# ---- Control ----
use_ps_list <- ps_all_methods
use_markers <- NULL
use_lakes   <- NULL

# ---- Filter ----
ps_filt <- filter_ps_list(use_ps_list, use_markers, use_lakes)
marker_names <- names(ps_filt)

# ---- Within-Marker Ordinations ----
for (mk in marker_names) {
  ps <- ps_filt[[mk]]

  for (dist_method in c("jaccard", "bray")) {
    dist_mat <- phyloseq::distance(ps, method = dist_method)
    ord <- ordinate(ps, method = "NMDS", distance = dist_mat, trymax = 50)

    p <- plot_ordination(ps, ord, color = "Lake") +
      geom_point(size = 2.5) +
      scale_color_manual(values = lake_colors) +
      theme_minimal(base_size = 10) +
      labs(title = paste(mk, "-", dist_method, "NMDS"))

    save_plot(p, file.path(outdir, "figures",
                            paste0("beta_NMDS_", mk, "_", dist_method, ".pdf")))

    # PERMANOVA
    sdf <- data.frame(sample_data(ps))
    perm <- vegan::adonis2(dist_mat ~ Lake, data = sdf, permutations = 999)

    save_stats(tidy_permanova(perm),
               file.path(outdir, "stats",
                          paste0("beta_permanova_", mk, "_", dist_method)),
               caption = paste("PERMANOVA:", mk, dist_method))

    # Betadisper
    bd <- run_betadisper(dist_mat, sdf$Lake)
    save_stats(bd, file.path(outdir, "stats",
                              paste0("beta_betadisper_", mk, "_",
                                     dist_method)),
               caption = paste("Betadisper:", mk, dist_method))
  }
}

# ---- Between-Marker Distance Jaccard PcoA ----
ps_mk <- build_marker_ps(ps_filt, tax_level = "Genus")
dist_mk <- phyloseq::distance(ps_mk, method = "jaccard")
ord_mk <- ordinate(ps_mk, "PCoA", distance = dist_mk)

p_mk <- plot_ordination(ps_mk, ord_mk, color = "Marker", shape = "Lake") +
  geom_point(size = 3) +
  scale_color_manual(values = marker_colors) +
  theme_minimal(base_size = 10) +
  labs(title = "Between-marker PCoA (Jaccard, Genus level)")
save_plot(p_mk, file.path(outdir, "figures",
                           "beta_PCoA_between_markers.pdf"),
          width = 10, height = 8)

# PERMANOVA (marker as predictor)
sdf_mk <- data.frame(sample_data(ps_mk))
perm_mk <- vegan::adonis2(dist_mk ~ Marker, data = sdf_mk,
                            permutations = 999)
save_stats(tidy_permanova(perm_mk),
           file.path(outdir, "stats", "beta_permanova_marker_effect"),
           caption = "PERMANOVA: marker effect on community composition")

# ---- Combined Markers Ordinations ----
for (ps_comb_name in c("ps_markers_combined", "ps_coi_combined")) {
  ps_c <- get(ps_comb_name, envir = .GlobalEnv)
  label <- gsub("ps_", "", ps_comb_name) %>% gsub("_combined", "", .)

  for (dist_method in c("jaccard", "bray")) {
    dist_c <- phyloseq::distance(ps_c, method = dist_method)
    ord_c <- ordinate(ps_c, "NMDS", distance = dist_c, trymax = 50)

    p_c <- plot_ordination(ps_c, ord_c, color = "Lake") +
      geom_point(size = 2.5) +
      scale_color_manual(values = lake_colors) +
      theme_minimal(base_size = 10) +
      labs(title = paste("Combined", label, "-", dist_method, "NMDS"))

    save_plot(p_c, file.path(outdir, "figures",
                              paste0("beta_NMDS_combined_", label, "_",
                                     dist_method, ".pdf")))

    sdf_c <- data.frame(sample_data(ps_c))
    perm_c <- vegan::adonis2(dist_c ~ Lake, data = sdf_c,
                              permutations = 999)
    save_stats(tidy_permanova(perm_c),
               file.path(outdir, "stats",
                          paste0("beta_permanova_combined_", label, "_",
                                 dist_method)),
               caption = paste("PERMANOVA: combined", label, dist_method))
  }
}
