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

# ---- Within-Marker Ordinations: Jaccard P/A. Use Bray-Curtis for abundance ----
ord_jac <- imap(ps_filt, ~ run_ordination(
  .x, "NMDS", "jaccard", TRUE, "Lake",
  title = paste(.y, "- NMDS Jaccard"),
  color_palette = lake_colors))

p_jac <- wrap_plots(map(ord_jac, "plot"), ncol = 2) +
  plot_annotation(title = "NMDS Jaccard (P/A)")
save_plot(p_jac, file.path(outdir, "figures",
                            "nmds_jaccard_all_markers.pdf"),
          width = 16, height = 14)

# ---- PERMANOVA: lake effect (Jaccard P/A) ----
perm_jac <- imap_dfr(ps_filt, function(ps, m) {
  res <- run_permanova(ps, "~ Lake", "jaccard", TRUE)
  as.data.frame(res) %>% rownames_to_column("Term") %>% mutate(Marker = m)
})
perm_jac_clean <- perm_jac %>% filter(!is.na(`Pr(>F)`))
save_stats(perm_jac_clean,
           file.path(outdir, "stats", "permanova_jaccard_by_marker"),
           caption = "PERMANOVA - lake effect (Jaccard P/A)")

# ---- Within-Marker Ordinations: Bray-Curtis Abundance ----
ord_bray <- imap(ps_filt, ~ run_ordination(
  .x, "NMDS", "bray", FALSE, "Lake",
  title = paste(.y, "- NMDS Bray-Curtis"),
  color_palette = lake_colors))

p_bray <- wrap_plots(map(ord_bray, "plot"), ncol = 2) +
  plot_annotation(title = "NMDS Bray-Curtis (abundance)")
save_plot(p_bray, file.path(outdir, "figures",
                             "nmds_bray_all_markers.pdf"),
          width = 16, height = 14)

# PERMANOVA - Bray-Curtis
perm_bray <- imap_dfr(ps_filt, function(ps, m) {
  res <- run_permanova(ps, "~ Lake", "bray", FALSE)
  as.data.frame(res) %>% rownames_to_column("Term") %>% mutate(Marker = m)
})
perm_bray_clean <- perm_bray %>% filter(!is.na(`Pr(>F)`))
save_stats(perm_bray_clean,
           file.path(outdir, "stats", "permanova_bray_by_marker"),
           caption = "PERMANOVA - lake effect (Bray-Curtis)")

# ---- Betadisper ----
bd_results <- imap_dfr(ps_filt, function(ps, m) {
  bd <- run_betadisper(ps, "Lake", "jaccard", TRUE)
  pt <- bd$permtest
  tibble(Marker = m, F_stat = pt$tab$F[1],
         p_value = pt$tab$`Pr(>F)`[1])
}) %>% mutate(sig = sig_stars(p_value))

save_stats(bd_results,
           file.path(outdir, "stats", "betadisper_jaccard"),
           caption = "Betadisper - homogeneity of dispersion (Jaccard)")
