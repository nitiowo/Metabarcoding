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
  as.data.frame(res) %>%
    rownames_to_column("Term") %>%
    mutate(Marker = m)
}) %>% filter(!is.na(`Pr(>F)`))

write.csv(perm_jac,
          file.path(outdir, "stats", "permanova_jaccard_by_marker.csv"),
          row.names = FALSE)
