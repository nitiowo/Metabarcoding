# heatmaps.R
# Top-N taxa heatmaps per marker

source("setup.R")

outdir <- file.path(output_root, "heatmaps")

# ---- Control ----
use_ps_list <- ps_markers
use_rank    <- "Genus"
top_n       <- 30

# ---- Heatmap Per Marker ----
for (m in names(use_ps_list)) {
  fname <- file.path(outdir, "figures",
                     paste0("heatmap_", tolower(use_rank), "_",
                            tolower(m), ".pdf"))
  dir.create(dirname(fname), recursive = TRUE, showWarnings = FALSE)
  pdf(fname, width = 14, height = 10)
  plot_top_heatmap(use_ps_list[[m]], rank = use_rank, top_n = top_n,
                   annotation_var = "Lake",
                   lake_colors = lake_colors, lake_order = lake_order)
  dev.off()
}
