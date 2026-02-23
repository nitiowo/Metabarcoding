# overlap.R
# Taxa overlap between markers

source("setup.R")

outdir <- file.path(output_root, "overlap")

# ---- Control ----
use_ps_list <- ps_markers
use_ranks   <- c("Species", "Genus", "Family")

# ---- Get overlap sets for each rank ----
for (rank in use_ranks) {

  taxa_sets <- imap(use_ps_list, ~ get_taxa_set(.x, rank))

  # Get set sizes
  size_df <- imap_dfr(taxa_sets, ~ tibble(Marker = .y, rank = rank,
                                           n_taxa = length(.x)))
  print(size_df)

  # ---- Venn Diagrams ----
  if (length(taxa_sets) >= 2 && length(taxa_sets) <= 5) {
    fc <- unname(marker_colors[names(taxa_sets)])
    venn_obj <- venn.diagram(
      x = taxa_sets, category.names = names(taxa_sets),
      filename = NULL, output = TRUE, imagetype = "none",
      col = fc, fill = adjustcolor(fc, alpha.f = 0.3),
      cat.col = fc, cat.cex = 1.2, margin = 0.1,
      main = paste("Marker overlap -", rank))

    venn_path <- file.path(outdir, "figures",
                           paste0("venn_", tolower(rank), ".pdf"))
    dir.create(dirname(venn_path), recursive = TRUE, showWarnings = FALSE)
    pdf(venn_path, width = 10, height = 8)
    grid::grid.newpage()
    grid::grid.draw(venn_obj)
    dev.off()
  }
}
