# focal_taxon.R
# Deep dive on Calanoida

source("setup.R")

outdir <- file.path(output_root, "focal_taxon")

# ---- Control ----
focal_rank  <- "Order"
focal_name  <- "Calanoida" # Rank and name should match
use_ps_list <- ps_all_methods

# ---- Detection Summary ----
focal_summary <- focal_detection_summary(use_ps_list, focal_rank,
                                          focal_name, "Species")
print(focal_summary)

save_stats(focal_summary,
           file.path(outdir, "stats",
                     paste0(tolower(focal_name), "_detection_summary")),
           caption = paste(focal_name, "detection by method"))

# ---- Distribution Across Lakes ----
focal_lake_data <- imap_dfr(use_ps_list, function(ps, m) {
  ps_sub <- subset_focal_taxon(ps, focal_rank, focal_name)
  if (is.null(ps_sub)) return(tibble())
  ps_agg <- agg_rank(ps_sub, "Speceis") %>% to_pa()
  sd <- data.frame(sample_data(ps_agg))
  sd$richness <- sample_sums(ps_agg)
  sd$Marker <- m
  sd$Sample_ID <- rownames(sd)
  sd
})

if (nrow(focal_lake_data) > 0) {
  focal_lake_data$Lake <- factor(focal_lake_data$Lake,
                                  levels = lake_order, ordered = TRUE)

  p_focal <- ggplot(focal_lake_data,
                    aes(Lake, richness, fill = Lake)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.4) +
    facet_wrap(~ Marker) +
    scale_fill_manual(values = lake_colors) +
    theme_minimal(base_size = 10) +
    labs(title = paste(focal_name, "richness by lake and marker"),
         y = paste(focal_name, "species per sample"))
  save_plot(p_focal,
            file.path(outdir, "figures", "Calanoida_richness_by_lake.pdf"))
}
