# alpha.R
# Alpha diversity: compute, compare, and plot

source("setup.R")

outdir <- file.path(output_root, "alpha")

# ---- Control ----
use_ps_list <- ps_all_methods
use_markers <- NULL
use_lakes   <- NULL
use_metrics <- c("Observed", "InvSimpson")

# ---- Filter ----
ps_filt <- filter_ps_list(use_ps_list, use_markers, use_lakes)

# ---- Compute ----
alpha_all <- compute_alpha_all(ps_filt, use_metrics, lake_order)

alpha_long <- alpha_all %>%
  pivot_longer(cols = all_of(use_metrics),
               names_to = "Metric", values_to = "Value")

# ---- Summary Table ----
alpha_summary <- alpha_all %>%
  group_by(Marker, Lake) %>%
  summarise(across(all_of(use_metrics),
                   list(mean = mean, sd = sd), .names = "{.col}_{.fn}"),
            n = n(), .groups = "drop")

save_stats(alpha_summary,
           file.path(outdir, "stats", "alpha_summary_by_marker_lake"),
           caption = "Alpha diversity by marker and lake (mean +/- SD)")

# ---- Between Markers: Kruskal-Wallis ----
kw_mk <- run_kruskal(alpha_long, "Marker")
pw_mk <- run_pairwise_wilcox(alpha_long, "Marker")

save_stats(kw_mk, file.path(outdir, "stats", "alpha_kruskal_markers"),
           caption = "Kruskal-Wallis: alpha diversity between markers")
save_stats(pw_mk, file.path(outdir, "stats", "alpha_pairwise_markers"),
           caption = "Pairwise Wilcoxon: alpha diversity between markers")

# ---- Between Lakes Per Marker ----
kw_lk <- run_kruskal(alpha_long, "Lake", group_by_vars = "Marker")
pw_lk <- run_pairwise_wilcox(alpha_long, "Lake", group_by_vars = "Marker")

save_stats(kw_lk, file.path(outdir, "stats",
                             "alpha_kruskal_lakes_per_marker"),
           caption = "Kruskal-Wallis: alpha between lakes (per marker)")
save_stats(pw_lk %>% filter(p.adj < 0.05),
           file.path(outdir, "stats",
                     "alpha_pairwise_lakes_significant"),
           caption = "Significant pairwise lake comparisons")

# ---- Boxplot: Between Markers ----
p_mk <- ggplot(alpha_long, aes(x = Marker, y = Value, fill = Marker)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  facet_wrap(~ Metric, scales = "free_y") +
  scale_fill_manual(values = marker_colors) +
  theme_minimal(base_size = 10) +
  labs(title = "Alpha diversity - between markers",
       y = "Diversity value") +
  theme(legend.position = "bottom")
save_plot(p_mk, file.path(outdir, "figures",
                           "alpha_boxplot_between_markers.pdf"))

# ---- Boxplot: Lake by Marker ----
p_lk <- ggplot(alpha_long, aes(x = Lake, y = Value, fill = Lake)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  facet_grid(Metric ~ Marker, scales = "free_y") +
  scale_fill_manual(values = lake_colors) +
  theme_minimal(base_size = 10) +
  labs(title = "Alpha diversity - between lakes (per marker)",
       y = "Diversity value") +
  theme(legend.position = "bottom")
save_plot(p_lk, file.path(outdir, "figures",
                           "alpha_boxplot_lake_by_marker.pdf"),
          width = 16, height = 10)
