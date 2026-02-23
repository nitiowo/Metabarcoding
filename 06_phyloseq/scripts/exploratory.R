# exploratory.R
# Dataset summaries

source("setup.R")

outdir <- file.path(output_root, "exploratory")

# ---- Summarize All Datasets ----
datasets <- ps_all_methods

summary_df <- imap_dfr(datasets, ~ summarise_ps(.x, .y, tax_ranks))

summary_df <- summary_df %>%
  select(dataset, samples, total_asvs, rank, unique_taxa,
         unassigned_asvs, pct_reads_unassigned)

# ---- Print to console and save----
for (nm in names(datasets)) {
  summarise_ps_print(datasets[[nm]], nm, tax_ranks)
}

save_stats(summary_df,
           file.path(outdir, "stats", "dataset_summary"),
           caption = "Dataset summary: taxa counts and % reads unassigned")

txt <- capture.output({
  for (nm in names(datasets))
    summarise_ps_print(datasets[[nm]], nm, tax_ranks)
})
save_summary(txt, file.path(outdir, "stats", "dataset_summary.txt"))
