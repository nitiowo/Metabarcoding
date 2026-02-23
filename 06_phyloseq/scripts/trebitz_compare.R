# trebitz_compare.R
# Compare detected taxa against Trebitz Great Lakes lists

source("setup.R")

outdir <- file.path(output_root, "trebitz_compare")

# ---- Control ----
use_ps_list <- ps_all_methods
use_rank    <- "Species"

trebitz_file <- file.path("../data/trebitz_lists",
                          "Trebitz_Zoops_2026_overall_taxfixed.csv")
lake_cols <- c(Superior = "Superior", Michigan = "Michigan",
               Huron = "Huron", Erie = "Erie", Ontario = "Ontario")

# ---- Load and split Lists per-Lake ----
trebitz <- list()
if (file.exists(trebitz_file)) {
  raw <- read.csv(trebitz_file, stringsAsFactors = FALSE, check.names = FALSE)
  raw <- raw[, colnames(raw) != "" & !is.na(colnames(raw)), drop = FALSE]

  for (lk in names(lake_cols)) {
    col <- lake_cols[[lk]]
    if (col %in% colnames(raw)) {
      lake_taxa <- raw %>% filter(.data[[col]] == "X") %>%
        select(Species) %>% distinct()
      trebitz[[lk]] <- lake_taxa
    }
  }

  trebitz[["Overall"]] <- raw %>% select(Species) %>% distinct()
}

# ---- Find unexpected taxa ----
if (length(trebitz) > 0) {
  unexpected_all <- list()
  for (lk in intersect(lake_order, names(trebitz))) {
    for (m in names(use_ps_list)) {
      unexpected <- compare_to_trebitz(
        use_ps_list[[m]], trebitz[[lk]], rank = use_rank, lake = lk)
      if (nrow(unexpected) > 0) {
        unexpected$Lake <- lk
        unexpected$Marker <- m
        unexpected_all[[paste(lk, m)]] <- unexpected
      }
    }
  }

  if (length(unexpected_all) > 0) {
    unexpected_df <- bind_rows(unexpected_all)
    save_stats(unexpected_df,
               file.path(outdir, "stats", "unexpected_detections_per_lake"),
               caption = "Taxa detected but absent from Trebitz list (per lake)")

    unexpected_summary <- unexpected_df %>%
      group_by(Lake, Marker) %>%
      summarise(n_unexpected = n(), .groups = "drop")
    save_stats(unexpected_summary,
               file.path(outdir, "stats", "unexpected_detections_summary"),
               caption = "Count of unexpected detections per lake and marker")
  }
}
