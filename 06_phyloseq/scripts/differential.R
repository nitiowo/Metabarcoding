# differential.R
# SIMPER analysis per marker

source("setup.R")

outdir <- file.path(output_root, "differential")

# ---- Control ----
use_ps_list <- ps_markers
use_markers <- NULL
use_rank    <- "Genus"
top_n       <- 15

# ---- Filter ----
ps_filt <- filter_ps_list(use_ps_list, use_markers)

# ---- SIMPER Per Marker ----
for (m in names(ps_filt)) {
  sim_res <- run_simper_analysis(ps_filt[[m]], "Lake",
                                 rank = use_rank, top_n = top_n)
  save_stats(sim_res$top,
             file.path(outdir, "stats",
                       paste0("simper_top_", tolower(m))),
             caption = paste("SIMPER top contributing", use_rank, "-", m))
}
