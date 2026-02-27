# varpart.R
# Variance partitoning using spatial data vs Lake

source("setup.R")

outdir <- file.path(output_root, "varpart")

# ---- Control ----
use_ps_list   <- ps_markers
use_markers   <- NULL
env_vars      <- c("Lake")
spatial_vars  <- c("Latitude", "Longitude")
use_binary    <- TRUE

# ---- Filter ----
ps_filt <- filter_ps_list(use_ps_list, use_markers)

# ---- Per-Marker Variance Partitioning ----
vp_results <- list()
for (m in names(ps_filt)) {
  vp <- tryCatch(
    run_varpart(ps_filt[[m]], env_vars, spatial_vars, use_binary),
    error = function(e) { cat("  Error:", e$message, "\n"); NULL })

  if (!is.null(vp)) {
    vp_results[[m]] <- vp

    pdf(file.path(outdir, "figures",
                  paste0("varpart_", tolower(m), ".pdf")),
        width = 8, height = 6)
    plot(vp, bg = c("steelblue", "tomato"),
         Xnames = c("Lake + Mesh", "Lat/Lon"))
    title(main = paste(m, "- variance partitioning"))
    dev.off()
  }
}
