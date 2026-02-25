# setup.R
# Load all data, build phyloseq objects, define palettes

library(phyloseq)
library(tidyverse)
library(vegan)
library(patchwork)
library(RColorBrewer)
library(flextable)
library(pheatmap)
library(viridis)
library(VennDiagram)

set.seed(42)
output_root <- "output"

source("zoop_functions.R")

# ---- Palettes and constants ----
lake_order <- c("Superior", "Michigan", "Huron", "Erie", "Ontario")

lake_colors <- c(
  Superior = "#1b9e77", Michigan = "#d95f02", Huron = "#7570b3",
  Erie     = "#e7298a", Ontario  = "#66a61e"
)

marker_colors <- c(
  Folmer = "#e41a1c", Leray = "#377eb8",
  `18S`  = "#4daf4a", Morphology = "#984ea3"
)

tax_ranks <- c("Kingdom", "Phylum", "Class", "Order",
               "Family", "Genus", "Species")

alpha_metrics <- c("Observed", "InvSimpson")

# ---- Load phyloseq objects ----
ps_folmer <- readRDS("setup_output/folmer_ps.RDS")
ps_leray  <- readRDS("setup_output/leray_ps.RDS")
ps_18S    <- readRDS("setup_output/ssu_ps.RDS")
ps_morph  <- readRDS("setup_output/morph_ps.RDS")

# ---- Set lake ordering ----
ps_folmer <- set_lake_order(ps_folmer, lake_order)
ps_leray  <- set_lake_order(ps_leray, lake_order)
ps_18S    <- set_lake_order(ps_18S, lake_order)
ps_morph  <- set_lake_order(ps_morph, lake_order)

# ---- Named lists ----
ps_markers     <- list(Folmer = ps_folmer, Leray = ps_leray, `18S` = ps_18S)
ps_all_methods <- c(ps_markers, list(Morphology = ps_morph))
ps_coi         <- list(Folmer = ps_folmer, Leray = ps_leray)

# ---- Output directories for each analysis step/script ----
out_dirs <- c("exploratory", "alpha", "beta", "composition",
              "overlap", "heatmaps", "trebitz_compare")
for (d in out_dirs) {
  dir.create(file.path(output_root, d, "figures"),
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_root, d, "stats"),
             recursive = TRUE, showWarnings = FALSE)
}
