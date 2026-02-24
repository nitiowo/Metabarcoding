# Taxonomic cleanup

Using `rgbif` R package to resolve DADA2 taxonomy against the GBIF backbone.

**`resolve_taxonomy.R:`** Resolves taxonomy from `04_dadaTax` `.Rds` files (matrix: rows = ASV sequences,
cols = Kingdom...Species) to standardized 7-rank taxonomy via `rgbif::name_backbone()`.

**Strategy:** tries from finest rank upward, checks `master_taxonomy.csv` before any API call so each name is only queried once.

**Outputs:**
- `xxx_taxa_resolved.Rds` — same matrix format as input, GBIF-resolved
- `master_taxonomy.csv` — name cache
- `resolution_log.csv` —  resolution details per input row