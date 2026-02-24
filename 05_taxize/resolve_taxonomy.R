# resolve_taxonomy.R
# Resolves DADA2 taxa matrices using rgbif::name_backbone(), checking a local
# cache first. Tries ranks from finest upward; accepts EXACT/HIGHERRANK ACCEPTED
# matches. Falls back to verbose=TRUE if needed.
#
# Input:  _taxa.Rds (matrix: rows=ASV seqs, cols=Kingdom...Species)
# Updated to resolve from any number of columns, starting from lowest rank.

# Output: _taxa_resolved.Rds (Kingdom, Phylum, Class, Order, Family, Genus, Species)
#         resolution_log.csv (resolution details for each input row)
#         master_taxonomy.csv (caches resolved taxonomy for future runs)

library(rgbif)
library(dplyr)
library(stringr)


# ---- Config ----

# Input .Rds files to resolve
RDS_FILES <- c(
  "04_dadaTax/dada_tax_output/LCO1490/Folmer_taxa.Rds",
  "04_dadaTax/dada_tax_output/mICOintF/Leray_taxa.Rds",
  "04_dadaTax/dada_tax_output/SSUF04/SSU_taxa.Rds"
)

# Output directory
OUT_DIR <- "05_taxize"


MASTER_SOURCES <- c("05_taxize/master_taxonomy.csv")
MASTER_OUT <- file.path(OUT_DIR, "master_taxonomy.csv")

# Resolution log
LOG_FILE <- file.path(OUT_DIR, "resolution_log.csv")

TARGET_RANKS   <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
RANK_COLS_GBIF <- tolower(TARGET_RANKS)

# ---- Master cache functions ----

MASTER_COLS <- c("query_name", TARGET_RANKS, "matchType", "confidence",
                 "usageKey", "resolved_at")

# Load and merge all master cache sources into one data frame
load_master <- function(sources) {
  dfs <- lapply(sources[file.exists(sources)], function(src) {
    df <- read.csv(src, stringsAsFactors = FALSE,
                   colClasses = c(usageKey = "character", confidence = "numeric"))
    missing <- setdiff(MASTER_COLS, colnames(df))
    if (length(missing) > 0) df[missing] <- NA
    df[, MASTER_COLS, drop = FALSE]
  })

  if (length(dfs) == 0) {
    df <- setNames(data.frame(matrix(ncol = length(MASTER_COLS), nrow = 0),
                              stringsAsFactors = FALSE), MASTER_COLS)
    df$confidence <- as.numeric(df$confidence)
    return(df)
  }

  merged <- bind_rows(dfs)
  merged[!duplicated(merged$query_name), ]
}

# Save one new entry to the cache
save_master <- function(master_df, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  write.csv(master_df, file = path, row.names = FALSE)
}

# ---- GBIF resolution ----

# Pull 7 rank fields from a name_backbone() result row
extract_gbif_ranks <- function(row) {
  ranks <- setNames(rep(NA_character_, 7), TARGET_RANKS)
  for (i in seq_along(TARGET_RANKS)) {
    col <- RANK_COLS_GBIF[i]
    if (col %in% colnames(row) && !is.na(row[[col]]) && nzchar(row[[col]]))
      ranks[TARGET_RANKS[i]] <- row[[col]]
  }
  ranks
}

# Try name_backbone(), fall back to verbose if needed
# Returns list(ranks, matchType, confidence, usageKey) or NULL
resolve_gbif <- function(query) {
  if (is.na(query) || !nzchar(trimws(query))) return(NULL)

  result <- tryCatch(name_backbone(name = query), error = function(e) NULL)

  if (!is.null(result) && is.data.frame(result) && nrow(result) > 0) {
    mt <- result$matchType[1]
    st <- if ("status" %in% colnames(result)) result$status[1] else NA
    if (!is.na(mt) && mt %in% c("EXACT", "HIGHERRANK") &&
        !is.na(st) && st == "ACCEPTED") {
      return(list(
        ranks = extract_gbif_ranks(result[1, ]),
        matchType = mt,
        confidence = as.numeric(result$confidence[1]),
        usageKey = as.character(result$usageKey[1])
      ))
    }
  }

  # Verbose fallback: pick best ACCEPTED candidate
  vb <- tryCatch(name_backbone(name = query, verbose = TRUE), error = function(e) NULL)
  if (is.null(vb) || !is.data.frame(vb) || nrow(vb) == 0) return(NULL)
  if (!"status" %in% colnames(vb)) return(NULL)

  accepted <- vb %>%
    filter(!is.na(matchType), matchType %in% c("EXACT", "HIGHERRANK"),
           !is.na(status), status == "ACCEPTED") %>%
    arrange(desc(matchType == "EXACT"), desc(confidence))

  if (nrow(accepted) == 0) return(NULL)

  best <- accepted[1, ]
  list(
    ranks = extract_gbif_ranks(best),
    matchType = best$matchType,
    confidence = as.numeric(best$confidence),
    usageKey = as.character(best$usageKey)
  )
}

# ---- Name cleaning ----

# Exception handling
# Decide whether to attempt resolution or clean the name
classify_name <- function(val) {
  if (is.na(val) || !nzchar(trimws(val)))
    return(list(query = NA, action = "skip"))

  s <- trimws(val)

  # Environmental / uncultured
  if (str_detect(s, "(?i)^(uncultured|unidentified|unknown|metagenome|environmental)"))
    return(list(query = NA, action = "skip"))

  # Numeric / accession-like
  if (str_detect(s, "^[A-Z]{0,2}[0-9]{3,}$") || str_detect(s, "^[0-9]+$"))
    return(list(query = NA, action = "skip"))

  # Species only (lowercase first letter) - needs genus prepended
  if (str_detect(s, "^[a-z]"))
    return(list(query = s, action = "epithet"))

  # Subfamily / tribe names ending in inae or ini
  if (str_detect(s, "(inae|ini)$"))
    return(list(query = s, action = "subfamily"))

  # Normal capitalized name
  list(query = s, action = "use")
}

# ---- Resolution loop ----

# Resolve unique taxonomy combos, update master cache, return resolved matrix
resolve_rds <- function(rds_path, master_df, log_rows) {

  marker <- tools::file_path_sans_ext(basename(rds_path))

  taxa_mat <- readRDS(rds_path)
  taxa_df  <- as.data.frame(taxa_mat, stringsAsFactors = FALSE)

  # Verify expected columns
  if (!all(TARGET_RANKS %in% colnames(taxa_df))) {
    missing <- setdiff(TARGET_RANKS, colnames(taxa_df))
    stop(sprintf("Missing columns in %s: %s", rds_path, paste(missing, collapse = ", ")))
  }

  # Build unique taxonomy combos
  combo_keys <- apply(taxa_df[, TARGET_RANKS], 1, function(r) paste(r, collapse = "|"))
  uniq_keys  <- unique(combo_keys)

  # Cache for this file's unique combos
  combo_cache <- vector("list", length(uniq_keys))
  names(combo_cache) <- uniq_keys

  for (i in seq_along(uniq_keys)) {
    key <- uniq_keys[i]
    idx <- which(combo_keys == key)[1]
    rv  <- setNames(as.character(taxa_df[idx, TARGET_RANKS]), TARGET_RANKS)

    # Try ranks finest to coarsest
    resolved <- FALSE
    final_ranks <- setNames(rep(NA_character_, 7), TARGET_RANKS)

    for (rank_idx in length(TARGET_RANKS):1) {
      rk  <- TARGET_RANKS[rank_idx]
      val <- rv[rk]
      cl  <- classify_name(val)

      if (cl$action == "skip") next
      
      # Build query
      query <- NULL
      if (cl$action == "epithet" && rank_idx == 7) {
        genus_cl <- classify_name(rv["Genus"])
        if (genus_cl$action == "use") {
          query <- paste(genus_cl$query, cl$query)
        } else {
          next
        }
      } else {
        query <- cl$query
      }

      if (is.null(query) || is.na(query)) next

      # Check master cache first
      cache_idx <- match(query, master_df$query_name)
      if (!is.na(cache_idx)) {
        mr <- master_df[cache_idx, ]
        final_ranks <- setNames(as.character(unlist(mr[, TARGET_RANKS])), TARGET_RANKS)
        resolved <- TRUE

        log_rows[[length(log_rows) + 1]] <- data.frame(
          marker = marker, combo_key = substr(key, 1, 80),
          attempt_rank = rk, query = query,
          outcome = "cache_hit", matchType = mr$matchType,
          confidence = mr$confidence,
          stringsAsFactors = FALSE
        )
        break
      }

      # Query GBIF
      gbif_result <- resolve_gbif(query)

      if (is.null(gbif_result)) {
        log_rows[[length(log_rows) + 1]] <- data.frame(
          marker = marker, combo_key = substr(key, 1, 80),
          attempt_rank = rk, query = query,
          outcome = "no_match", matchType = NA, confidence = NA,
          stringsAsFactors = FALSE
        )
        next
      }

      final_ranks <- gbif_result$ranks
      resolved <- TRUE

      log_rows[[length(log_rows) + 1]] <- data.frame(
        marker = marker, combo_key = substr(key, 1, 80),
        attempt_rank = rk, query = query,
        outcome = "resolved", matchType = gbif_result$matchType,
        confidence = gbif_result$confidence,
        stringsAsFactors = FALSE
      )

      # Add to master cache
      new_row <- data.frame(
        query_name = query,
        Kingdom = gbif_result$ranks["Kingdom"],
        Phylum = gbif_result$ranks["Phylum"],
        Class = gbif_result$ranks["Class"],
        Order = gbif_result$ranks["Order"],
        Family = gbif_result$ranks["Family"],
        Genus = gbif_result$ranks["Genus"],
        Species = gbif_result$ranks["Species"],
        matchType = gbif_result$matchType,
        confidence = gbif_result$confidence,
        usageKey = gbif_result$usageKey,
        resolved_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        stringsAsFactors = FALSE
      )
      master_df <- bind_rows(master_df, new_row)
      break
    }

    # If not resolved even at Kingdom level, keep original values and inspect later
    if (!resolved) {
      idx <- which(!is.na(rv) & nzchar(trimws(rv)))
      if (length(idx) > 0) final_ranks[idx] <- rv[TARGET_RANKS[idx]]
      log_rows[[length(log_rows) + 1]] <- data.frame(
        marker = marker, combo_key = substr(key, 1, 80),
        attempt_rank = "EXHAUSTED", query = NA,
        outcome = "unresolved", matchType = NA, confidence = NA,
        stringsAsFactors = FALSE
      )
    }

    combo_cache[[key]] <- final_ranks

    # Break after every 10 API calls
    if (i %% 10 == 0) Sys.sleep(0.5)
  }

  # Build resolved matrix
  resolved_mat <- t(sapply(combo_keys, function(k) combo_cache[[k]]))
  colnames(resolved_mat) <- TARGET_RANKS
  rownames(resolved_mat) <- rownames(taxa_df)

  resolved_mat[resolved_mat == "Metazoa"] <- "Animalia"

  out_name <- sub("(\\.Rds)$", "_resolved\\1", basename(rds_path), ignore.case = TRUE)
  out_path <- file.path(OUT_DIR, out_name)
  saveRDS(resolved_mat, out_path)
  message(sprintf("  Output: %s (%d ASVs)", out_path, nrow(resolved_mat)))

  list(master_df = master_df, log_rows = log_rows)
}

# ---- Main ----

master_df <- load_master(MASTER_SOURCES)

dir.create(OUT_DIR, showWarnings = FALSE)
log_rows <- list()

for (rds_file in RDS_FILES) {
  if (!file.exists(rds_file)) next
  result    <- resolve_rds(rds_file, master_df, log_rows)
  master_df <- result$master_df
  log_rows  <- result$log_rows
}

save_master(master_df, MASTER_OUT)

if (length(log_rows) > 0) {
  log_df <- bind_rows(log_rows)
  write.csv(log_df, LOG_FILE, row.names = FALSE)

  # Summary
  log_df %>%
    group_by(marker, outcome) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(marker, outcome) %>%
    print()
}
