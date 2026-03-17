#!/usr/bin/env Rscript
################################################################################
# METADATA MATCHING UTILITIES
# Multi-strategy genome ID reconciliation functions
################################################################################

#' Batch match genome IDs to metadata names
#'
#' This function uses multiple strategies to match genome IDs from analysis
#' files to sample names in metadata files. It tries progressively more
#' aggressive transformations until a match is found.
#'
#' @param genome_ids Character vector of genome IDs from analysis files
#' @param meta_names Character vector of sample names from metadata
#' @return Character vector of matched metadata names (NA for unmatched)
#'
#' @details
#' Matching strategies (in order):
#' 1. Direct match
#' 2. Underscore to space
#' 3. Underscore to dot
#' 4. Underscore to forward slash
#' 5. Underscore to pipe
#' 6. Underscore to plus
#' 7. Remove trailing underscores
#' 8. ESBL/AmpC normalization
#' 9. Three-part ID format (A_B_C -> A B:C)
#' 10. Double underscore format (A__B -> A (B))
#' 11. Numeric ID (strip leading zeros)
#' 12. Alphanumeric key matching (remove all punctuation)
#' 13. Final dot substitution
#'
#' @examples
#' genome_ids <- c("Sample_001", "Sample_002_X")
#' meta_names <- c("Sample 001", "Sample.002.X")
#' batch_match(genome_ids, meta_names)
batch_match <- function(genome_ids, meta_names) {
  n      <- length(genome_ids)
  result <- rep(NA_character_, n)
  todo   <- seq_len(n)

  # Helper: resolve matches and update todo list
  resolve <- function(idx, candidates) {
    hit <- !is.na(candidates) & (candidates %in% meta_names)
    result[idx[hit]] <<- candidates[hit]
    idx[!hit]
  }
  
  # Helper: get current unmatched IDs
  ids <- function() genome_ids[todo]

  # Strategy 1: Direct match
  todo <- resolve(todo, ids())
  
  # Strategy 2: Underscore to space
  todo <- resolve(todo, gsub("_", " ", ids(), fixed = TRUE))
  
  # Strategy 3: Underscore to dot
  todo <- resolve(todo, gsub("_", ".", ids(), fixed = TRUE))
  
  # Strategy 4: Underscore to forward slash
  todo <- resolve(todo, gsub("_", "/", ids(), fixed = TRUE))
  
  # Strategy 5: Underscore to pipe
  todo <- resolve(todo, gsub("_", "|", ids(), fixed = TRUE))
  
  # Strategy 6: Underscore to plus
  todo <- resolve(todo, gsub("_", "+", ids(), fixed = TRUE))
  
  # Strategy 7: Remove trailing underscores
  todo <- resolve(todo, sub("_+$", "", ids()))
  
  # Strategy 8: ESBL/AmpC normalization
  {
    v <- gsub("_", " ", ids(), fixed = TRUE)
    v <- gsub("ESBL AmpC", "ESBL/AmpC", v, ignore.case = TRUE)
    todo <- resolve(todo, v)
  }
  
  # Strategy 9: Three-part format (A_B_C -> A B:C)
  {
    v <- vapply(ids(), function(gid) {
      p <- strsplit(gid, "_", fixed = TRUE)[[1]]
      if (length(p) == 3L) paste0(p[1], " ", p[2], ":", p[3]) else NA_character_
    }, character(1L), USE.NAMES = FALSE)
    todo <- resolve(todo, v)
  }
  
  # Strategy 10: Double underscore format (A__B -> A (B))
  {
    v <- vapply(ids(), function(gid) {
      if (!grepl("__", gid, fixed = TRUE)) return(NA_character_)
      p    <- strsplit(gid, "__", fixed = TRUE)[[1]]
      rest <- gsub("_", ".", sub("_+$", "", p[2]), fixed = TRUE)
      paste0(p[1], " (", rest, ")")
    }, character(1L), USE.NAMES = FALSE)
    todo <- resolve(todo, v)
  }
  
  # Strategy 11: Numeric IDs (strip leading zeros)
  {
    v <- ifelse(grepl("^\\d+$", ids()), sub("^0+", "", ids()), NA_character_)
    v[v == ""] <- NA_character_
    todo <- resolve(todo, v)
  }
  
  # Strategy 12: Alphanumeric key matching
  {
    meta_key <- gsub("[^A-Za-z0-9]", "", meta_names)
    dup      <- duplicated(meta_key)
    meta_map <- setNames(meta_names[!dup], meta_key[!dup])
    v_keys   <- gsub("[^A-Za-z0-9]", "", ids())
    todo     <- resolve(todo, unname(meta_map[v_keys]))
  }
  
  # Strategy 13: Final dot substitution
  todo <- resolve(todo, sub("_([^_]*)$", ".\\1", ids()))
  
  result
}

#' Calculate matching statistics
#'
#' @param matched_names Character vector from batch_match
#' @param st Sequence type name
#' @return Tibble with matching statistics
get_matching_stats <- function(matched_names, st) {
  total   <- length(matched_names)
  matched <- sum(!is.na(matched_names))
  rate    <- if (total > 0) (matched / total * 100) else 0
  
  tibble::tibble(
    ST = st,
    Total = total,
    Matched = matched,
    Unmatched = total - matched,
    Match_Rate = round(rate, 1)
  )
}

#' Print matching summary
#'
#' @param stats_df Tibble with matching statistics
#' @return Invisible NULL
print_matching_summary <- function(stats_df) {
  cat("\n=== METADATA MATCHING SUMMARY ===\n\n")
  print(as.data.frame(stats_df))
  
  overall_rate <- with(stats_df, sum(Matched) / sum(Total) * 100)
  cat(sprintf("\nOverall match rate: %.1f%%\n", overall_rate))
  
  if (any(stats_df$Match_Rate < 90)) {
    low_match <- stats_df %>% 
      dplyr::filter(Match_Rate < 90) %>% 
      dplyr::pull(ST)
    cat(sprintf("\n⚠ Warning: Low match rate for: %s\n", 
                paste(low_match, collapse = ", ")))
  }
  
  invisible(NULL)
}
