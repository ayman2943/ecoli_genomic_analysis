#!/usr/bin/env Rscript
################################################################################
# AMR ANALYSIS MODULE 1: METADATA MATCHING
# Multi-strategy genome ID reconciliation and data integration
################################################################################

#' Run metadata matching for all sequence types
#'
#' @param STs Character vector of sequence types
#' @param CARD_DIR Path to CARD summary directory
#' @param META_DIR Path to metadata directory
#' @return List with matched data and summary statistics
run_metadata_matching <- function(STs, CARD_DIR, META_DIR) {
  
  all_card_data    <- list()
  matching_stats   <- list()
  
  for (st in STs) {
    cat(sprintf("### Processing %s ###\n", st))
    
    # Get CARD summary files
    st_dir <- file.path(CARD_DIR, st)
    if (!dir.exists(st_dir)) {
      cat(sprintf("  Warning: Directory not found: %s\n", st_dir))
      next
    }
    
    card_files <- list.files(st_dir, pattern = "_card\\.tsv$", full.names = TRUE)
    if (length(card_files) == 0) {
      cat(sprintf("  Warning: No CARD files found in %s\n", st_dir))
      next
    }
    
    cat(sprintf("  Found %d CARD files\n", length(card_files)))
    
    # Load metadata
    meta_file <- file.path(META_DIR, paste0(st, "_filtered.xlsx"))
    if (!file.exists(meta_file)) {
      cat(sprintf("  Warning: Metadata not found: %s\n", meta_file))
      next
    }
    
    meta <- read_excel(meta_file)
    cat(sprintf("  Loaded metadata: %d rows\n", nrow(meta)))
    
    # Standardize metadata column names
    meta <- meta %>%
      rename_with(~gsub(" ", "_", .x, fixed = TRUE)) %>%
      rename_with(~gsub("Collection_Year", "Collection_Year", .x, fixed = FALSE))
    
    # Get metadata names for matching
    meta_names <- if ("Name" %in% colnames(meta)) {
      meta$Name
    } else if ("Genome_ID" %in% colnames(meta)) {
      meta$Genome_ID
    } else {
      stop("Metadata must have 'Name' or 'Genome_ID' column")
    }
    
    # Extract genome IDs from CARD files
    genome_ids <- sapply(card_files, function(f) {
      extract_genome_id(f, st, "card")
    })
    
    # Perform batch matching
    matched_names <- batch_match(genome_ids, meta_names)
    
    # Calculate and store statistics
    stats <- get_matching_stats(matched_names, st)
    matching_stats[[st]] <- stats
    
    cat(sprintf("  Matched: %d / %d (%.1f%%)\n", 
                stats$Matched, stats$Total, stats$Match_Rate))
    
    # Load and combine CARD data
    card_data_list <- lapply(seq_along(card_files), function(i) {
      df <- tryCatch({
        read_tsv(card_files[i], show_col_types = FALSE)
      }, error = function(e) {
        cat(sprintf("  Error reading %s: %s\n", basename(card_files[i]), e$message))
        return(NULL)
      })
      
      if (is.null(df)) return(NULL)
      
      # Add genome identifier
      df$Genome_ID <- matched_names[i]
      df$Original_Filename <- basename(card_files[i])
      df
    })
    
    # Combine all CARD data for this ST
    card_data_list <- Filter(Negate(is.null), card_data_list)
    if (length(card_data_list) == 0) {
      cat(sprintf("  Warning: No valid CARD data for %s\n", st))
      next
    }
    
    st_card <- bind_rows(card_data_list)
    
    # Merge with metadata
    st_card <- st_card %>%
      left_join(meta, by = c("Genome_ID" = 
                              if ("Name" %in% colnames(meta)) "Name" else "Genome_ID")) %>%
      mutate(ST = st)
    
    all_card_data[[st]] <- st_card
    cat(sprintf("  Final dataset: %d genomes\n\n", nrow(st_card)))
  }
  
  # Combine all STs
  if (length(all_card_data) == 0) {
    stop("No data loaded for any sequence type")
  }
  
  combined_data <- bind_rows(all_card_data)
  matching_summary <- bind_rows(matching_stats)
  
  # Print summary
  print_matching_summary(matching_summary)
  
  # Return results
  list(
    data = combined_data,
    summary = matching_summary
  )
}
