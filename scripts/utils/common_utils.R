#!/usr/bin/env Rscript
################################################################################
# SHARED UTILITY FUNCTIONS
# Common functions used across all analysis modules
################################################################################

#' Modern publication-quality theme
#'
#' @param base_size Base font size
#' @return ggplot2 theme object
theme_modern <- function(base_size = 13) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title    = element_text(face = "bold", size = rel(1.3),
                                   margin = margin(0,0,10,0), hjust = 0),
      plot.subtitle = element_text(size = rel(0.9), color = "grey40",
                                   margin = margin(0,0,15,0), hjust = 0),
      plot.caption  = element_text(size = rel(0.75), color = "grey50",
                                   hjust = 1, margin = margin(15,0,0,0)),
      axis.title    = element_text(face = "bold", size = rel(0.95)),
      axis.text     = element_text(color = "grey30", size = rel(0.85)),
      axis.title.x  = element_text(margin = margin(10,0,0,0)),
      axis.title.y  = element_text(margin = margin(0,10,0,0)),
      axis.line     = element_line(color = "grey80", linewidth = 0.5),
      axis.ticks    = element_line(color = "grey80", linewidth = 0.4),
      panel.grid.major  = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.minor  = element_blank(),
      panel.background  = element_rect(fill = "white", color = NA),
      legend.position   = "right",
      legend.title      = element_text(face = "bold", size = rel(0.9)),
      legend.text       = element_text(size = rel(0.85)),
      legend.key        = element_blank(),
      legend.background = element_rect(fill = "white", color = NA),
      legend.spacing.y  = unit(0.3, "cm"),
      strip.text        = element_text(face = "bold", size = rel(0.95),
                                       margin = margin(5,5,5,5)),
      strip.background  = element_rect(fill = "grey95", color = NA),
      plot.margin       = margin(15,15,15,15),
      plot.background   = element_rect(fill = "white", color = NA)
    )
}

#' Publication theme (classic style)
#'
#' @param base_size Base font size
#' @return ggplot2 theme object
theme_pub <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title         = element_text(face = "bold", size = rel(1.15), hjust = 0),
      plot.subtitle      = element_text(color = "grey40", size = rel(0.88), hjust = 0,
                                        margin = margin(2, 0, 8, 0)),
      axis.title         = element_text(face = "bold", size = rel(0.95)),
      axis.text          = element_text(color = "grey20", size = rel(0.85)),
      axis.line          = element_line(color = "grey30", linewidth = 0.45),
      axis.ticks         = element_line(color = "grey30", linewidth = 0.35),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      strip.background   = element_rect(fill = "grey95", color = NA),
      strip.text         = element_text(face = "bold", size = rel(0.9)),
      legend.position    = "right",
      legend.title       = element_text(face = "bold", size = rel(0.88)),
      legend.text        = element_text(size = rel(0.82)),
      legend.key.size    = unit(0.38, "cm"),
      plot.margin        = margin(12, 12, 12, 12),
      plot.background    = element_rect(fill = "white", color = NA),
      panel.background   = element_rect(fill = "white", color = NA)
    )
}

#' Save plot to file
#'
#' @param plot ggplot object
#' @param name Output filename (without extension)
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param output_dir Output directory
#' @return Invisible NULL
save_plot <- function(plot, name, width = 12, height = 7, output_dir = "outputs/plots") {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(
    filename = file.path(output_dir, sprintf("%s.png", name)),
    plot = plot, width = width, height = height, dpi = 300, bg = "white"
  )
  cat(sprintf("  ✓ %s.png\n", name))
  invisible(NULL)
}

#' Convert p-value to significance stars
#'
#' @param p Numeric p-value
#' @return Character string with stars or "ns"
p_star <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

#' Safe odds ratio calculation using Fisher's exact test
#'
#' @param mat 2x2 contingency table
#' @return Tibble with OR, CI, and p-value
safe_or <- function(mat) {
  tryCatch({
    r <- epitools::oddsratio.fisher(mat)
    tibble::tibble(
      OR    = r$measure[2, "estimate"],
      CI_lo = r$measure[2, "lower"],
      CI_hi = r$measure[2, "upper"],
      p     = r$p.value[2, "fisher.exact"]
    )
  }, error = function(e) {
    tibble::tibble(OR = NA_real_, CI_lo = NA_real_, CI_hi = NA_real_, p = NA_real_)
  })
}

#' Detect directory from candidates
#'
#' @param candidates Character vector of potential directories
#' @param sentinel_glob Glob pattern to verify directory
#' @return First valid directory or first candidate
.detect_dir <- function(candidates, sentinel_glob) {
  for (d in candidates) {
    if (dir.exists(d) && length(Sys.glob(file.path(d, sentinel_glob))) > 0)
      return(d)
  }
  candidates[1]
}

#' Standard color palettes
ST_COLORS <- c(
  "ST10"  = "#FF6B6B", "ST131" = "#4ECDC4",
  "ST69"  = "#95E1D3", "ST73"  = "#F38181", "ST95"  = "#FFA07A"
)

ST_COLORS_ALT <- c(
  ST10 = "#E63946", ST131 = "#457B9D",
  ST69 = "#2A9D8F", ST73  = "#E9C46A", ST95 = "#F4A261"
)

HEATMAP_COLORS <- c("#440154", "#31688E", "#35B779", "#FDE724")

#' Extract genome ID from filename
#'
#' @param filename Character vector of filenames
#' @param st Sequence type (e.g., "ST10")
#' @param db Database type ("card" or "vfdb")
#' @return Cleaned genome IDs
extract_genome_id <- function(filename, st, db = "card") {
  prefix <- paste0(db, "/", st, "/Escherichia_coli_")
  suffix <- paste0("_", db, "\\.tsv$")
  
  cleaned <- sub(prefix, "", filename, fixed = TRUE)
  cleaned <- sub(suffix, "", cleaned)
  cleaned
}
