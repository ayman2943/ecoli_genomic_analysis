#!/usr/bin/env Rscript
# Per-cluster temporal trends + MK analysis: which clusters acquire/lose genes?

suppressPackageStartupMessages({
  library(tidyverse); library(patchwork); library(ggrepel)
})
source("config.R")
ts <- function(m) cat("[", Sys.time(), "] ", m, "\n", sep = "")

EARLY_YEARS <- c(2016, 2017, 2018)
LATE_YEARS <- c(2022, 2023, 2024, 2025)
OUT <- file.path(config$OUTPUT_DIR, "tables", "cluster_temporal_trends")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# Mann-Kendall test
mk_test <- function(x) {
  n <- length(x); if (n < 3) return(list(tau = NA, p = 1))
  s <- 0; for (i in 1:(n-1)) for (j in (i+1):n) s <- s + sign(x[j] - x[i])
  denom <- n * (n - 1) / 2; tau <- s / denom
  var_s <- n * (n - 1) * (2 * n + 5) / 18
  z <- if (var_s > 0) (s - sign(s)) / sqrt(var_s) else 0
  list(tau = tau, p = 2 * pnorm(-abs(z)))
}

cluster_colors <- function(n) {
  if (n <= 9) RColorBrewer::brewer.pal(max(n, 3), "Set1") else rainbow(n)
}

# Full analysis: table + yearly data for plotting
analyze <- function(label, df, total_col = NULL) {
  ts(paste("Analyzing:", label))
  df <- df %>% filter(!is.na(year), !is.na(shell_cluster))
  if (is.null(total_col)) {
    gene_cols <- setdiff(colnames(df), c("genome_id", "shell_cluster", "year",
      "clinical_binary", "source_niche", "country", "st", "genome",
      "niche", "continent", "country", "total_vf", "raw_cluster.x",
      "raw_cluster.y", "raw_cluster", "X.FILE", "X", "cluster"))
    gene_cols <- gene_cols[sapply(df[, gene_cols], is.numeric)]
    if (length(gene_cols) == 0) { ts("  No numeric gene columns"); return(NULL) }
    df <- df %>% mutate(total = rowSums(across(all_of(gene_cols)), na.rm = TRUE))
  } else {
    df <- df %>% mutate(total = as.numeric(.data[[total_col]]))
  }

  # Yearly mean per cluster
  yearly <- df %>%
    group_by(shell_cluster, year) %>%
    summarise(m = mean(total, na.rm = TRUE), n = n(), se = sd(total, na.rm = TRUE) / sqrt(n()),
              .groups = "drop")

  # Summary stats per cluster
  res <- df %>%
    group_by(shell_cluster) %>%
    summarise(
      n = n(),
      n_early = sum(year %in% EARLY_YEARS),
      n_late = sum(year %in% LATE_YEARS),
      mean_total = mean(total, na.rm = TRUE),
      mean_early = mean(total[year %in% EARLY_YEARS], na.rm = TRUE),
      mean_late = mean(total[year %in% LATE_YEARS], na.rm = TRUE),
      delta_mean = mean_late - mean_early,
      .groups = "drop"
    )

  # Linear regression per cluster
  clusts <- unique(df$shell_cluster)
  fit_lm <- function(cl) {
    d <- df %>% filter(shell_cluster == cl)
    if (n_distinct(d$year) < 3) return(tibble(shell_cluster = cl, slope = NA_real_, lm_p = 1))
    m <- tryCatch(lm(total ~ year, data = d), error = function(e) NULL)
    if (is.null(m)) return(tibble(shell_cluster = cl, slope = NA_real_, lm_p = 1))
    s <- summary(m)$coefficients
    tibble(shell_cluster = cl, slope = s["year", 1], lm_p = s["year", 4])
  }
  lm_res <- map_dfr(clusts, fit_lm)

  # MK per cluster (on yearly means)
  mk_res <- yearly %>%
    group_by(shell_cluster) %>%
    arrange(year) %>%
    summarise(mk = list(mk_test(m)), .groups = "drop") %>%
    mutate(tau = sapply(mk, `[[`, "tau"), mk_p = sapply(mk, `[[`, "p"))

  res <- res %>%
    left_join(lm_res, by = "shell_cluster") %>%
    left_join(mk_res %>% select(shell_cluster, tau, mk_p), by = "shell_cluster") %>%
    mutate(
      lm_p_adj = p.adjust(lm_p, method = "BH"),
      mk_p_adj = p.adjust(mk_p, method = "BH"),
      direction = case_when(
        slope > 0 & lm_p_adj < 0.05 ~ "Acquiring",
        slope < 0 & lm_p_adj < 0.05 ~ "Losing",
        TRUE ~ "Stable"
      ),
      lm_sig = case_when(
        lm_p_adj < 0.001 ~ "***", lm_p_adj < 0.01 ~ "**",
        lm_p_adj < 0.05 ~ "*", TRUE ~ "ns"
      ),
      mk_sig = case_when(
        mk_p_adj < 0.001 ~ "***", mk_p_adj < 0.01 ~ "**",
        mk_p_adj < 0.05 ~ "*", TRUE ~ "ns"
      )
    ) %>%
    arrange(lm_p, desc(abs(slope)))

  outfile <- file.path(OUT, paste0(gsub(" ", "_", label), "_cluster_temporal_trends.csv"))
  write_csv(res, outfile)
  ts(paste("  Saved:", outfile))

  list(table = res, yearly = yearly)
}

# Plot function: temporal lines + MK lollipop for one case
plot_case <- function(data, label, cols) {
  yearly <- data$yearly
  res <- data$table
  clusts <- unique(yearly$shell_cluster)
  clust_cols <- cluster_colors(length(clusts))
  names(clust_cols) <- sort(clusts)

  # Temporal lines (visual only — MK panel below carries significance)
  p_temp <- yearly %>% ggplot(aes(x = year, y = m, color = shell_cluster)) +
    geom_ribbon(aes(ymin = m - se, ymax = m + se, fill = shell_cluster), alpha = 0.1, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(aes(size = n), alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", linewidth = 0.4, alpha = 0.5) +
    scale_color_manual(values = clust_cols) +
    scale_fill_manual(values = clust_cols) +
    scale_size_continuous(range = c(1, 3.5), name = "n") +
    scale_x_continuous(breaks = 2016:2025) +
    labs(x = NULL, y = "Mean burden", color = "Cluster", fill = "Cluster",
         title = label) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 14),
          plot.margin = margin(5, 15, 5, 5),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 11))

  # MK lollipop (raw p used for both label and significance — consistent)
  mk_plot <- res %>%
    mutate(shell_cluster = factor(shell_cluster, levels = rev(sort(unique(shell_cluster)))),
           mk_sig_raw = case_when(mk_p < 0.001 ~ "***", mk_p < 0.01 ~ "**",
                                  mk_p < 0.05 ~ "*", TRUE ~ "ns"),
           mk_p_label = ifelse(mk_p < 0.001, sprintf("p < 0.001 %s", mk_sig_raw),
                               sprintf("p = %.3f %s", mk_p, mk_sig_raw)),
           tau_label = sprintf("tau = %.2f, %s", tau, mk_p_label),
           label_x = ifelse(tau >= 0, tau + 0.15, tau - 0.15)) %>%
    ggplot(aes(x = tau, y = shell_cluster, color = shell_cluster)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.7) +
    geom_segment(aes(xend = 0, yend = shell_cluster), linewidth = 0.8, alpha = 0.4) +
    geom_point(aes(size = mk_p < 0.05), stroke = 1.0) +
    geom_label_repel(aes(label = tau_label, x = label_x),
              size = 3.5, color = "grey20", fill = "white", label.size = 0,
              segment.color = NA, max.overlaps = 30,
              direction = "y", box.padding = 0.4, segment.padding = 0, fontface = "plain") +
    scale_color_manual(values = clust_cols, guide = "none") +
    scale_size_manual(values = c("TRUE" = 4, "FALSE" = 2.5), guide = "none") +
    labs(x = "Kendall tau", y = NULL) +
    xlim(c(min(res$tau, -0.1) - 0.6, max(res$tau, 0.1) + 0.7)) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = 11),
          axis.text.x = element_text(size = 10))

  p_temp / mk_plot + plot_layout(heights = c(1.2, 1))
}

# === ST69 VFDB ===
ts("Loading ST69 VFDB")
st69 <- read_csv("output/ST69/vfdb_analysis/04_master_shell_cluster_metadata_VFDB_table.csv",
  show_col_types = FALSE)
d1 <- analyze("ST69_VFDB", st69, total_col = "total_vf")

# === ST10 VirulenceFinder ===
ts("Loading ST10 VF")
st10_vf <- read_csv("output/ST10/virulencefinder_validation/04_master_shell_cluster_metadata_VF_table.csv",
  show_col_types = FALSE)
d2 <- analyze("ST10_VF", st10_vf, total_col = "total_vf")

# === ST10 ResFinder ===
ts("Loading ST10 ResFinder")
rf_raw <- read_tsv("finder_result/resfinder_summary/resfinder_binary_matrix.tsv",
  show_col_types = FALSE) %>%
  filter(st == "ST10") %>%
  mutate(genome_id = as.character(genome))
st10_vfdb <- read_csv("output/ST10/vfdb_analysis/04_master_shell_cluster_metadata_VFDB_table.csv",
  show_col_types = FALSE) %>% select(genome_id, shell_cluster, year)
rf_df <- rf_raw %>% inner_join(st10_vfdb, by = "genome_id")
d3 <- analyze("ST10_ResFinder", rf_df)

# === Figure: 3 panels ===
ts("Generating figure")
p1 <- plot_case(d1, "ST69 VFDB")
p2 <- plot_case(d2, "ST10 VirulenceFinder")
p3 <- plot_case(d3, "ST10 ResFinder")

fig <- wrap_plots(p1, p2, p3, ncol = 1)
out_png <- file.path(OUT, "cluster_temporal_trends_3panel.png")
out_pdf <- file.path(OUT, "cluster_temporal_trends_3panel.pdf")
ggsave(out_png, fig, width = 10, height = 16, dpi = 300, bg = "white")
ggsave(out_pdf, fig, width = 10, height = 16, bg = "white")
ts(paste("Saved figure:", out_png))
ts("DONE")
