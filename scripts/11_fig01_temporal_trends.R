#!/usr/bin/env Rscript
#
# Combined figure: temporal trends (top) + Mann-Kendall summary (bottom)
# All databases (VFDB, CARD, VF, ResFinder) across 5 STs.
#
suppressPackageStartupMessages({
  library(tidyverse); library(readxl); library(ggplot2); library(patchwork); library(scales)
})
source("config.R")

st_list <- c("ST10", "ST69", "ST73", "ST95", "ST131")
OUT <- file.path(config$OUTPUT_DIR, "figures_temporal_mk_combined")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

mk_test <- function(x) {
  n <- length(x); if (n < 3) return(list(tau = NA, p = 1))
  s <- 0; for (i in 1:(n-1)) for (j in (i+1):n) s <- s + sign(x[j] - x[i])
  denom <- n * (n - 1) / 2; tau <- s / denom
  var_s <- n * (n - 1) * (2 * n + 5) / 18
  z <- if (var_s > 0) (s - sign(s)) / sqrt(var_s) else 0
  list(tau = tau, p = 2 * pnorm(-abs(z)))
}

make_genome_id <- function(x) {
  x <- as.character(x); x <- trimws(x); x <- sub("[.]0$", "", x)
  case_when(grepl("^Escherichia_coli_", x) ~ x,
            grepl("^E\\.coli_", x) ~ sub("^E\\.coli_", "Escherichia_coli_", x),
            TRUE ~ paste0("Escherichia_coli_", x))
}

cat("Loading binary matrices...\n")
vf_binary <- read.delim(config$VF_BINARY, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
resf_binary <- read.delim(file.path(config$INPUT_DIR, "resfinder_summary", "resfinder_binary_matrix.tsv"),
                           header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

all_data <- list()
for (st in st_list) {
  cat("  ", st, "...\n", sep = "")
  meta_file <- config$st_metadata(st)
  if (!file.exists(meta_file)) { cat("    No metadata\n"); next }
  meta <- read_xlsx(meta_file)
  name_col <- grep("^(Name|genome|strain|isolate|assembly|sample)$", colnames(meta), value = TRUE)[1]
  if (is.na(name_col)) name_col <- colnames(meta)[2]
  year_col <- grep("year", colnames(meta), ignore.case = TRUE, value = TRUE)[1]
  meta <- meta %>% rename(genome_id = all_of(name_col)) %>%
    mutate(genome_id = make_genome_id(genome_id), year = as.integer(.data[[year_col]])) %>%
    filter(!is.na(year)) %>% select(genome_id, year)

  vfdb_file <- config$st_vfdb_summary(st)
  if (file.exists(vfdb_file)) {
    v <- read_tsv(vfdb_file, show_col_types = FALSE, col_types = cols(.default = "c")) %>%
      mutate(genome_id = make_genome_id(str_remove(basename(.data[["#FILE"]]), "_vfdb\\.tsv$")),
             burden = as.numeric(NUM_FOUND)) %>%
      select(genome_id, burden) %>% inner_join(meta, by = "genome_id") %>%
      group_by(year) %>% summarise(m = mean(burden, na.rm = TRUE), se = sd(burden, na.rm = TRUE) / sqrt(n()), n = n(), .groups = "drop") %>%
      filter(n >= 5) %>% mutate(st = st, db = "VFDB")
    all_data[[length(all_data) + 1]] <- v
  }

  card_file <- config$st_card_burden(st)
  if (file.exists(card_file)) {
    c <- read_tsv(card_file, show_col_types = FALSE, col_types = cols(.default = "c")) %>%
      mutate(genome_id = make_genome_id(str_remove(basename(.data[["#FILE"]]), "_card\\.tsv$")),
             burden = as.numeric(NUM_FOUND)) %>%
      select(genome_id, burden) %>% inner_join(meta, by = "genome_id") %>%
      group_by(year) %>% summarise(m = mean(burden, na.rm = TRUE), se = sd(burden, na.rm = TRUE) / sqrt(n()), n = n(), .groups = "drop") %>%
      filter(n >= 5) %>% mutate(st = st, db = "CARD")
    all_data[[length(all_data) + 1]] <- c
  }

  vf <- vf_binary %>% filter(st == !!st) %>%
    mutate(burden = rowSums(across(-c(st, genome), ~ as.integer(!is.na(.x) & .x != "" & .x != "0"))),
           genome_id = make_genome_id(genome)) %>%
    select(genome_id, burden) %>% distinct(genome_id, .keep_all = TRUE) %>%
    inner_join(meta, by = "genome_id") %>%
    group_by(year) %>% summarise(m = mean(burden, na.rm = TRUE), se = sd(burden, na.rm = TRUE) / sqrt(n()), n = n(), .groups = "drop") %>%
    filter(n >= 5) %>% mutate(st = st, db = "VF")
  all_data[[length(all_data) + 1]] <- vf

  resf <- resf_binary %>% filter(st == !!st) %>%
    mutate(burden = rowSums(across(-c(st, genome), ~ as.integer(!is.na(.x) & .x != "" & .x != "0"))),
           genome_id = make_genome_id(genome)) %>%
    select(genome_id, burden) %>% distinct(genome_id, .keep_all = TRUE) %>%
    inner_join(meta, by = "genome_id") %>%
    group_by(year) %>% summarise(m = mean(burden, na.rm = TRUE), se = sd(burden, na.rm = TRUE) / sqrt(n()), n = n(), .groups = "drop") %>%
    filter(n >= 5) %>% mutate(st = st, db = "ResFinder")
  all_data[[length(all_data) + 1]] <- resf
}

df <- bind_rows(all_data)

# --- Top panel: Temporal trends ---
st_colors <- c(ST10 = "#FF6B6B", ST131 = "#4ECDC4", ST69 = "#95E1D3", ST73 = "#F38181", ST95 = "#FFA07A")
st_levels <- c("ST10", "ST69", "ST73", "ST95", "ST131")
db_levels <- c("CARD", "ResFinder", "VFDB", "VF")
db_names <- c(CARD = "CARD (ARG)", ResFinder = "ResFinder (ARG)", VFDB = "VFDB (VF)", VF = "VirulenceFinder (VF)")

df <- df %>% mutate(st = factor(st, levels = st_levels), db = factor(db, levels = db_levels))

p_temporal <- df %>% ggplot(aes(x = year, y = m, color = st, fill = st)) +
  geom_ribbon(aes(ymin = m - se, ymax = m + se), alpha = 0.1, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(aes(size = n), alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  facet_wrap(~ db, nrow = 1, ncol = 4, labeller = labeller(db = db_names), scales = "free_y") +
  scale_color_manual(values = st_colors) +
  scale_fill_manual(values = st_colors) +
  scale_size_continuous(range = c(1, 4), name = "n") +
  scale_x_continuous(breaks = 2016:2025) +
  labs(x = "Collection Year", y = "Mean burden", color = "ST", fill = "ST") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

# --- Bottom panel: MK results ---
mk <- df %>% group_by(st, db) %>% arrange(year) %>%
  summarise(mk = list(mk_test(m)), .groups = "drop") %>%
  mutate(tau = sapply(mk, `[[`, "tau"), p = sapply(mk, `[[`, "p"),
         sig = p < 0.05,
         sig_label = ifelse(p < 0.05, ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", "*")), "ns"))

mk_plot <- mk %>%
  mutate(st = factor(st, levels = st_levels), db = factor(db, levels = db_levels),
         tau_label = sprintf("%+.2f", tau),
         cell_label = paste0(tau_label, sig_label))

# MK lollipop row
p_mk_lolli <- mk_plot %>%
  ggplot(aes(x = tau, y = st, color = st, alpha = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.7) +
  geom_segment(aes(xend = 0, yend = st), linewidth = 1.0) +
  geom_point(aes(size = sig, shape = tau > 0), stroke = 1.2) +
  facet_wrap(~ db, nrow = 1) +
  scale_color_manual(values = st_colors, guide = "none") +
  scale_alpha_manual(values = c("TRUE" = 1.0, "FALSE" = 0.40), guide = "none") +
  scale_size_manual(values = c("TRUE" = 5, "FALSE" = 3.5), guide = "none") +
  scale_shape_manual(values = c("TRUE" = 24, "FALSE" = 25), guide = "none") +
  labs(x = "Kendall's tau", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(face = "bold", size = 11),
        strip.text = element_text(face = "bold", size = 11))

# MK heatmap row
p_mk_heat <- mk_plot %>%
  ggplot(aes(x = db, y = st, fill = tau)) +
  geom_tile(color = "white", linewidth = 2.0, width = 0.95, height = 0.95) +
  geom_text(aes(label = cell_label, color = abs(tau) > 0.5), size = 4.5, fontface = "bold") +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "grey15"), guide = "none") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#D73027", midpoint = 0,
                       limits = c(-1, 1), name = "Kendall's tau") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(face = "bold", size = 11),
        axis.text.x = element_text(size = 11, face = "bold"),
        legend.position = "right")

fig <- (p_temporal / p_mk_lolli / p_mk_heat) + plot_layout(heights = c(2.5, 1.2, 1.2)) +
  plot_annotation(theme = theme(plot.margin = margin(5, 5, 5, 5)))

ggsave(file.path(OUT, "Fig_temporal_MK_combined.png"), fig, width = 10, height = 8, dpi = 300, bg = "white")
ggsave(file.path(OUT, "Fig_temporal_MK_combined.pdf"), fig, width = 10, height = 8, bg = "white")
cat("Saved:", file.path(OUT, "Fig_temporal_MK_combined.png"), "\n")
