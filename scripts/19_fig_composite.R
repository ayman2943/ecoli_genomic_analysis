#!/usr/bin/env Rscript
# Generate publication-ready figures for Word (one page each, 7.5x10 inches at 300 DPI)
# Uses png + grid for image composition
suppressPackageStartupMessages({
  library(png); library(grid); library(gridExtra)
})

IN <- "figures/new"
OUT <- "figures/new"
DPI <- 300
W <- 7.5  # inches width
H <- 10   # inches height

read_img <- function(f) {
  png::readPNG(file.path(IN, f))
}

# ======================================================================
# Fig 2+3 combined: silhouette stacked above content
# ======================================================================
cat("Fig 2+3 combined...\n")
img2 <- read_img("Fig_02_silhouette_3panel.png")
img3 <- read_img("Fig_03_content_by_cluster_3panel.png")

g2 <- rasterGrob(img2, interpolate = TRUE)
g3 <- rasterGrob(img3, interpolate = TRUE)

# Stack vertically with labels
lay <- rbind(c(1), c(2))
p <- grid.arrange(g2, g3, nrow = 2,
                   top = textGrob("Figure 2-3: Cluster silhouette and gene content", 
                                   gp = gpar(fontface = "bold", fontsize = 14)))
png(file.path(OUT, "Fig_02_03_combined.png"), width = W, height = H, units = "in", res = DPI)
grid.draw(p)
dev.off()
cat("  Saved Fig_02_03_combined.png\n")

# ======================================================================
# Single-figure re-saves (just resize to page, no recomposition needed)
# ======================================================================
resize_and_save <- function(in_name, out_name, label) {
  cat(label, "...\n")
  img <- read_img(in_name)
  g <- rasterGrob(img, interpolate = TRUE)
  png(file.path(OUT, out_name), width = W, height = H, units = "in", res = DPI)
  grid.draw(g)
  dev.off()
  cat("  Saved", out_name, "\n")
}

resize_and_save("Fig_01_temporal_MK_combined.png", "Fig_01_one_page.png", "Fig 1")
resize_and_save("Fig_05_decomposition_3panel_nosmooth.png", "Fig_05_one_page.png", "Fig 5")
resize_and_save("Fig_06_cluster_temporal_trends.png", "Fig_06_one_page.png", "Fig 6")
resize_and_save("Fig_07_combined_summary 2.png", "Fig_07_one_page.png", "Fig 7")
resize_and_save("Fig_8_sensitivity.png", "Fig_08_one_page.png", "Fig 8")

cat("Done.\n")
