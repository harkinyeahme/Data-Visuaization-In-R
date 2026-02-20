# ══════════════════════════════════════════════════════════════════════════════
# MULTI-PANEL FIGURE GENERATION SCRIPT
# Description: This script generates two multi-panel figures (Figure 1 and
#              Figure 3) from biological datasets. It combines pheatmap and
#              base R plots into publication-ready panels using grid graphics.
# Author:      AKINYEMI OLANREWAJU AKINWALE
# Date:        [20-02-2026]
# ══════════════════════════════════════════════════════════════════════════════


# ── SECTION 1: Install and Load Required Libraries ────────────────────────────
# These packages are required for plotting, layout, and data import.
# Run install.packages() only once; comment it out after first use.

# install.packages(c("gridExtra",    # arranging multiple plots in a grid
#                   "gridGraphics", # converting base R plots to grid objects
#                   "pheatmap",     # heatmap generation
#                   "png",          # reading PNG files into R
#                   "readxl",       # reading Excel files
#                   "igraph",       # network/graph plotting
#                   "RColorBrewer"  # color palettes for heatmaps
#))

library(grid)         # core grid graphics system
library(gridExtra)    # arrangeGrob and grid.arrange for panel layout
library(gridGraphics) # allows base R plots to be used in grid layout
library(pheatmap)     # heatmap with clustering
library(png)          # read PNG images
library(readxl)       # read Excel sheets
library(RColorBrewer) # color palettes
library(igraph)       # network graphs


# ── SECTION 2: Helper Functions ───────────────────────────────────────────────

# capture_base_plot():
# Base R plots cannot be directly used in grid.arrange() because they use
# a different graphics system. This function:
#   1. Renders the base R plot to a temporary PNG file at high resolution
#   2. Reads the PNG back into R as a raster image
#   3. Converts it to a grid-compatible rasterGrob object
# This makes base R plots behave like grid objects for panel assembly.

capture_base_plot <- function(expr) {
  tf <- tempfile(fileext = ".png")          # create a temporary file path
  png(tf, width = 1600, height = 1200, res = 150) # open PNG device at high res
  eval(expr)                                # evaluate and draw the plot
  dev.off()                                 # close the PNG device
  img <- png::readPNG(tf)                   # read the saved PNG back into R
  grid::rasterGrob(img,                     # convert to grid raster object
                   width  = unit(1, "npc"),
                   height = unit(1, "npc"))
}

# label_plot():
# Adds a bold label/title to the top-left of any grob (grid object).
# Works for both pheatmap gtables and rasterGrobs from capture_base_plot().
# Uses arrangeGrob() which is more compatible than grobTree() for gtables.

label_plot <- function(grob, label) {
  arrangeGrob(
    grob,
    top = textGrob(label,
                   x    = 0.02, just = "left",   # left-aligned position
                   gp   = gpar(fontsize = 14, fontface = "bold"))
  )
}


# ══════════════════════════════════════════════════════════════════════════════
#  FIGURE 1: Cancer and Gene Expression Analysis
#  Panels: Expression Heatmap | Volcano Plot | Scatter Plots | Density Plot
# ══════════════════════════════════════════════════════════════════════════════

# ── Plot 1A: Gene Expression Heatmap ─────────────────────────────────────────
# Dataset: Normalized counts from HBR/UHR RNA-seq experiment
# Goal: Visualize expression patterns of top differentially expressed genes
# Method: Hierarchical clustering on both rows (genes) and columns (samples)

link <- "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_top_deg_normalized_counts.csv"

expr_data <- read.delim(link, sep = ",", header = TRUE)

# Move gene names from first column into row names (required for matrix)
rownames(expr_data) <- expr_data[, 1]
expr_data <- expr_data[, -1]

# Convert data frame to numeric matrix (pheatmap requires a matrix)
mat <- as.matrix(expr_data)

# Generate heatmap; silent = TRUE suppresses auto-display so we can store it
# $gtable extracts the grid table object for use in panel assembly
pA <- pheatmap(mat,
               cluster_rows  = TRUE,
               cluster_cols  = TRUE,
               border_color  = "black",
               color         = colorRampPalette(c("white", "blue"))(100),
               fontsize_row  = 8,
               silent        = TRUE)$gtable


# ── Plot 1B: Volcano Plot ─────────────────────────────────────────────────────
# Dataset: DEG results with log2FoldChange and adjusted p-values
# Goal: Visualize significantly up/down regulated genes
# Dashed lines mark fold-change threshold (±1) and significance (p = 0.05)

expr_data2 <- read.delim(
  "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_deg_chr22_with_significance.csv",
  sep = ",", header = TRUE)

# Color map: green = upregulated, orange = downregulated, grey = not significant
vol_colors <- c("up" = "green", "down" = "orange", "ns" = "grey")

pB <- capture_base_plot(quote({
  plot(expr_data2$log2FoldChange, expr_data2$X.log10PAdj,
       col      = vol_colors[expr_data2$significance],
       pch      = 19,
       ylab     = "-log10PAdj",
       xlab     = "log2FoldChange",
       cex      = 1.2,
       cex.axis = 1.5,
       cex.lab  = 1.5)
  legend("topright", pch = 19,
         legend = c("down", "ns", "up"),
         col    = vol_colors[c("down", "ns", "up")],
         title  = "Significance", cex = 1.3)
  abline(v = c(-1, 1), lty = 2)           # fold change cutoff lines
  abline(h = -log10(0.05), lty = 2)       # significance threshold line
}))


# ── Plot 1C: Radius Mean vs Texture Mean Scatter ──────────────────────────────
# Dataset: Breast cancer diagnostic measurements
# Goal: Explore relationship between tumor radius and texture by diagnosis
# Color: Blue = Malignant (M), Orange = Benign (B)

breast_cancer_data <- read.delim(
  "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv",
  sep = ",", header = TRUE)

bc_cols <- c("B" = "orange", "M" = "blue")

pC <- capture_base_plot(quote({
  plot(breast_cancer_data$radius_mean, breast_cancer_data$texture_mean,
       col      = bc_cols[breast_cancer_data$diagnosis],
       pch      = 19,
       ylab     = "texture_mean",
       xlab     = "radius_mean",
       cex      = 1.2,
       cex.axis = 1.5,
       cex.lab  = 1.5)
  legend("topright", pch = 19,
         legend = c("M", "B"),
         col    = c("blue", "orange"),
         title  = "Diagnosis", cex = 1.3)
}))


# ── Plot 1D: Correlation Heatmap ──────────────────────────────────────────────
# Goal: Show pairwise Pearson correlations between 6 breast cancer features
# display_numbers = TRUE overlays the correlation values on each cell

features <- breast_cancer_data[, c("radius_mean", "texture_mean",
                                   "perimeter_mean", "area_mean",
                                   "smoothness_mean", "compactness_mean")]

# Compute pairwise correlation matrix
cor_mat <- cor(features, use = "complete.obs")

pD <- pheatmap(cor_mat,
               display_numbers = TRUE,
               number_format   = "%.1f",
               color           = colorRampPalette(c("white", "lightblue", "blue"))(100),
               border_color    = "black",
               cluster_rows    = FALSE,
               cluster_cols    = FALSE,
               legend_breaks   = seq(-1, 1, by = 0.2),
               legend_labels   = seq(-1, 1, by = 0.2),
               silent          = TRUE)$gtable


# ── Plot 1E: Smoothness vs Compactness Scatter ────────────────────────────────
# Goal: Compare smoothness and compactness between malignant and benign tumors
# xaxt = "n" suppresses default x-axis so we can draw a custom one with axis()

pE <- capture_base_plot(quote({
  plot(breast_cancer_data$smoothness_mean, breast_cancer_data$compactness_mean,
       xlim     = c(0.05, 0.17),
       ylim     = c(0.01, 0.35),
       col      = bc_cols[breast_cancer_data$diagnosis],
       pch      = 19,
       ylab     = "compactness_mean",
       xlab     = "smoothness_mean",
       xaxt     = "n",              # suppress default x-axis
       cex      = 1.2,
       cex.axis = 1.5,
       cex.lab  = 1.5)
  axis(side = 1, seq(0.05, 0.15, 0.025), cex.axis = 1.5)  # custom x-axis
  legend("topleft", pch = 19,
         legend = c("M", "B"),
         col    = c("blue", "orange"),
         title  = "Diagnosis", cex = 1.3)
}))


# ── Plot 1F: Density Plot of Area Mean ────────────────────────────────────────
# Goal: Compare distribution of area_mean between malignant and benign tumors
# polygon() fills the area under the density curve with transparent color

# Compute density estimates for each diagnosis group
dens_M <- density(breast_cancer_data$area_mean[breast_cancer_data$diagnosis == "M"])
dens_B <- density(breast_cancer_data$area_mean[breast_cancer_data$diagnosis == "B"])

pF <- capture_base_plot(quote({
  plot(dens_M,
       col      = "blue",
       lwd      = 2,
       xlab     = "Area Mean",
       ylab     = "Density",
       ylim     = c(0, 0.0035),
       main     = "",
       cex.axis = 1.5,
       cex.lab  = 1.5)
  polygon(dens_M, col = rgb(0, 0, 1,   0.3), border = NA)  # blue fill, M
  polygon(dens_B, col = rgb(1, 0.5, 0, 0.3), border = NA)  # orange fill, B
  lines(dens_M, col = "blue",   lwd = 2)   # redraw lines on top of polygon
  lines(dens_B, col = "orange", lwd = 2)
  legend("topright", title = "Diagnosis",
         legend = c("M", "B"),
         fill   = c(rgb(0, 0, 1, 0.3), rgb(1, 0.5, 0, 0.3)),
         cex    = 1.3)
}))


# ── Label Figure 1 Plots ──────────────────────────────────────────────────────
# label_plot() bakes the title into each grob so it appears in both
# individual saves and the combined panel without extra steps

pA_l <- label_plot(pA, "a) Expression Heatmap")
pB_l <- label_plot(pB, "b) Volcano Plot")
pC_l <- label_plot(pC, "c) Radius vs Texture")
pD_l <- label_plot(pD, "d) Correlation Heatmap")
pE_l <- label_plot(pE, "e) Smoothness vs Compactness")
pF_l <- label_plot(pF, "f) Area Mean Density")


# ── Save Individual Figure 1 Plots ───────────────────────────────────────────
# Each plot is saved separately at 300 dpi for publication quality.
# grid.newpage() clears the canvas before drawing to avoid overlap.

png("figure1_A_expression_heatmap.png", width = 2400, height = 1800, res = 300)
grid.newpage(); grid.draw(pA_l); dev.off()

png("figure1_B_volcano.png", width = 2400, height = 1800, res = 300)
grid.newpage(); grid.draw(pB_l); dev.off()

png("figure1_C_radius_texture.png", width = 2400, height = 1800, res = 300)
grid.newpage(); grid.draw(pC_l); dev.off()

png("figure1_D_correlation_heatmap.png", width = 2400, height = 1800, res = 300)
grid.newpage(); grid.draw(pD_l); dev.off()

png("figure1_E_smoothness_compactness.png", width = 2400, height = 1800, res = 300)
grid.newpage(); grid.draw(pE_l); dev.off()

png("figure1_F_density.png", width = 2400, height = 1800, res = 300)
grid.newpage(); grid.draw(pF_l); dev.off()


# ── Save Figure 1 Panel ───────────────────────────────────────────────────────
# Arrange all 6 plots in a 3-column grid panel and save as one PNG

png("figure1_panel.png", width = 3600, height = 2400, res = 300)
grid.arrange(pA_l, pB_l, pC_l,
             pD_l, pE_l, pF_l,
             ncol = 3)
dev.off()


# ══════════════════════════════════════════════════════════════════════════════
#  FIGURE 3: Immune Cell and Gene Regulation Analysis
#  Panels: Boxplot | Quadrant Scatter | Heatmaps | Bubble | Barplot | Network
# ══════════════════════════════════════════════════════════════════════════════

# Download the Excel workbook containing all Figure 3 datasets
# Each sheet corresponds to a different sub-analysis

url  <- "https://github.com/HackBio-Internship/2025_project_collection/raw/refs/heads/main/hb_stage_2.xlsx"
dest <- tempfile(fileext = ".xlsx")    # save to a temporary file
download.file(url, dest, mode = "wb") # mode "wb" = write binary (needed for Excel)


# ── Plot 3A: Cell Type Ratio Boxplot ─────────────────────────────────────────
# Dataset: Sheet "a" — cell type ratios across samples
# Goal: Compare ratio distributions across immune cell types
# las = 2 rotates x-axis labels vertically for readability

cell_type_ratio_data <- read_excel(dest, sheet = "a")
cell_type_ratio_data$cell_type <- as.factor(cell_type_ratio_data$cell_type)

# Named color vector for each cell type
box_cols <- c(B          = "darkblue",
              Basophil   = "lightgreen",
              DC         = "brown",
              HSPC       = "pink",
              Macrophage = "lightblue",
              Monocyte   = "darkgreen",
              Neutrophil = "purple",
              NK         = "orange",
              "Pre B"    = "darkorange",
              T          = "yellow")

p3A <- capture_base_plot(quote({
  boxplot(new_ratio ~ cell_type,
          data     = cell_type_ratio_data,
          ylab     = "Ratio",
          xlab     = "",
          las      = 2,
          col      = box_cols[levels(cell_type_ratio_data$cell_type)],
          cex.axis = 1.3,
          cex.lab  = 1.5)
}))


# ── Plot 3B: Half Life vs Alpha Quadrant Scatter ──────────────────────────────
# Dataset: Sheet "b" — gene half-life and alpha values
# Goal: Segment genes into 4 quadrants based on log2 thresholds
# Ccr2 labels Q1 (top-left), Camp labels Q4 (bottom-right)

alpha_half_life_data <- read_excel(dest, sheet = "b")

# Log2 transform both axes
x <- log2(alpha_half_life_data$half_life)
y <- log2(alpha_half_life_data$alpha)

# Assign each point to a quadrant based on threshold values
quad <- ifelse(x <= 2.5 & y >= -3.6, "Q1",
               ifelse(x >  2.5 & y >= -3.6, "Q2",
                      ifelse(x <= 2.5 & y <  -3.6, "Q3", "Q4")))

quad_cols <- c(Q1 = "green", Q2 = "red", Q3 = "darkgrey", Q4 = "blue")

p3B <- capture_base_plot(quote({
  plot(x, y,
       ylim     = c(-10, 0),
       col      = quad_cols[quad],
       pch      = 19,
       cex      = 0.8,
       xlab     = "log2(Half Life)",
       ylab     = "log2(Alpha)",
       cex.axis = 1.5,
       cex.lab  = 1.5)
  abline(v = 2.5,  lty = 2)   # vertical threshold
  abline(h = -3.6, lty = 2)   # horizontal threshold
  # Label Q1 (top-left quadrant) with gene name Ccr2
  text((min(x) + 2.5) / 2, (-3.6 + 0) / 2,
       "Ccr2", cex = 1.8, font = 2, col = "black")
  # Label Q4 (bottom-right quadrant) with gene name Camp
  text((2.5 + max(x)) / 2, (-10 + -3.6) / 2,
       "Camp", cex = 1.8, font = 2, col = "black")
}))


# ── Plot 3C: Time-series Expression Heatmap ───────────────────────────────────
# Dataset: Sheet "c" — gene expression across cell types and time points
# Goal: Visualize expression dynamics over time without column clustering
# Annotation bar added at top to show Time and CellType groupings

data_3C     <- read_excel(dest, sheet = "c")

# Convert tibble slice to numeric matrix (pheatmap requires numeric matrix)
data_3C_mat <- as.matrix(data_3C[, 2:50])
class(data_3C_mat) <- "numeric"

# Build annotation data frame from column names
# Assumes last 3 characters = time point, rest = cell type
annotation_col <- data.frame(
  Time     = substr(colnames(data_3C_mat),
                    nchar(colnames(data_3C_mat)) - 2,
                    nchar(colnames(data_3C_mat))),
  CellType = substr(colnames(data_3C_mat),
                    1,
                    nchar(colnames(data_3C_mat)) - 4),
  row.names = colnames(data_3C_mat)
)

p3C <- pheatmap(data_3C_mat,
                cluster_cols   = FALSE,         # preserve time order
                color          = colorRampPalette(brewer.pal(9, "Blues"))(200),
                show_colnames  = FALSE,
                show_rownames  = FALSE,
                annotation_col = annotation_col,
                silent         = TRUE)$gtable


# ── Plot 3D: Pathway Activity Heatmap ────────────────────────────────────────
# Dataset: Sheet "d_1" — pathway scores across cell types
# Goal: Compare pathway activity levels across conditions
# No clustering applied to preserve biological ordering

cell_data     <- read_excel(dest, sheet = "d_1")

# Convert to matrix and assign pathway names as row names
cell_data_mat <- as.matrix(cell_data[, 2:8])
rownames(cell_data_mat) <- cell_data$pathway
class(cell_data_mat) <- "numeric"

p3D <- pheatmap(cell_data_mat,
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                color        = colorRampPalette(c("firebrick", "white", "blue"))(100),
                silent       = TRUE)$gtable


# ── Plot 3E: Bubble Plot ──────────────────────────────────────────────────────
# Dataset: Sheet "e" — half life, alpha, count, and stage per gene
# Goal: Show relationship between half_life and alpha sized by count
# Point size scaled by count; color distinguishes 6h vs 72h stage

bubble_data <- read_excel(dest, sheet = "e")
bub_cols    <- c("6h" = "blue", "72h" = "green")

p3E <- capture_base_plot(quote({
  plot(bubble_data$half_life, bubble_data$alpha,
       col      = bub_cols[bubble_data$stage],
       pch      = 19,
       cex      = bubble_data$count * 0.1,   # scale point size by count
       xlab     = "Half Life",
       ylab     = "Alpha Life",
       cex.axis = 1.5,
       cex.lab  = 1.5)
  legend("topright", pch = 19,
         legend = c("72h", "6h"),
         col    = c("green", "blue"),
         cex    = 1.3, bty = "n")
  legend("bottomright", pch = 19,
         legend  = c(10, 20, 30),
         pt.cex  = c(1, 2, 3),
         title   = "Count",
         cex     = 1.3, bty = "n")
}))


# ── Plot 3F: Cell Type Proportion Stacked Barplot ────────────────────────────
# Dataset: Sheet "f" — cell type proportions at s00h and s72h stages
# Goal: Compare cell type composition between two time points
# legend = TRUE inside barplot() causes overlap, so we draw it manually

cell_type_stage_data <- read_excel(dest, sheet = "f")

# Keep only the two time points of interest
cell_type_stage_data <- cell_type_stage_data[
  cell_type_stage_data$stage %in% c("s00h", "s72h"), ]

# Create a cross-tabulation matrix: rows = cell types, cols = stages
freq <- tapply(cell_type_stage_data$proportion,
               list(cell_type_stage_data$cell_type,
                    cell_type_stage_data$stage), sum)

p3F <- capture_base_plot(quote({
  bp <- barplot(freq,
                col       = c("pink", "darkblue"),
                ylim      = c(0, 0.35),
                cex.axis  = 1.5,
                cex.lab   = 1.5,
                cex.names = 1.3)
  # Manual legend avoids overlap with bars (bty="n" removes legend border box)
  legend("topright",
         legend = rownames(freq),
         fill   = c("pink", "darkblue"),
         cex    = 1.3,
         bty    = "n")
}))


# ── Plot 3G: Gene Regulatory Network ─────────────────────────────────────────
# Dataset: Sheet "g" — adjacency matrix of gene interactions
# Goal: Visualize directed gene regulatory relationships
# Edge weight drives arrow size; zero-weight edges are removed

data_3G     <- read_excel(dest, sheet = "g")

# Extract row labels from first column and build adjacency matrix
rownames_3G <- data_3G[[1]]
data_3G.mat <- as.matrix(data_3G[, -1])
rownames(data_3G.mat) <- rownames_3G

# Build directed weighted graph from adjacency matrix
g <- graph_from_adjacency_matrix(data_3G.mat,
                                 mode     = "directed",
                                 weighted = TRUE,
                                 diag     = FALSE)

# Remove self-loops and unconnected edges (weight == 0)
g <- delete_edges(g, E(g)[weight == 0])

p3G <- capture_base_plot(quote({
  plot(g,
       edge.arrow.size  = E(g)$weight * 2,  # arrow size proportional to weight
       vertex.size      = 15,
       vertex.color     = "lightpink",
       vertex.label.cex = 1.2,
       layout           = layout_with_fr(g)) # Fruchterman-Reingold layout
}))


# ── Label Figure 3 Plots ──────────────────────────────────────────────────────
# Titles are baked into each grob so they appear in both individual and panel saves

p3A_l <- label_plot(p3A, "a) Cell Type Ratio")
p3B_l <- label_plot(p3B, "b) Half Life vs Alpha")
p3C_l <- label_plot(p3C, "c) Time-series Expression Heatmap")
p3D_l <- label_plot(p3D, "d) Pathway Activity Heatmap")
p3E_l <- label_plot(p3E, "e) Bubble Plot")
p3F_l <- label_plot(p3F, "f) Cell Type Proportion")
p3G_l <- label_plot(p3G, "g) Gene Regulatory Network")


# ── Save Individual Figure 3 Plots ───────────────────────────────────────────
# Each labelled plot saved at 300 dpi for publication quality

png("figure3_A_celltype_ratio.png", width = 2400, height = 1800, res = 300)
grid.newpage(); grid.draw(p3A_l); dev.off()

png("figure3_B_halflife_alpha.png", width = 2400, height = 1800, res = 300)
grid.newpage(); grid.draw(p3B_l); dev.off()

png("figure3_C_timeseries_heatmap.png", width = 2400, height = 1800, res = 300)
grid.newpage(); grid.draw(p3C_l); dev.off()

png("figure3_D_pathway_heatmap.png", width = 2400, height = 1800, res = 300)
grid.newpage(); grid.draw(p3D_l); dev.off()

png("figure3_E_bubble_plot.png", width = 2400, height = 1800, res = 300)
grid.newpage(); grid.draw(p3E_l); dev.off()

png("figure3_F_celltype_proportion.png", width = 2400, height = 1800, res = 300)
grid.newpage(); grid.draw(p3F_l); dev.off()

png("figure3_G_network.png", width = 2400, height = 1800, res = 300)
grid.newpage(); grid.draw(p3G_l); dev.off()


# ── Save Figure 3 Panel ───────────────────────────────────────────────────────
# layout_matrix uses a 6-column grid so each plot gets equal width:
#   Row 1: 3A spans cols 1-3, 3B spans cols 4-6
#   Row 2: 3C spans cols 1-3, 3D spans cols 4-6
#   Row 3: 3E spans cols 1-2, 3F spans cols 3-4, 3G spans cols 5-6

png("figure3_panel.png", width = 4800, height = 3600, res = 300)
grid.arrange(
  p3A_l, p3B_l,
  p3C_l, p3D_l,
  p3E_l, p3F_l, p3G_l,
  layout_matrix = rbind(c(1, 1, 1, 2, 2, 2),
                        c(3, 3, 3, 4, 4, 4),
                        c(5, 5, 6, 6, 7, 7))
)
dev.off()

# ── End of Script ─────────────────────────────────────────────────────────────
# Output files saved to working directory. Run getwd() to confirm location.
# Figure 1: 6 individual plots + 1 panel = 7 files
# Figure 3: 7 individual plots + 1 panel = 8 files
# Total: 15 PNG files# ══════════════════════════════════════════════════════════════════════════════
# MULTI-PANEL FIGURE GENERATION SCRIPT
# Description: This script generates two multi-panel figures (Figure 1 and
#              Figure 3) from biological datasets. It combines pheatmap and
#              base R plots into publication-ready panels using grid graphics.
# Author:      [Your Name]
# Date:        [Date]
# ══════════════════════════════════════════════════════════════════════════════


# ── SECTION 1: Install and Load Required Libraries ────────────────────────────
# These packages are required for plotting, layout, and data import.
# Run install.packages() only once; comment it out after first use.

install.packages(c("gridExtra",    # arranging multiple plots in a grid
                   "gridGraphics", # converting base R plots to grid objects
                   "pheatmap",     # heatmap generation
                   "png",          # reading PNG files into R
                   "readxl",       # reading Excel files
                   "igraph",       # network/graph plotting
                   "RColorBrewer"  # color palettes for heatmaps
))

library(grid)         # core grid graphics system
library(gridExtra)    # arrangeGrob and grid.arrange for panel layout
library(gridGraphics) # allows base R plots to be used in grid layout
library(pheatmap)     # heatmap with clustering
library(png)          # read PNG images
library(readxl)       # read Excel sheets
library(RColorBrewer) # color palettes
library(igraph)       # network graphs


# ── SECTION 2: Helper Functions ───────────────────────────────────────────────

# capture_base_plot():
# Base R plots cannot be directly used in grid.arrange() because they use
# a different graphics system. This function:
#   1. Renders the base R plot to a temporary PNG file at high resolution
#   2. Reads the PNG back into R as a raster image
#   3. Converts it to a grid-compatible rasterGrob object
# This makes base R plots behave like grid objects for panel assembly.

capture_base_plot <- function(expr) {
  tf <- tempfile(fileext = ".png")          # create a temporary file path
  png(tf, width = 1600, height = 1200, res = 150) # open PNG device at high res
  eval(expr)                                # evaluate and draw the plot
  dev.off()                                 # close the PNG device
  img <- png::readPNG(tf)                   # read the saved PNG back into R
  grid::rasterGrob(img,                     # convert to grid raster object
                   width  = unit(1, "npc"),
                   height = unit(1, "npc"))
}

# label_plot():
# Adds a bold label/title to the top-left of any grob (grid object).
# Works for both pheatmap gtables and rasterGrobs from capture_base_plot().
# Uses arrangeGrob() which is more compatible than grobTree() for gtables.

label_plot <- function(grob, label) {
  arrangeGrob(
    grob,
    top = textGrob(label,
                   x    = 0.02, just = "left",   # left-aligned position
                   gp   = gpar(fontsize = 14, fontface = "bold"))
  )
}


# ══════════════════════════════════════════════════════════════════════════════
#  FIGURE 1: Cancer and Gene Expression Analysis
#  Panels: Expression Heatmap | Volcano Plot | Scatter Plots | Density Plot
# ══════════════════════════════════════════════════════════════════════════════

# ── Plot 1A: Gene Expression Heatmap ─────────────────────────────────────────
# Dataset: Normalized counts from HBR/UHR RNA-seq experiment
# Goal: Visualize expression patterns of top differentially expressed genes
# Method: Hierarchical clustering on both rows (genes) and columns (samples)

link <- "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_top_deg_normalized_counts.csv"

expr_data <- read.delim(link, sep = ",", header = TRUE)

# Move gene names from first column into row names (required for matrix)
rownames(expr_data) <- expr_data[, 1]
expr_data <- expr_data[, -1]

# Convert data frame to numeric matrix (pheatmap requires a matrix)
mat <- as.matrix(expr_data)

# Generate heatmap; silent = TRUE suppresses auto-display so we can store it
# $gtable extracts the grid table object for use in panel assembly
pA <- pheatmap(mat,
               cluster_rows  = TRUE,
               cluster_cols  = TRUE,
               border_color  = "black",
               color         = colorRampPalette(c("white", "blue"))(100),
               fontsize_row  = 8,
               silent        = TRUE)$gtable


# ── Plot 1B: Volcano Plot ─────────────────────────────────────────────────────
# Dataset: DEG results with log2FoldChange and adjusted p-values
# Goal: Visualize significantly up/down regulated genes
# Dashed lines mark fold-change threshold (±1) and significance (p = 0.05)

expr_data2 <- read.delim(
  "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_deg_chr22_with_significance.csv",
  sep = ",", header = TRUE)

# Color map: green = upregulated, orange = downregulated, grey = not significant
vol_colors <- c("up" = "green", "down" = "orange", "ns" = "grey")

pB <- capture_base_plot(quote({
  plot(expr_data2$log2FoldChange, expr_data2$X.log10PAdj,
       col      = vol_colors[expr_data2$significance],
       pch      = 19,
       ylab     = "-log10PAdj",
       xlab     = "log2FoldChange",
       cex      = 1.2,
       cex.axis = 1.5,
       cex.lab  = 1.5)
  legend("topright", pch = 19,
         legend = c("down", "ns", "up"),
         col    = vol_colors[c("down", "ns", "up")],
         title  = "Significance", cex = 1.3)
  abline(v = c(-1, 1), lty = 2)           # fold change cutoff lines
  abline(h = -log10(0.05), lty = 2)       # significance threshold line
}))


# ── Plot 1C: Radius Mean vs Texture Mean Scatter ──────────────────────────────
# Dataset: Breast cancer diagnostic measurements
# Goal: Explore relationship between tumor radius and texture by diagnosis
# Color: Blue = Malignant (M), Orange = Benign (B)

breast_cancer_data <- read.delim(
  "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv",
  sep = ",", header = TRUE)

bc_cols <- c("B" = "orange", "M" = "blue")

pC <- capture_base_plot(quote({
  plot(breast_cancer_data$radius_mean, breast_cancer_data$texture_mean,
       col      = bc_cols[breast_cancer_data$diagnosis],
       pch      = 19,
       ylab     = "texture_mean",
       xlab     = "radius_mean",
       cex      = 1.2,
       cex.axis = 1.5,
       cex.lab  = 1.5)
  legend("topright", pch = 19,
         legend = c("M", "B"),
         col    = c("blue", "orange"),
         title  = "Diagnosis", cex = 1.3)
}))


# ── Plot 1D: Correlation Heatmap ──────────────────────────────────────────────
# Goal: Show pairwise Pearson correlations between 6 breast cancer features
# display_numbers = TRUE overlays the correlation values on each cell

features <- breast_cancer_data[, c("radius_mean", "texture_mean",
                                   "perimeter_mean", "area_mean",
                                   "smoothness_mean", "compactness_mean")]

# Compute pairwise correlation matrix
cor_mat <- cor(features, use = "complete.obs")

pD <- pheatmap(cor_mat,
               display_numbers = TRUE,
               number_format   = "%.1f",
               color           = colorRampPalette(c("white", "lightblue", "blue"))(100),
               border_color    = "black",
               main            = "",
               cluster_rows    = FALSE,
               cluster_cols    = FALSE,
               legend_breaks   = seq(-1, 1, by = 0.2),
               legend_labels   = seq(-1, 1, by = 0.2),
               silent          = TRUE)$gtable


# ── Plot 1E: Smoothness vs Compactness Scatter ────────────────────────────────
# Goal: Compare smoothness and compactness between malignant and benign tumors
# xaxt = "n" suppresses default x-axis so we can draw a custom one with axis()

pE <- capture_base_plot(quote({
  plot(breast_cancer_data$smoothness_mean, breast_cancer_data$compactness_mean,
       xlim     = c(0.05, 0.17),
       ylim     = c(0.01, 0.35),
       col      = bc_cols[breast_cancer_data$diagnosis],
       pch      = 19,
       ylab     = "compactness_mean",
       xlab     = "smoothness_mean",
       xaxt     = "n",              # suppress default x-axis
       cex      = 1.2,
       cex.axis = 1.5,
       cex.lab  = 1.5)
  axis(side = 1, seq(0.05, 0.15, 0.025), cex.axis = 1.5)  # custom x-axis
  legend("topleft", pch = 19,
         legend = c("M", "B"),
         col    = c("blue", "orange"),
         title  = "Diagnosis", cex = 1.3)
}))


# ── Plot 1F: Density Plot of Area Mean ────────────────────────────────────────
# Goal: Compare distribution of area_mean between malignant and benign tumors
# polygon() fills the area under the density curve with transparent color

# Compute density estimates for each diagnosis group
dens_M <- density(breast_cancer_data$area_mean[breast_cancer_data$diagnosis == "M"])
dens_B <- density(breast_cancer_data$area_mean[breast_cancer_data$diagnosis == "B"])

pF <- capture_base_plot(quote({
  plot(dens_M,
       col      = "blue",
       lwd      = 2,
       xlab     = "Area Mean",
       ylab     = "Density",
       ylim     = c(0, 0.0035),
       main     = "",
       cex.axis = 1.5,
       cex.lab  = 1.5)
  polygon(dens_M, col = rgb(0, 0, 1,   0.3), border = NA)  # blue fill, M
  polygon(dens_B, col = rgb(1, 0.5, 0, 0.3), border = NA)  # orange fill, B
  lines(dens_M, col = "blue",   lwd = 2)   # redraw lines on top of polygon
  lines(dens_B, col = "orange", lwd = 2)
  legend("topright", title = "Diagnosis",
         legend = c("M", "B"),
         fill   = c(rgb(0, 0, 1, 0.3), rgb(1, 0.5, 0, 0.3)),
         cex    = 1.3)
}))


# ── Label Figure 1 Plots ──────────────────────────────────────────────────────
# label_plot() bakes the title into each grob so it appears in both
# individual saves and the combined panel without extra steps

pA_l <- label_plot(pA, "a) Expression Heatmap")
pB_l <- label_plot(pB, "b) Expression Volcano Plot")
pC_l <- label_plot(pC, "c) Radius vs Texture")
pD_l <- label_plot(pD, "d) Correlation Heatmap")
pE_l <- label_plot(pE, "e) Smoothness vs Compactness")
pF_l <- label_plot(pF, "f) Area Distribution")


# ── Save Individual Figure 1 Plots ───────────────────────────────────────────
# Each plot is saved separately at 300 dpi for publication quality.
# grid.newpage() clears the canvas before drawing to avoid overlap.

png("figure1_A_expression_heatmap.png", width = 2400, height = 1800, res = 300)
grid.newpage(); grid.draw(pA_l); dev.off()

png("figure1_B_volcano.png", width = 2400, height = 1800, res = 300)
grid.newpage(); grid.draw(pB_l); dev.off()

png("figure1_C_radius_texture.png", width = 2400, height = 1800, res = 300)
grid.newpage(); grid.draw(pC_l); dev.off()

png("figure1_D_correlation_heatmap.png", width = 2400, height = 1800, res = 300)
grid.newpage(); grid.draw(pD_l); dev.off()

png("figure1_E_smoothness_compactness.png", width = 2400, height = 1800, res = 300)
grid.newpage(); grid.draw(pE_l); dev.off()

png("figure1_F_density.png", width = 2400, height = 1800, res = 300)
grid.newpage(); grid.draw(pF_l); dev.off()


# ── Save Figure 1 Panel ───────────────────────────────────────────────────────
# Arrange all 6 plots in a 3-column grid panel and save as one PNG

png("figure1_panel.png", width = 3600, height = 2400, res = 300)
grid.arrange(pA_l, pB_l, pC_l,
             pD_l, pE_l, pF_l,
             ncol = 3)
dev.off()


# ══════════════════════════════════════════════════════════════════════════════
#  FIGURE 3: Immune Cell and Gene Regulation Analysis
#  Panels: Boxplot | Quadrant Scatter | Heatmaps | Bubble | Barplot | Network
# ══════════════════════════════════════════════════════════════════════════════

# Download the Excel workbook containing all Figure 3 datasets
# Each sheet corresponds to a different sub-analysis

url  <- "https://github.com/HackBio-Internship/2025_project_collection/raw/refs/heads/main/hb_stage_2.xlsx"
dest <- tempfile(fileext = ".xlsx")    # save to a temporary file
download.file(url, dest, mode = "wb") # mode "wb" = write binary (needed for Excel)


# ── Plot 3A: Cell Type Ratio Boxplot ─────────────────────────────────────────
# Dataset: Sheet "a" — cell type ratios across samples
# Goal: Compare ratio distributions across immune cell types
# las = 2 rotates x-axis labels vertically for readability

cell_type_ratio_data <- read_excel(dest, sheet = "a")
cell_type_ratio_data$cell_type <- as.factor(cell_type_ratio_data$cell_type)

# Named color vector for each cell type
box_cols <- c(B          = "darkblue",
              Basophil   = "lightgreen",
              DC         = "brown",
              HSPC       = "pink",
              Macrophage = "lightblue",
              Monocyte   = "darkgreen",
              Neutrophil = "purple",
              NK         = "orange",
              "Pre B"    = "darkorange",
              T          = "yellow")

p3A <- capture_base_plot(quote({
  boxplot(new_ratio ~ cell_type,
          data     = cell_type_ratio_data,
          ylab     = "Ratio",
          xlab     = "",
          las      = 2,
          col      = box_cols[levels(cell_type_ratio_data$cell_type)],
          cex.axis = 1.3,
          cex.lab  = 1.5)
}))


# ── Plot 3B: Half Life vs Alpha Quadrant Scatter ──────────────────────────────
# Dataset: Sheet "b" — gene half-life and alpha values
# Goal: Segment genes into 4 quadrants based on log2 thresholds
# Ccr2 labels Q1 (top-left), Camp labels Q4 (bottom-right)

alpha_half_life_data <- read_excel(dest, sheet = "b")

# Log2 transform both axes
x <- log2(alpha_half_life_data$half_life)
y <- log2(alpha_half_life_data$alpha)

# Assign each point to a quadrant based on threshold values
quad <- ifelse(x <= 2.5 & y >= -3.6, "Q1",
               ifelse(x >  2.5 & y >= -3.6, "Q2",
                      ifelse(x <= 2.5 & y <  -3.6, "Q3", "Q4")))

quad_cols <- c(Q1 = "green", Q2 = "red", Q3 = "darkgrey", Q4 = "blue")

p3B <- capture_base_plot(quote({
  plot(x, y,
       ylim     = c(-10, 0),
       col      = quad_cols[quad],
       pch      = 19,
       cex      = 0.8,
       xlab     = "log2(Half Life)",
       ylab     = "log2(Alpha)",
       cex.axis = 1.5,
       cex.lab  = 1.5)
  abline(v = 2.5,  lty = 2)   # vertical threshold
  abline(h = -3.6, lty = 2)   # horizontal threshold
  # Label Q1 (top-left quadrant) with gene name Ccr2
  text((min(x) + 2.5) / 2, (-3.6 + 0) / 2,
       "Ccr2", cex = 1.8, font = 2, col = "black")
  # Label Q4 (bottom-right quadrant) with gene name Camp
  text((2.5 + max(x)) / 2, (-10 + -3.6) / 2,
       "Camp", cex = 1.8, font = 2, col = "black")
}))


# ── Plot 3C: Time-series Expression Heatmap ───────────────────────────────────
# Dataset: Sheet "c" — gene expression across cell types and time points
# Goal: Visualize expression dynamics over time without column clustering
# Annotation bar added at top to show Time and CellType groupings

data_3C     <- read_excel(dest, sheet = "c")

# Convert tibble slice to numeric matrix (pheatmap requires numeric matrix)
data_3C_mat <- as.matrix(data_3C[, 2:50])
class(data_3C_mat) <- "numeric"

# Build annotation data frame from column names
# Assumes last 3 characters = time point, rest = cell type
annotation_col <- data.frame(
  Time     = substr(colnames(data_3C_mat),
                    nchar(colnames(data_3C_mat)) - 2,
                    nchar(colnames(data_3C_mat))),
  CellType = substr(colnames(data_3C_mat),
                    1,
                    nchar(colnames(data_3C_mat)) - 4),
  row.names = colnames(data_3C_mat)
)

p3C <- pheatmap(data_3C_mat,
                cluster_cols   = FALSE,         # preserve time order
                color          = colorRampPalette(brewer.pal(9, "Blues"))(200),
                show_colnames  = FALSE,
                show_rownames  = FALSE,
                annotation_col = annotation_col,
                silent         = TRUE)$gtable


# ── Plot 3D: Pathway Activity Heatmap ────────────────────────────────────────
# Dataset: Sheet "d_1" — pathway scores across cell types
# Goal: Compare pathway activity levels across conditions
# No clustering applied to preserve biological ordering

cell_data     <- read_excel(dest, sheet = "d_1")

# Convert to matrix and assign pathway names as row names
cell_data_mat <- as.matrix(cell_data[, 2:8])
rownames(cell_data_mat) <- cell_data$pathway
class(cell_data_mat) <- "numeric"

p3D <- pheatmap(cell_data_mat,
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                color        = colorRampPalette(c("firebrick", "white", "blue"))(100),
                silent       = TRUE)$gtable


# ── Plot 3E: Bubble Plot ──────────────────────────────────────────────────────
# Dataset: Sheet "e" — half life, alpha, count, and stage per gene
# Goal: Show relationship between half_life and alpha sized by count
# Point size scaled by count; color distinguishes 6h vs 72h stage

bubble_data <- read_excel(dest, sheet = "e")
bub_cols    <- c("6h" = "blue", "72h" = "green")

p3E <- capture_base_plot(quote({
  plot(bubble_data$half_life, bubble_data$alpha,
       col      = bub_cols[bubble_data$stage],
       pch      = 19,
       cex      = bubble_data$count * 0.1,   # scale point size by count
       xlab     = "Half Life",
       ylab     = "Alpha Life",
       cex.axis = 1.5,
       cex.lab  = 1.5)
  legend("topright", pch = 19,
         legend = c("72h", "6h"),
         col    = c("green", "blue"),
         cex    = 1.3, bty = "n")
  legend("bottomright", pch = 19,
         legend  = c(10, 20, 30),
         pt.cex  = c(1, 2, 3),
         title   = "Count",
         cex     = 1.3, bty = "n")
}))


# ── Plot 3F: Cell Type Proportion Stacked Barplot ────────────────────────────
# Dataset: Sheet "f" — cell type proportions at s00h and s72h stages
# Goal: Compare cell type composition between two time points
# legend = TRUE inside barplot() causes overlap, so we draw it manually

cell_type_stage_data <- read_excel(dest, sheet = "f")

# Keep only the two time points of interest
cell_type_stage_data <- cell_type_stage_data[
  cell_type_stage_data$stage %in% c("s00h", "s72h"), ]

# Create a cross-tabulation matrix: rows = cell types, cols = stages
freq <- tapply(cell_type_stage_data$proportion,
               list(cell_type_stage_data$cell_type,
                    cell_type_stage_data$stage), sum)

p3F <- capture_base_plot(quote({
  bp <- barplot(freq,
                col       = c("pink", "darkblue"),
                ylim      = c(0, 0.35),
                cex.axis  = 1.5,
                cex.lab   = 1.5,
                cex.names = 1.3)
  # Manual legend avoids overlap with bars (bty="n" removes legend border box)
  legend("topright",
         legend = rownames(freq),
         fill   = c("pink", "darkblue"),
         cex    = 1.3,
         bty    = "n")
}))


# ── Plot 3G: Gene Regulatory Network ─────────────────────────────────────────
# Dataset: Sheet "g" — adjacency matrix of gene interactions
# Goal: Visualize directed gene regulatory relationships
# Edge weight drives arrow size; zero-weight edges are removed

data_3G     <- read_excel(dest, sheet = "g")

# Extract row labels from first column and build adjacency matrix
rownames_3G <- data_3G[[1]]
data_3G.mat <- as.matrix(data_3G[, -1])
rownames(data_3G.mat) <- rownames_3G

# Build directed weighted graph from adjacency matrix
g <- graph_from_adjacency_matrix(data_3G.mat,
                                 mode     = "directed",
                                 weighted = TRUE,
                                 diag     = FALSE)

# Remove self-loops and unconnected edges (weight == 0)
g <- delete_edges(g, E(g)[weight == 0])

p3G <- capture_base_plot(quote({
  plot(g,
       edge.arrow.size  = E(g)$weight * 2,  # arrow size proportional to weight
       vertex.size      = 15,
       vertex.color     = "lightpink",
       vertex.label.cex = 1.2,
       layout           = layout_with_fr(g)) # Fruchterman-Reingold layout
}))


# ── Label Figure 3 Plots ──────────────────────────────────────────────────────
# Titles are baked into each grob so they appear in both individual and panel saves

p3A_l <- label_plot(p3A, "a) Fig 2a: Cell Type Ratio Distribution")
p3B_l <- label_plot(p3B, "b) Fig 2b: Half Life vs Alpha")
p3C_l <- label_plot(p3C, "c) Fig 2c: Time-series Expression Heatmap")
p3D_l <- label_plot(p3D, "d) Fig 2d: Pathway Enrichment Heatmap")
p3E_l <- label_plot(p3E, "e) Fig 2e: Kinetic Regimes")
p3F_l <- label_plot(p3F, "f) Fig 2f: Cell Type Proportion")
p3G_l <- label_plot(p3G, "g) Fig 2g: Directed cell–cell interaction network")


# ── Save Individual Figure 3 Plots ───────────────────────────────────────────
# Each labelled plot saved at 300 dpi for publication quality

png("figure3_A_celltype_ratio.png", width = 2400, height = 1800, res = 300)
grid.newpage(); grid.draw(p3A_l); dev.off()

png("figure3_B_halflife_alpha.png", width = 2400, height = 1800, res = 300)
grid.newpage(); grid.draw(p3B_l); dev.off()

png("figure3_C_timeseries_heatmap.png", width = 2400, height = 1800, res = 300)
grid.newpage(); grid.draw(p3C_l); dev.off()

png("figure3_D_pathway_heatmap.png", width = 2400, height = 1800, res = 300)
grid.newpage(); grid.draw(p3D_l); dev.off()

png("figure3_E_bubble_plot.png", width = 2400, height = 1800, res = 300)
grid.newpage(); grid.draw(p3E_l); dev.off()

png("figure3_F_celltype_proportion.png", width = 2400, height = 1800, res = 300)
grid.newpage(); grid.draw(p3F_l); dev.off()

png("figure3_G_network.png", width = 2400, height = 1800, res = 300)
grid.newpage(); grid.draw(p3G_l); dev.off()


# ── Save Figure 3 Panel ───────────────────────────────────────────────────────
# layout_matrix uses a 6-column grid so each plot gets equal width:
#   Row 1: 3A spans cols 1-3, 3B spans cols 4-6
#   Row 2: 3C spans cols 1-3, 3D spans cols 4-6
#   Row 3: 3E spans cols 1-2, 3F spans cols 3-4, 3G spans cols 5-6

png("figure3_panel.png", width = 4800, height = 3600, res = 300)
grid.arrange(
  p3A_l, p3B_l,
  p3C_l, p3D_l,
  p3E_l, p3F_l, p3G_l,
  layout_matrix = rbind(c(1, 1, 1, 2, 2, 2),
                        c(3, 3, 3, 4, 4, 4),
                        c(5, 5, 6, 6, 7, 7))
)
dev.off()

# ── End of Script ─────────────────────────────────────────────────────────────
# Output files saved to working directory. Run getwd() to confirm location.
# Figure 1: 6 individual plots + 1 panel = 7 files
# Figure 3: 7 individual plots + 1 panel = 8 files
# Total: 15 PNG files