# --- 1. INSTALL AND LOAD NECESSARY LIBRARIES ---

# Check if BiocManager is installed, if not, install it
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# List of required packages
required_packages <- c("tidyverse", "pheatmap", "RColorBrewer", "DESeq2", "ggplot2", "ggrepel")

# Check and install CRAN packages
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    # Only install DESeq2 via BiocManager, others via install.packages
    if (pkg == "DESeq2") {
      message("Installing DESeq2 using BiocManager...")
      BiocManager::install("DESeq2")
    } else if (!pkg %in% rownames(installed.packages())) {
      message(paste("Installing", pkg, "from CRAN..."))
      install.packages(pkg)
    }
  }
}

# Load libraries
library(tidyverse)    # For data manipulation (dplyr, tidyr) and ggplot2
library(pheatmap)     # For pretty heatmaps
library(RColorBrewer) # For color palettes
library(DESeq2)       # For RNA-Seq data analysis, including normalization and transformations
library(ggplot2)      # For plotting
library(ggrepel)      # For non-overlapping text labels in ggplot

# --- 2. IMPORT DATA ---
# Define the file path
file_path <- "../public/GeneCountMatrix_Original.csv"

# Read the CSV file
# Check the first few lines to confirm structure if needed: readLines(file_path, n=5)
count_data <- read.csv(file_path, row.names = "gene_id", stringsAsFactors = FALSE) # Use gene_id as row names

# Display the first few rows and columns to verify
print("Original Count Data Head:")
print(head(count_data))

# Separate gene symbols (if needed later, though not used for outlier detection itself)
gene_symbols <- count_data$gene_symbol
count_matrix <- count_data[, -1] # Remove the gene_symbol column to keep only counts

# Ensure count matrix contains only numeric integer counts (DESeq2 requires integers)
# The example shows fractional counts, which is unusual for raw counts (might be TPMs/FPKMs?).
# Assuming they *should* be counts, we'll round them for DESeq2 compatibility.
# If they are *not* raw counts, DESeq2's transformations might not be the most appropriate,
# but PCA/clustering on log-transformed normalized values (like TPM/FPKM) is still valid.
# Let's proceed assuming they represent counts and need rounding.
count_matrix <- round(count_matrix)

print("Count Matrix Head (numeric, rounded):")
print(head(count_matrix))
print(paste("Dimensions of count matrix:", dim(count_matrix)[1], "genes,", dim(count_matrix)[2], "samples"))


# --- 3. PREPARE METADATA (SAMPLE INFORMATION) ---

# Create a data frame describing the samples
# Infer conditions from column names
sample_names <- colnames(count_matrix)
conditions <- gsub("_R[0-9]+$", "", sample_names) # Remove replicate number suffix

sample_info <- data.frame(
  row.names = sample_names,
  condition = factor(conditions) # Factor is important for DESeq2 and plotting
)

print("Sample Metadata:")
print(sample_info)

# Check if column names in count_matrix match row names in sample_info
if (!all(rownames(sample_info) == colnames(count_matrix))) {
  stop("Mismatch between count matrix column names and sample info row names!")
}


# --- 4. PERFORM DATA TRANSFORMATION FOR VISUALIZATION ---
# It's crucial to transform count data to stabilize variance across the range of mean expression
# DESeq2's variance stabilizing transformation (VST) or regularized logarithm (rlog) are recommended.
# VST is generally faster for larger datasets.

# Create DESeqDataSet object (requires integer counts)
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_info,
                              design = ~ condition) # Basic design formula

# Optional: Pre-filter low-count genes (can sometimes improve visualization)
# Keep rows with at least (e.g.) 10 reads total
# keep <- rowSums(counts(dds)) >= 10
# dds <- dds[keep,]

# Apply Variance Stabilizing Transformation
# Use blind=TRUE for QC plots where you don't want prior condition info biasing the transformation
vsd <- vst(dds, blind = TRUE)
# Alternatively, use rlog: rld <- rlog(dds, blind = TRUE) # Can be slow for many samples

print("Transformed Data Head (using VST):")
print(head(assay(vsd)))


# --- 5. VISUALIZATION-BASED OUTLIER DETECTION ---

## Method 5.1: Principal Component Analysis (PCA) ##
# PCA reduces dimensionality and helps visualize sample clustering. Outliers often appear distant from their group.

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  geom_text_repel(size = 3, max.overlaps = Inf) + # Add sample labels without overlap
  ggtitle("PCA Plot of Samples (VST Data)") +
  theme_bw() +
  scale_color_brewer(palette = "Set1") # Use a distinct color palette

print(pca_plot)
cat("\n--- PCA Interpretation ---\n")
cat("Look for samples that cluster far away from other replicates of the same condition.\n")
cat("Samples lying distant from the main cluster(s) might be outliers.\n\n")


## Method 5.2: Hierarchical Clustering based on Sample Distances ##
# Calculates distances between samples and clusters them. Outliers may not cluster tightly with replicates.

# Calculate pairwise sample distances using Euclidean distance on VST data
sample_dists <- dist(t(assay(vsd))) # Transpose needed: rows=samples, cols=genes

# Create distance matrix object
sample_dist_matrix <- as.matrix(sample_dists)

# Define colors for annotation
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# Plot heatmap of sample distances
dist_heatmap <- pheatmap(sample_dist_matrix,
                         clustering_distance_rows = sample_dists,
                         clustering_distance_cols = sample_dists,
                         col = colors,
                         annotation_col = sample_info,
                         annotation_row = sample_info,
                         main = "Heatmap of Sample-to-Sample Euclidean Distances (VST Data)")

print(dist_heatmap)
cat("\n--- Sample Distance Heatmap Interpretation ---\n")
cat("This heatmap visualizes the distance matrix calculated above.\n")
cat("Look for samples (rows/columns) that show large distances (lighter colors) to their expected replicates.\n")
cat("An outlier might appear as a row/column that is distinctly different in color compared to its replicate group.\n")
cat("The dendrograms also show clustering: check if any sample branches off much earlier than its replicates.\n\n")


# --- 6. STATISTICAL ANALYSIS-BASED OUTLIER DETECTION ---

## Method 6.1: Sample Correlation Heatmap ##
# Calculates pairwise correlations between samples. Outliers often show lower correlation with replicates.

# Calculate pairwise Pearson correlations on VST data
sample_cor <- cor(assay(vsd), method = "pearson")

# Plot heatmap of sample correlations
cor_heatmap <- pheatmap(sample_cor,
                        annotation_col = sample_info,
                        annotation_row = sample_info,
                        main = "Heatmap of Sample-to-Sample Pearson Correlation (VST Data)",
                        annotation_names_row = FALSE, # Avoid duplicating labels if annotation_row is used
                        annotation_names_col = FALSE)

print(cor_heatmap)
cat("\n--- Sample Correlation Heatmap Interpretation ---\n")
cat("This heatmap visualizes the Pearson correlation coefficients between samples.\n")
cat("Look for samples that exhibit lower correlation values (lighter colors or specified by legend) with their biological replicates compared to the correlations *among* other replicates in that group.\n")
cat("A block of high correlation (darker colors) is expected within each condition group. Samples outside this pattern might be outliers.\n\n")


## Method 6.2 (Optional but informative): Inter-sample distance distribution ##
# Calculate the mean distance of each sample to all *other* samples. Outliers might have a higher average distance.

avg_distances <- apply(sample_dist_matrix, 1, mean)
avg_dist_df <- data.frame(sample = names(avg_distances), avg_dist = avg_distances)
avg_dist_df <- merge(avg_dist_df, sample_info, by.x="sample", by.y="row.names") # Add condition info

# Simple bar plot of average distances
dist_barplot <- ggplot(avg_dist_df, aes(x = reorder(sample, avg_dist), y = avg_dist, fill = condition)) +
  geom_bar(stat = "identity") +
  coord_flip() + # Flip for better readability of sample names
  labs(title = "Average Euclidean Distance to Other Samples (VST Data)",
       x = "Sample",
       y = "Average Distance") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1")

print(dist_barplot)
cat("\n--- Average Sample Distance Interpretation ---\n")
cat("This plot shows the average Euclidean distance of each sample to all other samples.\n")
cat("While not definitive on its own, samples with notably higher average distances than others (especially within their condition group) warrant closer inspection as potential outliers.\n\n")


# --- 7. SUMMARY AND CONCLUSION ---

cat("\n--- Overall Outlier Assessment ---\n")
cat("Review the PCA plot, the sample distance heatmap, the correlation heatmap, and the average distance plot.\n")
cat("Identify samples that consistently appear separated or dissimilar from their replicates across multiple visualizations.\n")
cat("For example:\n")
cat("  - Is a sample far from its group on the PCA plot?\n")
cat("  - Does it show high distance / low correlation to its replicates in the heatmaps?\n")
cat("  - Does it have a notably high average distance to other samples?\n")
cat("Based on the combined evidence, decide which samples, if any, are strong candidates for being outliers.\n")
cat("The decision to remove an outlier often depends on the experimental context and the potential reason for the outlier behavior.\n")

# Note: This script identifies potential outliers based on sample similarity.
# Further investigation (e.g., checking library size, mapping rates, experimental notes)
# is often needed to confirm and understand why a sample might be an outlier.
# Packages like DESeq2 also have internal outlier detection (e.g., Cook's distance) during
# differential expression analysis, which provides another layer of checks.

# End of Script