# --- 0. SETUP OUTPUT DIRECTORY ---
# Create a directory to store the output tables if it doesn't exist
output_dir <- "output_tables"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  message(paste("Created directory:", output_dir))
} else {
  message(paste("Output directory", output_dir, "already exists."))
}

# Create a directory to store the public output figures and tables if it doesn't exist
public_output_dir <- "../public" # Define the public output directory path
if (!dir.exists(public_output_dir)) {
  dir.create(public_output_dir, recursive = TRUE) # Use recursive = TRUE if parent dirs might not exist
  message(paste("Created directory:", public_output_dir))
} else {
  message(paste("Public output directory", public_output_dir, "already exists."))
}


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
library(tidyverse)# For data manipulation (dplyr, tidyr) and ggplot2
library(pheatmap) # For pretty heatmaps
library(RColorBrewer) # For color palettes
library(DESeq2) # For RNA-Seq data analysis, including normalization and transformations
library(ggplot2) # For plotting
library(ggrepel) # For non-overlapping text labels in ggplot

# --- 2. IMPORT DATA ---
# Define the file path
file_path <- "../public/GeneCountMatrix_Original.csv" # Adjust if your file is elsewhere

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
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
message(paste("Number of genes after filtering low counts:", nrow(dds)))


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

# --- 5.1 EXTRACT: PCA DATA ---
cat("\n--- PCA Data Table ---\n")
print(pca_data)
pca_data_filename <- file.path(output_dir, "pca_coordinates.csv")
write.csv(pca_data, pca_data_filename, row.names = TRUE) # Save full table
cat(paste("PCA coordinates saved to:", pca_data_filename, "\n"))
# --- NEW ---
# Export PCA table to public directory
pca_data_public_filename <- file.path(public_output_dir, "Outlier_PCA.csv")
write.csv(pca_data, pca_data_public_filename, row.names = TRUE)
cat(paste("PCA coordinates also saved to public directory:", pca_data_public_filename, "\n"))
# --- END NEW ---

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  geom_text_repel(size = 3, max.overlaps = Inf) + # Add sample labels without overlap
  ggtitle("PCA Plot of Samples (VST Data)") +
  theme_bw() +
  scale_color_brewer(palette = "Set1") # Use a distinct color palette

print(pca_plot)
# --- NEW ---
# Export PCA figure to public directory
pca_figure_public_filename <- file.path(public_output_dir, "Outlier_PCA.png")
ggsave(pca_figure_public_filename, plot = pca_plot, width = 8, height = 6)
cat(paste("PCA plot saved to public directory:", pca_figure_public_filename, "\n"))
# --- END NEW ---

cat("\n--- PCA Interpretation ---\n")
cat("Look for samples that cluster far away from other replicates of the same condition.\n")
cat("Samples lying distant from the main cluster(s) might be outliers.\n\n")


## Method 5.2: Hierarchical Clustering based on Sample Distances ##
# Calculates distances between samples and clusters them. Outliers may not cluster tightly with replicates.

# Calculate pairwise sample distances using Euclidean distance on VST data
sample_dists <- dist(t(assay(vsd))) # Transpose needed: rows=samples, cols=genes

# Create distance matrix object
sample_dist_matrix <- as.matrix(sample_dists)

# --- 5.2 EXTRACT: SAMPLE DISTANCE MATRIX ---
cat("\n--- Sample Distance Matrix ---\n")
print(round(sample_dist_matrix, 2)) # Print rounded matrix for preview
dist_matrix_filename <- file.path(output_dir, "sample_distance_matrix.csv")
write.csv(sample_dist_matrix, dist_matrix_filename, row.names = TRUE)
cat(paste("Sample distance matrix saved to:", dist_matrix_filename, "\n"))
# --- NEW ---
# Export Euclidean distance table to public directory
dist_matrix_public_filename <- file.path(public_output_dir, "Outlier_Heatmap_Euclidean_Distance.csv")
write.csv(sample_dist_matrix, dist_matrix_public_filename, row.names = TRUE)
cat(paste("Sample distance matrix also saved to public directory:", dist_matrix_public_filename, "\n"))
# --- END NEW ---

# Define colors for annotation
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# --- NEW ---
# Define filename for the heatmap plot
dist_heatmap_public_filename <- file.path(public_output_dir, "Outlier_Heatmap_Euclidean_Distance.png")
# --- END NEW ---

# Plot heatmap of sample distances
dist_heatmap <- pheatmap(sample_dist_matrix,
                         clustering_distance_rows = sample_dists,
                         clustering_distance_cols = sample_dists,
                         col = colors,
                         annotation_col = sample_info,
                         annotation_row = sample_info,
                         main = "Heatmap of Sample-to-Sample Euclidean Distances (VST Data)",
                         # --- NEW ---
                         filename = dist_heatmap_public_filename, # Save the plot directly
                         width = 8, height = 7, # Adjust dimensions as needed
                         silent = TRUE # Prevent pheatmap from printing the plot object details
                         # --- END NEW ---
)

# Note: pheatmap returns the plot object invisibly when filename is not used.
# print(dist_heatmap) # Uncomment if running interactively and need explicit print when not saving to file

# --- NEW ---
cat(paste("Euclidean distance heatmap saved to public directory:", dist_heatmap_public_filename, "\n"))
# --- END NEW ---

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

# --- 6.1 EXTRACT: SAMPLE CORRELATION MATRIX ---
cat("\n--- Sample Correlation Matrix ---\n")
print(round(sample_cor, 3)) # Print rounded matrix for preview
cor_matrix_filename <- file.path(output_dir, "sample_correlation_matrix.csv")
write.csv(sample_cor, cor_matrix_filename, row.names = TRUE)
cat(paste("Sample correlation matrix saved to:", cor_matrix_filename, "\n"))
# --- NEW ---
# Export Pearson correlation table to public directory
cor_matrix_public_filename <- file.path(public_output_dir, "Outlier_Heatmap_Pearson_Correlation.csv")
write.csv(sample_cor, cor_matrix_public_filename, row.names = TRUE)
cat(paste("Sample correlation matrix also saved to public directory:", cor_matrix_public_filename, "\n"))
# --- END NEW ---

# --- NEW ---
# Define filename for the correlation heatmap plot
cor_heatmap_public_filename <- file.path(public_output_dir, "Outlier_Heatmap_Pearson_Correlation.png")
# --- END NEW ---

# Plot heatmap of sample correlations
cor_heatmap <- pheatmap(sample_cor,
                        annotation_col = sample_info,
                        annotation_row = sample_info,
                        main = "Heatmap of Sample-to-Sample Pearson Correlation (VST Data)",
                        annotation_names_row = FALSE, # Avoid duplicating labels if annotation_row is used
                        annotation_names_col = FALSE,
                        # --- NEW ---
                        filename = cor_heatmap_public_filename, # Save the plot directly
                        width = 8, height = 7, # Adjust dimensions as needed
                        silent = TRUE # Prevent pheatmap from printing the plot object details
                        # --- END NEW ---
)

# print(cor_heatmap) # Uncomment if running interactively when not saving to file

# --- NEW ---
cat(paste("Pearson correlation heatmap saved to public directory:", cor_heatmap_public_filename, "\n"))
# --- END NEW ---

cat("\n--- Sample Correlation Heatmap Interpretation ---\n")
cat("This heatmap visualizes the Pearson correlation coefficients between samples.\n")
cat("Look for samples that exhibit lower correlation values (lighter colors or specified by legend) with their biological replicates compared to the correlations *among* other replicates in that group.\n")
cat("A block of high correlation (darker colors) is expected within each condition group. Samples outside this pattern might be outliers.\n\n")


## Method 6.2 (Optional but informative): Inter-sample distance distribution ##
# Calculate the mean distance of each sample to all *other* samples. Outliers might have a higher average distance.

# Use the sample_dist_matrix calculated earlier in section 5.2
avg_distances <- apply(sample_dist_matrix, 1, mean)
avg_dist_df <- data.frame(sample = names(avg_distances), avg_dist = avg_distances)
avg_dist_df <- merge(avg_dist_df, sample_info, by.x="sample", by.y="row.names") # Add condition info

# --- 6.2 EXTRACT: AVERAGE SAMPLE DISTANCE DATA ---
cat("\n--- Average Sample Distance Table ---\n")
print(avg_dist_df) # Print the data frame
avg_dist_filename <- file.path(output_dir, "average_sample_distances.csv")
# Save without R's default row numbers, as sample ID is already a column
write.csv(avg_dist_df, avg_dist_filename, row.names = FALSE)
cat(paste("Average sample distance data saved to:", avg_dist_filename, "\n"))
# --- NEW ---
# Export average distance table to public directory
avg_dist_public_filename <- file.path(public_output_dir, "Outlier_Bar_Euclidean_Distance.csv")
write.csv(avg_dist_df, avg_dist_public_filename, row.names = FALSE)
cat(paste("Average sample distance data also saved to public directory:", avg_dist_public_filename, "\n"))
# --- END NEW ---

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
# --- NEW ---
# Export Euclidean distance bar plot figure to public directory
dist_barplot_public_filename <- file.path(public_output_dir, "Outlier_Bar_Euclidean_Distance.png")
ggsave(dist_barplot_public_filename, plot = dist_barplot, width = 7, height = 6)
cat(paste("Average Euclidean distance bar plot saved to public directory:", dist_barplot_public_filename, "\n"))
# --- END NEW ---

cat("\n--- Average Sample Distance Interpretation ---\n")
cat("This plot shows the average Euclidean distance of each sample to all other samples.\n")
cat("While not definitive on its own, samples with notably higher average distances than others (especially within their condition group) warrant closer inspection as potential outliers.\n\n")


# --- 7. SUMMARY AND CONCLUSION ---

cat("\n--- Overall Outlier Assessment ---\n")
cat("Review the PCA plot, the sample distance heatmap, the correlation heatmap, and the average distance plot.\n")
cat("Identify samples that consistently appear separated or dissimilar from their replicates across multiple visualizations.\n")
cat("For example:\n")
cat("  - Is a sample far from its group on the PCA plot?\n")
cat("  - Does it show high distance / low correlation to its replicates in the heatmaps?\n")
cat("  - Does it have a notably high average distance to other samples?\n")
cat(paste("Numerical data underlying these plots has been saved to the '", output_dir, "' directory.\n", sep=""))
# --- NEW ---
cat(paste("Visualizations and corresponding data tables have also been saved to the '", public_output_dir, "' directory.\n", sep=""))
# --- END NEW ---
cat("Based on the combined evidence, decide which samples, if any, are strong candidates for being outliers.\n")
cat("The decision to remove an outlier often depends on the experimental context and the potential reason for the outlier behavior.\n")

# Note: This script identifies potential outliers based on sample similarity.
# Further investigation (e.g., checking library size, mapping rates, experimental notes)
# is often needed to confirm and understand why a sample might be an outlier.
# Packages like DESeq2 also have internal outlier detection (e.g., Cook's distance) during
# differential expression analysis, which provides another layer of checks.

# End of Script