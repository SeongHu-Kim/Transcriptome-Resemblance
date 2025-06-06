# --- 1. Install and Load Required Libraries ---
# Ensure necessary packages are installed
# Run these lines in the R console if you haven't installed them before:
#
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#
# BiocManager::install("DESeq2")
# BiocManager::install("limma")
# install.packages("ggplot2")
# install.packages("pheatmap")
# install.packages("RColorBrewer")
# install.packages("tibble")

library(DESeq2)
library(limma)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(tibble)

# --- 2. Load Data ---
# Set the path to your combined count matrix file
# <<< ADJUST THIS PATH TO YOUR FILE >>>
count_file <- "../public/GeneCountMatrix_Combined.csv"

# Read the count data
# Assuming gene IDs are in the first column and should be row names
count_data <- read.csv(count_file, header = TRUE, row.names = 1, check.names = FALSE)

# Separate gene symbols (optional, but good practice)
gene_symbols <- count_data$gene_symbol
count_matrix <- count_data[, -1] # Remove the gene_symbol column

# Ensure counts are integers (DESeq2 prefers integers)
if(any(count_matrix != floor(count_matrix))){
  print("Warning: Non-integer values found in count matrix. Rounding for DESeq2.")
  count_matrix <- round(count_matrix)
}
count_matrix <- as.matrix(count_matrix)

# --- 3. Prepare Metadata (Sample Information) - REVISED AGAIN ---
# Create a data frame describing the samples using explicit logic based on prefixes

sample_names <- colnames(count_matrix)
conditions <- character(length(sample_names)) # Initialize empty vector

# Assign conditions based on known sample name prefixes
conditions[grepl("^Control_", sample_names)] <- "Control"
conditions[grepl("^Experimental_", sample_names)] <- "Experimental"
conditions[grepl("^SLE\\.active", sample_names)] <- "SLE.active"      # Escaped '.'
conditions[grepl("^SLE\\.less\\.active", sample_names)] <- "SLE.less.active" # Escaped '.'
conditions[grepl("^Ctrl_LHC", sample_names)] <- "Ctrl_LHC"

# Check if any sample didn't match a pattern
if(any(conditions == "")){
  unmatched_samples <- sample_names[conditions == ""]
  stop("ERROR: Some sample names did not match expected condition patterns: ",
       paste(unmatched_samples, collapse=", "))
}

# Define the desired order of factor levels for plotting consistency
condition_levels <- c("Control", "Experimental", "Ctrl_LHC", "SLE.less.active", "SLE.active")

# Create the initial sample_info data frame
sample_info <- data.frame(
  row.names = sample_names,
  condition = factor(conditions, levels = condition_levels), # Apply factor with correct levels
  batch = factor(ifelse(grepl("Ctrl_LHC|SLE", sample_names), "Public", "Personal"))
)

# <<< NEW: Create a combined group factor for the DESeq2 design >>>
sample_info$group <- factor(paste0(sample_info$batch, "_", sample_info$condition))

# Verify sample names match between metadata and count matrix
if(!all(rownames(sample_info) == colnames(count_matrix))) {
  stop("Mismatch between sample info rownames and count matrix colnames!")
}

# Check for NAs specifically in the condition column
if(any(is.na(sample_info$condition))) {
  na_samples <- rownames(sample_info)[is.na(sample_info$condition)]
  stop("ERROR: NA values were generated in the 'condition' column for samples: ",
       paste(na_samples, collapse=", "),
       ". Check the condition assignment logic in Step 3.")
}
# Check for NAs in the new group column
if(any(is.na(sample_info$group))) {
  na_samples <- rownames(sample_info)[is.na(sample_info$group)]
  stop("ERROR: NA values were generated in the 'group' column. Check metadata creation.")
}


print("Sample Information (colData) - Including 'group':")
print(sample_info)
print("Summary of conditions:")
print(summary(sample_info$condition))
print("Summary of batch:")
print(summary(sample_info$batch))
print("Summary of group (used for DESeq2 design):")
print(summary(sample_info$group))

# --- 4. Create DESeqDataSet Object - REVISED DESIGN ---
# Use the combined 'group' factor in the design to avoid the non-full rank error
# We still keep 'batch' and 'condition' in colData for later use (plotting, batch correction)
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_info,
                              design = ~ group) # Use combined factor

print("DESeqDataSet object created successfully using design = ~ group.")

# --- 5. Pre-filtering (Optional but Recommended) ---
# Remove genes with very low counts across all samples
# keep <- rowSums(counts(dds)) >= 10 # Example threshold
# if(length(keep) > 0) {
#   dds <- dds[keep,]
#   print(paste("Number of genes after filtering:", nrow(dds)))
# } else {
#   print("No genes filtered or pre-filtering skipped.")
# }
# Optional: Filter low count genes (uncomment to use)
keep_genes <- rowSums(counts(dds) >= 10) >= min(table(dds$group)) # Keep if >= 10 counts in at least the smallest group number of samples
dds <- dds[keep_genes,]
print(paste("Number of genes after filtering:", nrow(dds)))


# --- 6. Variance Stabilizing Transformation (VST) ---
# Apply VST. blind=FALSE allows VST to use the design formula (~ group)
# This helps stabilize variance more effectively.
# This step can take a few minutes for large datasets.
print("Starting Variance Stabilizing Transformation (VST)...")
# Using blind=FALSE is generally recommended when you have a design
vsd <- vst(dds, blind = FALSE)
print("VST completed.")

# --- 7. Assess Similarity BEFORE Batch Correction (PCA) ---
# Visualize the data *before* explicit correction to see the batch effect
# Use the original 'condition' and 'batch' for coloring/shaping
print("Generating PCA plot before batch correction...")
pcaData_raw <- plotPCA(vsd, intgroup = c("condition", "batch"), returnData = TRUE)
percentVar_raw <- round(100 * attr(pcaData_raw, "percentVar"))

pca_plot_raw <- ggplot(pcaData_raw, aes(x = PC1, y = PC2, color = condition, shape = batch)) +
  geom_point(size = 3, alpha = 0.8) +
  xlab(paste0("PC1: ", percentVar_raw[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_raw[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA Plot - Before Batch Correction (VST Data)") +
  theme_bw() +
  scale_shape_manual(values=c(Personal=16, Public=17))

# Save or print the plot
# ggsave("PCA_before_correction.png", pca_plot_raw)
print(pca_plot_raw)
print("PCA plot before batch correction generated.")

# --- 8. Explicit Batch Correction using limma::removeBatchEffect ---
# Apply removeBatchEffect to the VST-transformed data
# IMPORTANT: Use the original 'batch' factor here, NOT the 'group' factor
print("Applying limma::removeBatchEffect using the 'batch' variable...")
vst_matrix <- assay(vsd)
batch_vector <- vsd$batch # Get the original batch factor from the object
corrected_vst_matrix <- removeBatchEffect(vst_matrix, batch = batch_vector)
print("Batch correction using removeBatchEffect completed.")

# --- 9. Assess Similarity AFTER Batch Correction (PCA) ---
# Create a new DESeqTransform object with the corrected data for plotting convenience
print("Generating PCA plot after batch correction...")
vsd_corrected <- vsd # Copy structure
assay(vsd_corrected) <- corrected_vst_matrix # Replace data with corrected data

# Use the original 'condition' and 'batch' for coloring/shaping
pcaData_corrected <- plotPCA(vsd_corrected, intgroup = c("condition", "batch"), returnData = TRUE)
percentVar_corrected <- round(100 * attr(pcaData_corrected, "percentVar"))

pca_plot_corrected <- ggplot(pcaData_corrected, aes(x = PC1, y = PC2, color = condition, shape = batch)) +
  geom_point(size = 3, alpha = 0.8) +
  xlab(paste0("PC1: ", percentVar_corrected[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_corrected[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA Plot - After Batch Correction (limma::removeBatchEffect)") +
  theme_bw() +
  scale_shape_manual(values=c(Personal=16, Public=17))

# Save or print the plot
# ggsave("PCA_after_correction.png", pca_plot_corrected)
print(pca_plot_corrected)
print("PCA plot after batch correction generated.")

# --- 10. Assess Similarity AFTER Batch Correction (Heatmap) ---
# Select top N most variable genes based on the corrected data
print("Generating heatmap...")
num_variable_genes <- 500 # Adjust as needed
# Calculate row variances on the corrected matrix
rv <- rowVars(corrected_vst_matrix)
# Ensure num_variable_genes doesn't exceed the number of genes available
num_genes_available <- length(rv)
select <- order(rv, decreasing = TRUE)[seq_len(min(num_variable_genes, num_genes_available))]
heatmap_matrix <- corrected_vst_matrix[select, ]

# Z-score scaling for heatmap visualization (scale rows)
# Add small constant to avoid issues with zero variance rows if any exist (unlikely post-VST)
heatmap_matrix_scaled <- t(scale(t(heatmap_matrix + 1e-6)))

# Create annotation for heatmap columns using original factors
df_annotation <- as.data.frame(colData(vsd_corrected)[, c("condition", "batch")])

# Define colors
num_conditions <- length(levels(df_annotation$condition))
color_palette_conditions <- if (num_conditions <= 8) "Set1" else if (num_conditions <= 12) "Paired" else "Set3"
colors_condition <- colorRampPalette(brewer.pal(min(num_conditions, 8), color_palette_conditions))(num_conditions) # Handle more colors if needed
names(colors_condition) <- levels(df_annotation$condition)
colors_batch <- c(Personal = "lightblue", Public = "darkorange")
anno_colors <- list(
  condition = colors_condition,
  batch = colors_batch # Corrected variable name
)


# Generate the heatmap
heatmap_plot <- pheatmap(heatmap_matrix_scaled,
                         annotation_col = df_annotation,
                         annotation_colors = anno_colors,
                         show_rownames = FALSE,
                         show_colnames = TRUE,
                         cluster_rows = TRUE,
                         cluster_cols = TRUE,
                         main = paste("Heatmap of Top", length(select), "Variable Genes (After Batch Correction)"),
                         fontsize_col = 8,
                         clustering_distance_cols = "euclidean",
                         clustering_method = "ward.D2")

# Save or print the plot
# png("Heatmap_after_correction.png", width=10, height=8, units="in", res=300)
# print(heatmap_plot)
# dev.off()
print(heatmap_plot)
print("Heatmap generated.")


# --- 11. Interpretation ---
cat("\n--- Interpretation Guidance --- \n")
cat("1. Examine 'PCA Plot - Before Batch Correction': Samples likely cluster strongly by 'batch' (shape).\n")
cat("2. Examine 'PCA Plot - After Batch Correction': Focus on the proximity of 'Experimental' (Personal batch) and 'SLE.active' (Public batch) points (colors) relative to each other and other conditions, now that the major batch effect is reduced.\n")
cat("3. Examine the 'Heatmap': Look at the column dendrogram (clustering tree above the heatmap). Do the columns corresponding to 'Experimental' and 'SLE.active' samples cluster closely together, suggesting similar expression profiles across the most variable genes after correction?\n")
cat("4. Conclusion: Based on the visual proximity in the corrected PCA and the clustering in the heatmap, assess the degree of resemblance. Are 'Experimental' and 'SLE.active' highly similar (clustering tightly together, separate from others), somewhat similar (clustering nearer to each other than to distant groups), or still distinct despite batch correction?\n")
cat("\nAnalysis complete.\n")