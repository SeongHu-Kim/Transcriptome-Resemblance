# --- Load Required Libraries ---
# Load DESeq2 for differential expression analysis
# Load tidyverse for data manipulation (specifically dplyr and readr)
library(DESeq2)
library(tidyverse)

# --- Set File Paths ---
# Define the input file path for the gene count matrix
input_file <- "../public/GeneCountMatrix_filtered.csv"
# Define the output file path for the full DESeq2 results
output_file <- "../public/DESeq2_results.csv"
# Define output file paths for filtered differentially expressed genes (DEGs)
up_deg_file <- "../public/UpDEG.csv"
down_deg_file <- "../public/DownDEG.csv"

# --- 1. Import Data ---
# Read the gene count matrix from the CSV file
# 'check.names=FALSE' prevents R from modifying column names (e.g., replacing '-' with '.')
count_data_raw <- read.csv(input_file, header = TRUE, row.names = "gene_id", check.names = FALSE)

# Display the first few rows and dimensions to verify import
print("Raw count data dimensions:")
print(dim(count_data_raw))
print("First few rows of raw count data:")
print(head(count_data_raw))

# --- 2. Prepare Data for DESeq2 ---

# Separate gene symbols (if needed later, though not directly used by DESeq2)
gene_symbols <- count_data_raw %>% select(gene_symbol)
count_matrix_raw <- count_data_raw %>% select(-gene_symbol) # Keep only count columns

# DESeq2 requires integer counts. Round the count data to the nearest integer.
count_matrix <- round(count_matrix_raw)

# Display the first few rows of the rounded count matrix
print("First few rows of rounded count matrix:")
print(head(count_matrix))

# Create the metadata (colData) describing the samples
# Extract sample names from the column names of the count matrix
sample_names <- colnames(count_matrix)

# Determine the condition (e.g., 'Control' or 'Experimental') for each sample
# This assumes column names clearly indicate the group (e.g., start with 'Control_' or 'Experimental_')
conditions <- sub("_R\\d+", "", sample_names) # Remove replicate number (e.g., _R1)
colData <- data.frame(row.names = sample_names, condition = factor(conditions))

# Ensure the order of rows in colData matches the order of columns in count_matrix
if (!all(rownames(colData) == colnames(count_matrix))) {
  stop("Mismatch between column names of count matrix and row names of metadata!")
}

print("Sample Metadata (colData):")
print(colData)

# --- 3. Create DESeqDataSet Object ---
# Construct the DESeqDataSet object from the count matrix and metadata
# The design formula (~ condition) tells DESeq2 to model counts based on the 'condition' variable
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = colData,
                              design = ~ condition)

# --- 4. Run DESeq2 Analysis ---
# This function performs normalization, dispersion estimation, and statistical testing
print("Running DESeq2 analysis...")
dds <- DESeq(dds)
print("DESeq2 analysis complete.")

# --- 5. Extract and Format Results ---
# Extract results for the comparison between 'Experimental' and 'Control'
# 'contrast' specifies the comparison: condition, numerator level, denominator level
# By default, levels are compared alphabetically, so 'Experimental' vs 'Control'
res <- results(dds, contrast=c("condition", "Experimental", "Control"))

# Order results by adjusted p-value (padj)
resOrdered <- res[order(res$padj),]

# Convert the results object to a data frame
resOrdered_df <- as.data.frame(resOrdered)

# Add gene IDs (which are the rownames) as a column
resOrdered_df <- resOrdered_df %>%
  rownames_to_column(var = "gene_id")

# Merge with the original gene symbols
# Ensure gene_symbols has 'gene_id' as a column for merging
gene_symbols_df <- gene_symbols %>%
  rownames_to_column(var = "gene_id")

# Perform the merge
results_final <- merge(gene_symbols_df, resOrdered_df, by = "gene_id", all.y = TRUE)

# Reorder columns for clarity
results_final <- results_final %>%
  select(gene_id, gene_symbol, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)

# Display the first few rows of the final results table
print("First few rows of DESeq2 results:")
print(head(results_final))

# --- 6. Output Full Results ---
# Save the final results table to a CSV file
write.csv(results_final, file = output_file, row.names = FALSE)

print(paste("DESeq2 full results successfully saved to:", output_file))

# --- 7. Filter and Save Upregulated DEGs ---
print("Filtering for upregulated genes...")
# Define thresholds
lfc_threshold_up <- 1 # log2(2)
padj_threshold <- 0.05

# Filter for upregulated genes (log2FC > threshold and padj < threshold)
# Also remove rows where padj or log2FoldChange is NA to avoid errors in sorting/filtering
up_degs <- results_final %>%
  filter(!is.na(padj) & !is.na(log2FoldChange)) %>% # Remove NAs first
  filter(log2FoldChange > lfc_threshold_up & padj < padj_threshold) %>%
  mutate(abs_log2FoldChange = abs(log2FoldChange)) %>% # Add absolute LFC for sorting
  arrange(desc(abs_log2FoldChange)) %>% # Sort by absolute LFC descending
  select(-abs_log2FoldChange) # Remove the temporary column

# Save the upregulated DEGs to a CSV file
write.csv(up_degs, file = up_deg_file, row.names = FALSE)

print(paste("Upregulated DEGs (log2FC >", round(lfc_threshold_up, 3), ", padj <", padj_threshold, ") saved to:", up_deg_file))
print(paste("Number of upregulated DEGs:", nrow(up_degs)))

# --- 8. Filter and Save Downregulated DEGs ---
print("Filtering for downregulated genes...")
# Define threshold (negative for downregulation)
lfc_threshold_down <- -1 # log2(1/2)

# Filter for downregulated genes (log2FC < threshold and padj < threshold)
# Also remove rows where padj or log2FoldChange is NA
down_degs <- results_final %>%
  filter(!is.na(padj) & !is.na(log2FoldChange)) %>% # Remove NAs first
  filter(log2FoldChange < lfc_threshold_down & padj < padj_threshold) %>%
  mutate(abs_log2FoldChange = abs(log2FoldChange)) %>% # Add absolute LFC for sorting
  arrange(desc(abs_log2FoldChange)) %>% # Sort by absolute LFC descending
  select(-abs_log2FoldChange) # Remove the temporary column

# Save the downregulated DEGs to a CSV file
write.csv(down_degs, file = down_deg_file, row.names = FALSE)

print(paste("Downregulated DEGs (log2FC <", round(lfc_threshold_down, 3), ", padj <", padj_threshold, ") saved to:", down_deg_file))
print(paste("Number of downregulated DEGs:", nrow(down_degs)))

# --- Optional: Save Normalized Counts ---
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts_df <- as.data.frame(normalized_counts) %>%
 rownames_to_column(var = "gene_id") %>%
 merge(gene_symbols_df, by = "gene_id") %>%
 select(gene_id, gene_symbol, everything()) # Reorder columns
write.csv(normalized_counts_df, file = "../public/DESeq2_normalized_counts.csv", row.names = FALSE)
print("Normalized counts saved to ../public/DESeq2_normalized_counts.csv")

