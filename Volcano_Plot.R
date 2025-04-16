# Install necessary packages if you haven't already
# install.packages("ggplot2")
# install.packages("dplyr")

# Load required libraries
library(ggplot2)
library(dplyr)

# --- Configuration ---
# Define input file path
input_file <- "../public/DESeq2_results.csv"
# Define output file path
output_file <- "../public/volcano_plot.jpg"
# Define thresholds
padj_threshold <- 0.05
# log2FoldChange threshold (equivalent to fold change of 1.5)
lfc_threshold <- 0.58496

# --- Data Loading and Preparation ---

# Check if the input file exists
if (!file.exists(input_file)) {
  stop(paste("Error: Input file not found at", input_file))
}

# Read the DESeq2 results CSV file
# Using tryCatch to handle potential errors during file reading
deseq_results <- tryCatch({
  read.csv(input_file)
}, error = function(e) {
  stop(paste("Error reading CSV file:", e$message))
})

# Check if essential columns exist
required_cols <- c("log2FoldChange", "padj")
if (!all(required_cols %in% names(deseq_results))) {
  stop(paste("Error: Input CSV must contain columns:", paste(required_cols, collapse=", ")))
}

# Add a column to classify genes based on significance and fold change
# Handle potential NA values in padj - replace NA with 1 (least significant)
# Handle potential NA values in log2FoldChange - replace NA with 0
deseq_results <- deseq_results %>%
  mutate(
    padj = ifelse(is.na(padj), 1, padj), # Replace NA padj with 1
    log2FoldChange = ifelse(is.na(log2FoldChange), 0, log2FoldChange), # Replace NA log2FC with 0
    gene_status = case_when(
      padj < padj_threshold & log2FoldChange > lfc_threshold  ~ "Upregulated",
      padj < padj_threshold & log2FoldChange < -lfc_threshold ~ "Downregulated",
      TRUE                                                    ~ "Not Significant" # Default case
    ),
    # Ensure gene_status is a factor for consistent coloring
    gene_status = factor(gene_status, levels = c("Upregulated", "Downregulated", "Not Significant"))
  )

# --- Volcano Plot Creation ---

# Create the volcano plot using ggplot2
volcano_plot <- ggplot(deseq_results, aes(x = log2FoldChange, y = -log10(padj), color = gene_status)) +
  geom_point(alpha = 0.6, size = 1.5) + # Add points with some transparency
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "Not Significant" = "grey")) + # Define colors
  theme_minimal(base_size = 14) + # Use a minimal theme
  labs(
    title = "Volcano Plot of Differential Gene Expression",
    x = expression(Log[2]~"Fold Change"), # Label x-axis
    y = expression(-Log[10]~"Adjusted p-value"), # Label y-axis
    color = "Gene Status" # Legend title
  ) +
  # Add horizontal line for padj threshold
  geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "black", linewidth = 0.5) +
  # Add vertical lines for log2FoldChange thresholds
  geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", color = "black", linewidth = 0.5) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), # Center and bold title
    legend.position = "bottom", # Place legend at the bottom
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5) # Add border
  )

# Print the plot to the R graphics device (optional, good for interactive sessions)
# print(volcano_plot)

# --- Save Plot ---

# Save the plot to the specified file path
# Using tryCatch to handle potential errors during saving
tryCatch({
  ggsave(output_file, plot = volcano_plot, width = 8, height = 6, dpi = 300, device = "jpeg")
  print(paste("Volcano plot successfully saved to:", output_file))
}, error = function(e) {
  warning(paste("Warning: Could not save the plot.", e$message))
  # Attempt to save in the current working directory if ../public fails
  alt_output_file <- "volcano_plot.jpg"
  tryCatch({
    ggsave(alt_output_file, plot = volcano_plot, width = 8, height = 6, dpi = 300, device = "jpeg")
    print(paste("Plot saved to current working directory instead:", alt_output_file))
  }, error = function(e_alt) {
    stop(paste("Error: Failed to save plot to both specified and current directory.", e_alt$message))
  })
})

# Clean up variables (optional)
# rm(deseq_results, volcano_plot, input_file, output_file, padj_threshold, lfc_threshold, required_cols)
