# --- 0. Install Packages (Run these lines in the R console if needed) ---

# install.packages("readr")
# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("stringr")
# install.packages("tidyr") # Often used with dplyr/ggplot2, good to have
#
# # biomaRt requires BiocManager for installation
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("biomaRt")

# --- 1. Load Required Libraries ---
library(readr)
library(dplyr)
library(biomaRt)
library(ggplot2)
library(stringr)
library(tidyr) # For complete() function

# --- 2. Define File Paths ---
# !!! IMPORTANT: Adjust these paths to match your file locations !!!
# Example assumes your script is in a 'scripts' folder and data is in 'public' sibling folder
# Use absolute paths if necessary: e.g., "/path/to/your/project/public/UpDEG.csv"
up_file <- "../public/UpDEG.csv"
down_file <- "../public/DownDEG.csv"
output_plot_file <- "../public/chromosome_distribution_plot.jpg"
output_up_chr_file <- "../public/UpDEG_chr.csv"
output_down_chr_file <- "../public/DownDEG_chr.csv"
output_dir <- "../public/"

# --- 2b. Create Output Directory ---
# Create output directory if it doesn't exist (do this early)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message("Created output directory: ", output_dir)
}

# --- 3. Import DEG Files ---
message("Importing DEG files...")
# Check if files exist
if (!file.exists(up_file)) stop("UpDEG file not found at: ", up_file, ". Please check the path.")
if (!file.exists(down_file)) stop("DownDEG file not found at: ", down_file, ". Please check the path.")

up_deg <- read_csv(up_file, show_col_types = FALSE)
down_deg <- read_csv(down_file, show_col_types = FALSE)

# --- 4. Prepare Gene IDs and Get Chromosome Annotations ---
message("Fetching gene annotations from Ensembl using biomaRt...")

# Remove version suffix from Ensembl IDs for matching with biomaRt
up_deg <- up_deg %>% mutate(ensembl_gene_id_base = str_remove(gene_id, "\\..*"))
down_deg <- down_deg %>% mutate(ensembl_gene_id_base = str_remove(gene_id, "\\..*"))

# Combine unique gene IDs from both files
all_gene_ids <- unique(c(up_deg$ensembl_gene_id_base, down_deg$ensembl_gene_id_base))

# Connect to Ensembl (uses the current version by default)
# Using https is recommended
# Retry mechanism basic example (could be more sophisticated)
ensembl_mart <- NULL
attempt <- 1
max_attempts <- 3
while(is.null(ensembl_mart) && attempt <= max_attempts) {
  try({
    # Use host = "https://www.ensembl.org" for the main site
    # Or use a mirror like "https://useast.ensembl.org", "https://uswest.ensembl.org", etc.
    ensembl_mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "https://www.ensembl.org")
    # Explicitly select the human dataset immediately after connecting
    human_dataset <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl_mart)
  }, silent = TRUE) # silent=TRUE suppresses the default error message
  if (is.null(ensembl_mart)) {
    message("Attempt ", attempt, " to connect to Ensembl failed. Retrying in 5 seconds...")
    Sys.sleep(5)
    attempt <- attempt + 1
  }
}
if (is.null(ensembl_mart)) {
  stop("Failed to connect to Ensembl after ", max_attempts, " attempts. Check connection or Ensembl status (www.ensembl.org). Try using a mirror host like 'https://useast.ensembl.org'.")
} else {
  message("Successfully connected to Ensembl mart and selected human dataset.")
}

# Check if human_dataset was successfully assigned
if (!exists("human_dataset") || is.null(human_dataset)) {
  stop("Failed to select the 'hsapiens_gene_ensembl' dataset.")
}

# Fetch chromosome information for the genes
# Handle potential errors during biomaRt query
gene_annotations <- tryCatch({
  getBM(
    attributes = c("ensembl_gene_id", "chromosome_name"),
    filters    = "ensembl_gene_id",
    values     = all_gene_ids,
    mart       = human_dataset # Use the specific dataset object
  )
}, error = function(e) {
  stop("Error fetching data from biomaRt: ", e$message,
       "\nCheck your internet connection and if Ensembl services are available.")
})

message("Annotation fetching complete. ", nrow(gene_annotations), " annotations retrieved.")

# --- 5. Merge Annotations and Count Genes per Chromosome ---
message("Processing gene counts per chromosome...")

# Merge annotations with DEG data
# Using 'ensembl_gene_id_base' from DEG files and 'ensembl_gene_id' from annotations
up_deg_annotated <- left_join(up_deg, gene_annotations, by = c("ensembl_gene_id_base" = "ensembl_gene_id"))
down_deg_annotated <- left_join(down_deg, gene_annotations, by = c("ensembl_gene_id_base" = "ensembl_gene_id"))

# --- 5b. Export Annotated DEG Lists ---
message("Exporting annotated DEG lists with chromosome information...")

# Select desired columns and arrange them for clarity
# Keep original columns + chromosome_name, remove the temporary base ID column
# Use dplyr::select to avoid namespace conflicts
up_deg_to_save <- up_deg_annotated %>%
  dplyr::select(gene_id, gene_symbol, chromosome_name, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)

down_deg_to_save <- down_deg_annotated %>%
  dplyr::select(gene_id, gene_symbol, chromosome_name, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)

# Write the annotated data frames to new CSV files
write_csv(up_deg_to_save, output_up_chr_file)
write_csv(down_deg_to_save, output_down_chr_file)

message("Annotated DEG lists saved to:")
message("  ", output_up_chr_file)
message("  ", output_down_chr_file)

# --- Processing for Plotting ---

# Define standard chromosomes (adjust if you need others like MT)
valid_chromosomes <- c(as.character(1:22), "X", "Y")
# Ensure correct order, making it a factor *before* using complete
chromosome_order <- factor(valid_chromosomes, levels = valid_chromosomes)

# Ensure 'regulation' is a factor with levels in the desired order for plotting/legend
regulation_levels <- c("Up", "Down")

# Count Up-regulated genes per valid chromosome
up_counts <- up_deg_annotated %>%
  filter(chromosome_name %in% valid_chromosomes) %>%
  mutate(chromosome_name = factor(chromosome_name, levels = levels(chromosome_order))) %>% # Ensure factor levels before count
  count(chromosome_name, name = "count", .drop = FALSE) %>% # use .drop=FALSE to keep all chromosomes
  mutate(regulation = factor("Up", levels = regulation_levels)) # Set regulation as factor

# Count Down-regulated genes per valid chromosome
# Make counts negative for plotting downwards
down_counts <- down_deg_annotated %>%
  filter(chromosome_name %in% valid_chromosomes) %>%
  mutate(chromosome_name = factor(chromosome_name, levels = levels(chromosome_order))) %>% # Ensure factor levels before count
  count(chromosome_name, name = "count", .drop = FALSE) %>% # use .drop=FALSE to keep all chromosomes
  mutate(regulation = factor("Down", levels = regulation_levels), # Set regulation as factor
         count = -count) # Negative count

# Combine counts
all_counts <- bind_rows(up_counts, down_counts)

# Ensure all chromosomes are present for both Up/Down, filling missing with 0 counts
# The chromosome_name factor is already set correctly from the previous steps
# The regulation factor levels are also set
all_counts <- all_counts %>%
  complete(chromosome_name, regulation, fill = list(count = 0)) %>%
  # Re-apply factor levels just in case they were lost (shouldn't be with factors set beforehand)
  mutate(chromosome_name = factor(chromosome_name, levels = levels(chromosome_order)),
         regulation = factor(regulation, levels = regulation_levels)) %>%
  filter(!is.na(regulation)) # Remove any rows created with NA regulation (shouldn't happen here)


# --- 6. Create Chromosome Distribution Plot ---
message("Creating the chromosome distribution plot...")

# Find the maximum absolute count for symmetrical y-axis scaling
max_abs_count <- max(abs(all_counts$count), na.rm = TRUE)
# Handle case where there might be no counts at all
y_limit <- if (max_abs_count > 0) ceiling(max_abs_count * 1.1) else 1 # Add 10% padding, minimum limit 1

# Create the plot using ggplot2
chromosome_plot <- ggplot(all_counts, aes(x = chromosome_name, y = count, fill = regulation)) +
  geom_col(position = "identity", width = 0.7) + # Use geom_col as counts are pre-calculated
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) + # Use linewidth instead of size for lines
  scale_fill_manual(
    name = NULL, # No legend title
    values = c("Up" = "red", "Down" = "blue"), # Assigns colors to factor levels
    # *** CHANGE HERE: Use a named vector for labels to ensure correct mapping ***
    labels = c("Up" = "UpDEGs", "Down" = "DownDEGs"),
    # Ensure breaks match the factor levels if needed (usually automatic with factors)
    breaks = c("Up", "Down")
  ) +
  scale_y_continuous(
    limits = c(-y_limit, y_limit), # Symmetrical limits
    breaks = scales::pretty_breaks(n = 5), # Auto-calculate nice breaks
    labels = abs # Show absolute values on y-axis labels
  ) +
  labs(
    # title = "Chromosome Distribution of DEGs", # Optional title
    x = "Chromosome",
    y = "Gene Count"
  ) +
  theme_classic() + # A clean theme
  theme(
    # plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 12),
    axis.text.x = element_text(angle = 0, hjust = 0.5), # Ensure chromosome labels are horizontal
    legend.position = "top", # Position legend at the top
    legend.justification = "center", # Center the legend
    legend.text = element_text(size=10)
  )

# --- 7. Export the Plot ---
message("Exporting the plot to: ", output_plot_file)

# Save the plot as JPEG (Directory was already created earlier)
ggsave(
  filename = output_plot_file,
  plot = chromosome_plot,
  device = "jpeg",
  width = 8,       # Inches
  height = 5,      # Inches
  dpi = 300        # Resolution
)

message("Script finished successfully.")