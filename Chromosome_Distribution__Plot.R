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
up_file <- "../public/UpDEG.csv"
down_file <- "../public/DownDEG.csv"
output_plot_file <- "../public/chromosome_distribution_plot.jpg"

# --- 3. Import DEG Files ---
message("Importing DEG files...")
# Check if files exist
if (!file.exists(up_file)) stop("UpDEG file not found at: ", up_file)
if (!file.exists(down_file)) stop("DownDEG file not found at: ", down_file)

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
ensembl_mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "https://www.ensembl.org")
human_dataset <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl_mart)

# Fetch chromosome information for the genes
# Handle potential errors during biomaRt query
gene_annotations <- tryCatch({
  getBM(
    attributes = c("ensembl_gene_id", "chromosome_name"),
    filters    = "ensembl_gene_id",
    values     = all_gene_ids,
    mart       = human_dataset
  )
}, error = function(e) {
  stop("Error fetching data from biomaRt: ", e$message,
       "\nCheck your internet connection and if Ensembl services are available.")
})

message("Annotation fetching complete. ", nrow(gene_annotations), " annotations retrieved.")

# --- 5. Merge Annotations and Count Genes per Chromosome ---
message("Processing gene counts per chromosome...")

# Merge annotations with DEG data
up_deg_annotated <- left_join(up_deg, gene_annotations, by = c("ensembl_gene_id_base" = "ensembl_gene_id"))
down_deg_annotated <- left_join(down_deg, gene_annotations, by = c("ensembl_gene_id_base" = "ensembl_gene_id"))

# Define standard chromosomes (adjust if you need others like MT)
valid_chromosomes <- c(as.character(1:22), "X", "Y")
chromosome_order <- factor(valid_chromosomes, levels = valid_chromosomes) # Ensure correct order

# Count Up-regulated genes per valid chromosome
up_counts <- up_deg_annotated %>%
  filter(chromosome_name %in% valid_chromosomes) %>%
  count(chromosome_name, name = "count") %>%
  mutate(regulation = "Up")

# Count Down-regulated genes per valid chromosome
# Make counts negative for plotting downwards
down_counts <- down_deg_annotated %>%
  filter(chromosome_name %in% valid_chromosomes) %>%
  count(chromosome_name, name = "count") %>%
  mutate(regulation = "Down", count = -count) # Negative count

# Combine counts
all_counts <- bind_rows(up_counts, down_counts)

# Ensure all chromosomes are present for both Up/Down, filling missing with 0 counts
# Convert chromosome_name to a factor with the specified order first
all_counts <- all_counts %>%
  mutate(chromosome_name = factor(chromosome_name, levels = levels(chromosome_order))) %>%
  complete(chromosome_name, regulation, fill = list(count = 0)) %>%
  filter(!is.na(regulation)) # Remove any rows created with NA regulation (shouldn't happen here)

# --- 6. Create Chromosome Distribution Plot ---
message("Creating the chromosome distribution plot...")

# Find the maximum absolute count for symmetrical y-axis scaling
max_abs_count <- max(abs(all_counts$count), na.rm = TRUE)
y_limit <- ceiling(max_abs_count * 1.1) # Add 10% padding

# Create the plot using ggplot2
chromosome_plot <- ggplot(all_counts, aes(x = chromosome_name, y = count, fill = regulation)) +
  geom_col(position = "identity", width = 0.7) + # Use geom_col as counts are pre-calculated
  geom_hline(yintercept = 0, color = "black", size = 0.5) + # Add line at y=0
  scale_fill_manual(
    name = NULL, # No legend title
    values = c("Up" = "red", "Down" = "blue"),
    labels = c("UpDEGs", "DownDEGs") # Legend labels matching example
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

# Create output directory if it doesn't exist
output_dir <- dirname(output_plot_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message("Created output directory: ", output_dir)
}

# Save the plot as JPEG
ggsave(
  filename = output_plot_file,
  plot = chromosome_plot,
  device = "jpeg",
  width = 8,    # Inches
  height = 5,   # Inches
  dpi = 300     # Resolution
)

message("Script finished successfully.")