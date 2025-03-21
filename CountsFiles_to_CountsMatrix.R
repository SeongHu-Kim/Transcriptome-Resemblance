# Define file paths
file_paths <- c(
  "../public/Control_R1.genes.results",
  "../public/Control_R2.genes.results",
  "../public/Control_R3.genes.results",
  "../public/Control_R4.genes.results",
  "../public/Experimental_R1.genes.results",
  "../public/Experimental_R2.genes.results",
  "../public/Experimental_R3.genes.results",
  "../public/Experimental_R4.genes.results"
)

# Read "gene_id" column from all files
gene_id_lists <- lapply(file_paths, function(file) {
  read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)$gene_id
})

# Check if all "gene_id" columns are identical
all_identical <- all(sapply(gene_id_lists, identical, gene_id_lists[[1]]))

# Print results
if (all_identical) {
  cat("All 8 files have identical rows for the gene_id column\n")
} else {
  cat("The 8 files does not have identical rows for the gene_id column\n")
}

# If gene_id columns are identical, create and export GeneCountMatrix_Original.csv
if (all_identical) {
  gene_count_matrix <- data.frame(gene_id = gene_id_lists[[1]])
  sample_labels <- c("Control_R1", "Control_R2", "Control_R3", "Control_R4", 
                     "Experimental_R1", "Experimental_R2", "Experimental_R3", "Experimental_R4")
  
  for (i in seq_along(file_paths)) {
    counts <- read.table(file_paths[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)$expected_count
    gene_count_matrix[[sample_labels[i]]] <- counts
  }
  
  write.csv(gene_count_matrix, file = "../public/GeneCountMatrix_Original.csv", row.names = FALSE)
}
