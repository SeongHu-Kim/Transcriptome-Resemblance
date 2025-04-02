# --- Existing Code Starts ---

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
  # Handle potential errors during file reading
  tryCatch({
    read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)$gene_id
  }, error = function(e) {
    warning("Error reading file: ", file, " - ", e$message)
    return(NULL) # Return NULL if a file cannot be read
  })
})

# Remove NULL elements from the list if any file failed to read
gene_id_lists <- Filter(Negate(is.null), gene_id_lists)

# Check if we have lists to compare and if they are all identical
all_identical <- FALSE
if (length(gene_id_lists) == length(file_paths) && length(gene_id_lists) > 0) {
  all_identical <- all(sapply(gene_id_lists[-1], identical, gene_id_lists[[1]]))
} else if (length(gene_id_lists) < length(file_paths)) {
  cat("Warning: Not all files could be read successfully. Cannot guarantee gene_id identity.\n")
} else {
  cat("Warning: No gene_id lists were successfully read.\n")
}


# Print results
if (all_identical) {
  cat("All", length(file_paths), "files have identical rows for the gene_id column\n")
} else {
  cat("The", length(file_paths), "files do not have identical rows for the gene_id column, or some files could not be read.\n")
  # Optional: Add more detailed reporting here if needed
  # For example, which files differ or which failed to read.
}

# If gene_id columns are identical, create and export GeneCountMatrix_Original.csv
if (all_identical) {
  
  # --- New Code Starts Here ---
  
  # Load necessary Bioconductor annotation libraries
  # Ensure BiocManager, AnnotationDbi, and org.Hs.eg.db are installed
  # if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  # BiocManager::install("AnnotationDbi")
  # BiocManager::install("org.Hs.eg.db") # Use the correct organism package if not Human
  
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) || !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("Please install required Bioconductor packages: BiocManager::install(c('AnnotationDbi', 'org.Hs.eg.db'))")
  }
  library(AnnotationDbi)
  library(org.Hs.eg.db) # Assumes Human data (Homo sapiens)
  
  # Get the list of gene IDs (assuming they are all identical, use the first list)
  gene_ids_with_version <- gene_id_lists[[1]]
  
  # Remove version numbers from Ensembl IDs for mapping (e.g., ".16")
  gene_ids_no_version <- gsub("\\.[0-9]+$", "", gene_ids_with_version)
  
  # Map Ensembl IDs (without version) to Gene Symbols
  # Handle potential multiple mappings or unmapped IDs
  cat("Mapping Ensembl IDs to Gene Symbols...\n")
  gene_symbols <- mapIds(org.Hs.eg.db,
                         keys = gene_ids_no_version,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first") # Takes the first symbol if multiple exist
  
  cat("Mapping complete.\n")
  
  # Create the initial data frame with gene_id (original with version) and gene_symbol
  gene_count_matrix <- data.frame(
    gene_id = gene_ids_with_version,
    gene_symbol = gene_symbols,
    stringsAsFactors = FALSE # Good practice
  )
  
  # --- End of New Code ---
  
  
  # --- Existing Code Resumes (Modified Slightly) ---
  
  # Define sample labels (remains the same)
  sample_labels <- c("Control_R1", "Control_R2", "Control_R3", "Control_R4",
                     "Experimental_R1", "Experimental_R2", "Experimental_R3", "Experimental_R4")
  
  # Loop through files to add count data (remains the same logic)
  cat("Reading count data...\n")
  for (i in seq_along(file_paths)) {
    # Read the expected_count column for the current file
    counts <- tryCatch({
      read.table(file_paths[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)$expected_count
    }, error = function(e) {
      warning("Error reading counts from file: ", file_paths[i], " - ", e$message)
      return(rep(NA, nrow(gene_count_matrix))) # Return NAs if counts can't be read
    })
    
    # Add the counts as a new column to the data frame
    # Ensure the sample label is valid before assigning
    if (nchar(sample_labels[i]) > 0) {
      gene_count_matrix[[sample_labels[i]]] <- counts
    } else {
      warning("Invalid sample label generated for index: ", i)
    }
  }
  cat("Count data added.\n")
  
  # Write the combined data frame to a CSV file (remains the same)
  output_file <- "../public/GeneCountMatrix_Original.csv"
  cat("Writing output file:", output_file, "\n")
  write.csv(gene_count_matrix, file = output_file, row.names = FALSE, quote = TRUE) # Use quote=TRUE for safety
  cat("Successfully created:", output_file, "\n")
  
} else {
  cat("GeneCountMatrix_Original.csv was not created because gene_id columns were not identical across all files or file reading errors occurred.\n")
}

# --- Existing Code Ends ---