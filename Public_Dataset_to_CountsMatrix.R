# --- Configuration ---

# Base directory containing the subdirectories
base_dir <- "../public/GSE97263" 

# Output file path and name
output_file <- "../public/GeneCountMatrix_GSE97263.csv"

# --- Step 0: Define the Sample Name Mapping ---

# Create a data frame for mapping GSM IDs to desired sample names
# Sourced from the provided list in the request
mapping_data <- data.frame(
  GSM_ID = c(
    "GSM2560285", "GSM2560286", "GSM2560287", "GSM2560288", "GSM2560289",
    "GSM2560290", "GSM2560291", "GSM2560292", "GSM2560293", "GSM2560294",
    "GSM2560295", "GSM2560296", "GSM2560297", "GSM2560298", "GSM2560299",
    "GSM2560300", "GSM2560301", "GSM2560302", "GSM2560303", "GSM2560304",
    "GSM2560305", "GSM2560306", "GSM2560307", "GSM2560308", "GSM2560309",
    "GSM2560310", "GSM2560311", "GSM2560312", "GSM2560313", "GSM2560314",
    "GSM2560315", "GSM2560316", "GSM2560317", "GSM2560318", "GSM2560319",
    "GSM2560320", "GSM2560321", "GSM2560322", "GSM2560323", "GSM2560324",
    "GSM2560325", "GSM2560326", "GSM2560327", "GSM2560328"
  ),
  SampleName = c(
    "SLE active_CDR-2", "SLE active_CDR-6", "SLE active_CDR-20", "SLE active_CDR-21", 
    "SLE active_CDR-24", "SLE active_CDR-25", "SLE active_CDR-27", "SLE active_CDR-35",
    "SLE active_CDR-36", "SLE active_CDR-60", "SLE active_CDR-65", "SLE active_CDR-67", 
    "SLE active_CDR-69", "SLE active_CDR-23", "SLE less active_CDR-29", "SLE less active_CDR-37",
    "SLE less active_CDR-41", "SLE less active_CDR-45", "SLE less active_CDR-46", "SLE less active_CDR-50",
    "SLE less active_CDR-51", "SLE less active_CDR-52", "SLE less active_CDR-53", "SLE less active_CDR-54", 
    "SLE less active_CDR-55", "SLE less active_CDR-56", "SLE less active_CDR-57", "SLE less active_CDR-59",
    "SLE less active_CDR-61", "SLE less active_CDR-64", "Ctrl_LHC-1", "Ctrl_LHC-4", 
    "Ctrl_LHC-5", "Ctrl_LHC-6", "Ctrl_LHC-7", "Ctrl_LHC-8", "Ctrl_LHC-9", 
    "Ctrl_LHC-12", "Ctrl_LHC-13", "Ctrl_LHC-14", "Ctrl_LHC-15", "Ctrl_LHC-16", 
    "Ctrl_LHC-17", "Ctrl_LHC-19"
  ),
  stringsAsFactors = FALSE # Important to keep names as characters
)

# --- Step 1: Find and Import Files ---

# List the subdirectories (which have .txt extension in their name)
# full.names = TRUE gives the full path
# recursive = FALSE ensures we only look at the top level within base_dir
subdirs <- list.dirs(path = base_dir, full.names = TRUE, recursive = FALSE)

# Filter out any unexpected directories/files if necessary (e.g., if base_dir contained other items)
# This pattern assumes directories start with GSM and end with .txt
subdirs <- subdirs[grepl("GSM.*\\.txt$", basename(subdirs))]

# Check if any directories were found
if (length(subdirs) == 0) {
  stop("No subdirectories matching the expected pattern found in: ", base_dir)
}

# Construct the full path to the actual text file within each subdirectory
# The text file name is identical to the subdirectory name
file_paths <- file.path(subdirs, basename(subdirs))

# Verify that the files exist
existing_files <- file.exists(file_paths)
if (!all(existing_files)) {
  warning("Some expected files are missing: \n", paste(file_paths[!existing_files], collapse = "\n"))
  # Filter to keep only paths to files that actually exist
  file_paths <- file_paths[existing_files] 
  subdirs <- subdirs[existing_files] # Keep subdirs and file_paths in sync
  if (length(file_paths) == 0) {
    stop("No valid files found to import.")
  }
}

# Extract GSM IDs from the subdirectory names to use for mapping later
# This assumes the GSM ID is the first part before the first underscore
# Example: GSM2560285_IGF0005584.count.RF.txt -> GSM2560285
gsm_ids_from_files <- sub("\\_.*", "", basename(subdirs))

# --- Step 2 & 3: Read data, Rename Columns (implicitly), and Combine ---

# Read the first file to get the gene IDs (assuming all files have the same genes in the same order)
if (length(file_paths) > 0) {
  message("Reading gene IDs from: ", file_paths[1])
  first_file_data <- read.table(file_paths[1], sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  # Check if the file has at least two columns
  if (ncol(first_file_data) < 2) {
    stop("File ", file_paths[1], " does not have at least two columns.")
  }
  gene_ids <- first_file_data[, 1] 
} else {
  stop("Cannot read gene IDs, no files to process.")
}

# Read the count data (second column) from all files
# Use lapply to iterate through file paths and read each one
# We only need the second column (counts)
message("Reading count data from ", length(file_paths), " files...")
all_counts_list <- lapply(file_paths, function(file) {
  message("Reading: ", file)
  tryCatch({
    # Read file, expecting 2 tab-separated columns, no header
    data <- read.table(file, sep = "\t", header = FALSE, colClasses = c(NA, "numeric")) 
    # Check if the file has exactly two columns
    if (ncol(data) != 2) {
      warning("File ", file, " does not have exactly two columns. Skipping.")
      return(NULL) # Return NULL if format is wrong
    }
    # Return only the second column (counts)
    return(data[, 2])
  }, error = function(e) {
    warning("Error reading file ", file, ": ", e$message, ". Skipping.")
    return(NULL) # Return NULL in case of read errors
  })
})

# Remove any NULL entries from the list (resulting from read errors or incorrect columns)
original_length <- length(all_counts_list)
all_counts_list <- all_counts_list[!sapply(all_counts_list, is.null)]
if (length(all_counts_list) < original_length) {
  warning(original_length - length(all_counts_list), " files were skipped due to errors or incorrect format.")
  # Adjust gsm_ids_from_files to match the successfully read files
  # This assumes the order is preserved and NULLs correspond to skipped files
  gsm_ids_from_files <- gsm_ids_from_files[!sapply(all_counts_list, is.null)] # Re-check, this line might be redundant if list is already filtered
  # Need to find which *indices* were kept to filter gsm_ids_from_files correctly
  indices_kept <- which(!sapply(all_counts_list, is.null)) # Apply this *before* filtering the list
  gsm_ids_from_files <- gsm_ids_from_files[indices_kept]
  # Now filter the list
  all_counts_list <- all_counts_list[indices_kept]
}


# Check if the number of rows (genes) is consistent across all imported files
first_file_gene_count <- length(gene_ids)
row_counts_consistent <- all(sapply(all_counts_list, length) == first_file_gene_count)

if (!row_counts_consistent) {
  stop("Files do not have a consistent number of genes (rows). Cannot safely combine.")
}

# Combine the list of count vectors into a single matrix/data frame
# do.call(cbind, ...) efficiently binds columns from a list
gene_counts_matrix <- do.call(cbind, all_counts_list)

# Assign gene IDs as row names
rownames(gene_counts_matrix) <- gene_ids

# --- Step 2 (Actual Renaming): Map GSM IDs to Sample Names and Assign Column Names ---

# Find the matching SampleName for each GSM ID based on the order files were read
# Use the 'match' function for safe lookup
matched_indices <- match(gsm_ids_from_files, mapping_data$GSM_ID)

# Check for any GSM IDs that weren't found in the mapping table
if (any(is.na(matched_indices))) {
  missing_gsms <- gsm_ids_from_files[is.na(matched_indices)]
  warning("The following GSM IDs from file names were not found in the mapping table:\n",
          paste(missing_gsms, collapse = "\n"),
          "\nColumns corresponding to these files will have NA names.")
}

# Get the corresponding SampleNames in the correct order
new_colnames <- mapping_data$SampleName[matched_indices]

# Assign the new sample names as column names to the matrix
colnames(gene_counts_matrix) <- new_colnames

# --- Step 3 (Final Formatting): Convert to Data Frame with gene_id column ---

# Convert the matrix to a data frame
gene_counts_df <- as.data.frame(gene_counts_matrix)

# Add the gene IDs (currently row names) as the first column named "gene_id"
gene_counts_df <- cbind(gene_id = rownames(gene_counts_df), gene_counts_df)

# Reset row names to default sequential integers
rownames(gene_counts_df) <- NULL

# --- Step 4: Export Gene Counts Matrix ---

# Create the output directory if it doesn't exist
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Write the data frame to a CSV file
# row.names = FALSE prevents writing the default integer row names
# quote = TRUE (default for columns with spaces) matches the example output format
message("Exporting combined matrix to: ", output_file)
write.csv(gene_counts_df, file = output_file, row.names = FALSE, quote = TRUE)

message("Script finished successfully.")

# Display the first few rows and columns of the resulting data frame as a check
print("Preview of the combined gene count matrix:")
print(head(gene_counts_df[, 1:min(ncol(gene_counts_df), 6)])) # Show gene_id + up to 5 samples