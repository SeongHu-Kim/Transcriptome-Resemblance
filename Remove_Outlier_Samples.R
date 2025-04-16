# R script to import a CSV, remove specific columns, and save the result.

# Define the input file path
# Assumes the script is run from a directory where '../public/' is accessible
input_file_path <- "../public/GeneCountMatrix_Original.csv"

# Define the output file path
output_file_path <- "../public/GeneCountMatrix_filtered.csv"

# Define the columns to remove
columns_to_remove <- c("Control_R2", "Experimental_R1")

# 1. Import the CSV file
# Using tryCatch for basic error handling during file reading
tryCatch({
  # Read the CSV file into a data frame
  # header=TRUE assumes the first row contains column names
  # check.names=FALSE prevents R from modifying column names (e.g., replacing spaces with dots)
  gene_data <- read.csv(input_file_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  
  # Print a message indicating successful import
  cat("Successfully imported:", input_file_path, "\n")
  cat("Original dimensions:", dim(gene_data)[1], "rows,", dim(gene_data)[2], "columns\n")
  cat("Original column names:", paste(colnames(gene_data), collapse = ", "), "\n")
  
  # 2. Remove the specified columns
  # Check which columns to remove actually exist in the dataframe
  cols_exist <- columns_to_remove %in% colnames(gene_data)
  if (!all(cols_exist)) {
    warning("The following columns specified for removal were not found in the data: ",
            paste(columns_to_remove[!cols_exist], collapse = ", "))
  }
  
  # Select columns that are NOT in the 'columns_to_remove' list
  gene_data_filtered <- gene_data[, !(colnames(gene_data) %in% columns_to_remove)]
  
  # Print a message about the columns removed
  cat("Removed columns:", paste(columns_to_remove[cols_exist], collapse = ", "), "\n")
  cat("New dimensions:", dim(gene_data_filtered)[1], "rows,", dim(gene_data_filtered)[2], "columns\n")
  cat("Remaining column names:", paste(colnames(gene_data_filtered), collapse = ", "), "\n")
  
  # 3. Save the new data frame to a CSV file
  # Using tryCatch for basic error handling during file writing
  tryCatch({
    # Write the filtered data frame to the output CSV file
    # row.names = FALSE prevents R from writing row numbers as a column
    # quote = FALSE prevents strings from being quoted (often preferred for bioinformatics data)
    write.csv(gene_data_filtered, file = output_file_path, row.names = FALSE, quote = FALSE)
    
    # Print a success message for saving
    cat("Successfully saved filtered data to:", output_file_path, "\n")
    
  }, error = function(e) {
    # Handle errors during file writing
    cat("Error saving file:", output_file_path, "\n")
    cat("Error message:", conditionMessage(e), "\n")
  })
  
}, error = function(e) {
  # Handle errors during file reading
  cat("Error reading file:", input_file_path, "\n")
  cat("Error message:", conditionMessage(e), "\n")
  cat("Please ensure the file exists at the specified path and is readable.\n")
})

# End of script
