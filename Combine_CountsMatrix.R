# Load necessary libraries
# Install dplyr and stringr if you haven't already
# install.packages("dplyr")
# install.packages("stringr")
library(dplyr)
library(stringr)

# --- 1. Import the files ---

# Define file paths
file1_path <- "../public/GeneCountMatrix_filtered.csv"
file2_path <- "../public/GeneCountMatrix_GSE97263.csv"
output_path <- "../public/GeneCountMatrix_Combined.csv"

# Check if files exist before attempting to read
if (!file.exists(file1_path)) {
  stop("Error: File not found at ", file1_path)
}
if (!file.exists(file2_path)) {
  stop("Error: File not found at ", file2_path)
}

# Read the CSV files into data frames
# Use stringsAsFactors = FALSE to prevent automatic conversion of character columns to factors
# Use check.names = FALSE for the second file to handle potentially problematic column names gracefully,
# although read.csv often mangles them anyway if they aren't syntactically valid.
# dplyr handles non-standard names better later on.
tryCatch({
  df1 <- read.csv(file1_path, stringsAsFactors = FALSE, check.names = TRUE)
  df2 <- read.csv(file2_path, stringsAsFactors = FALSE, check.names = TRUE)
  
  # A quick check of the first few rows and column names
  print("--- Head of GeneCountMatrix_filtered.csv ---")
  print(head(df1))
  print("Column names:")
  print(colnames(df1))
  
  print("--- Head of GeneCountMatrix_GSE97263.csv ---")
  print(head(df2))
  print("Column names:")
  print(colnames(df2))
  
}, error = function(e) {
  stop("Error reading CSV files: ", e$message)
})


# --- 2. Identify common gene_ids (ignoring version numbers) ---

# Create a new column in each data frame with the base gene ID (without version)
# The pattern "\\.\\d+$" matches a dot followed by one or more digits at the end of the string
df1 <- df1 %>%
  mutate(gene_id_base = str_replace(gene_id, "\\.\\d+$", ""))

df2 <- df2 %>%
  mutate(gene_id_base = str_replace(gene_id, "\\.\\d+$", "")) # Assuming first column is gene_id

# Check the base gene IDs created
print("--- Head of df1 with gene_id_base ---")
print(head(select(df1, gene_id, gene_id_base)))

print("--- Head of df2 with gene_id_base ---")
# Ensure the first column name is indeed 'gene_id' or adjust if needed
# Let's assume the first column name in df2 is 'gene_id' based on the example structure
if (!"gene_id" %in% colnames(df2)) {
  warning("First column of GSE97263 file not named 'gene_id'. Assuming it is the gene ID column.")
  colnames(df2)[1] <- "gene_id" # Rename first column if necessary
  df2 <- df2 %>%
    mutate(gene_id_base = str_replace(gene_id, "\\.\\d+$", ""))
}
print(head(select(df2, gene_id, gene_id_base)))


# --- 3. Combine the two csv files based on common base gene_ids ---

# Identify the count columns in each data frame
# Exclude the original gene_id, gene_symbol (only in df1), and the new base gene_id column
count_cols_df1 <- setdiff(names(df1), c("gene_id", "gene_symbol", "gene_id_base"))
count_cols_df2 <- setdiff(names(df2), c("gene_id", "gene_id_base")) # Exclude potential gene_symbol if it exists

# Perform an inner join to keep only genes present in both files
# The join is done on the 'gene_id_base' column
# Suffixes ".x" and ".y" are automatically added to columns with the same name
# that are not used in the 'by' argument (e.g., gene_id.x, gene_id.y)
combined_df <- inner_join(df1, df2, by = "gene_id_base")

# Check the intermediate combined dataframe structure
# print("--- Head of intermediate combined_df ---")
# print(head(combined_df))
# print("Column names of intermediate combined_df:")
# print(colnames(combined_df))

# Select and arrange the final columns:
# - Keep gene_id and gene_symbol from df1 (which become gene_id.x, gene_symbol.x)
# - Keep all count columns from df1
# - Keep all count columns from df2
# - Discard gene_id.y and gene_id_base
final_combined_df <- combined_df %>%
  select(
    gene_id = gene_id.x,          # Rename gene_id.x from df1 back to gene_id
    gene_symbol = gene_symbol,    # Keep gene_symbol (it comes from df1)
    all_of(count_cols_df1),       # Include all original count columns from df1
    all_of(count_cols_df2)        # Include all original count columns from df2
  )

# Check the first few rows and columns of the final dataframe
print("--- Head of final combined data frame ---")
print(head(final_combined_df))
print("--- Dimensions of final combined data frame ---")
print(dim(final_combined_df))
print("--- Column names of final combined data frame ---")
print(colnames(final_combined_df))


# --- 4. Output the combined gene counts matrix ---

# Create the output directory if it doesn't exist (optional, adjust as needed)
# output_dir <- dirname(output_path)
# if (!dir.exists(output_dir)) {
#   dir.create(output_dir, recursive = TRUE)
# }

# Write the final combined data frame to a new CSV file
# row.names = FALSE prevents writing row numbers as a separate column
# quote = TRUE ensures that column names with special characters (like '-')
# and potentially string values are quoted, matching the desired output header format.
tryCatch({
  write.csv(final_combined_df, output_path, row.names = FALSE, quote = TRUE)
  print(paste("Successfully wrote combined matrix to:", output_path))
}, error = function(e) {
  stop("Error writing output CSV file: ", e$message)
})