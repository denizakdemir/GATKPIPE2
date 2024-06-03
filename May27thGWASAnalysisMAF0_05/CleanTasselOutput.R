# Required packages
#install.packages(c("data.table", "stringr"))
library(data.table)
library(stringr)

# Set the input and output directories
input_dir <- "6_tassel_analysis"
output_dir <- "cleaned_tassel_analysis"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Function to process files
process_file <- function(file_path, output_dir) {
  # Read the file
  data <- fread(file_path, header = TRUE)

  # Remove rows containing "combined_"
  rows_to_remove <- apply(data, 1, function(row) any(grepl("combined_", row)))
  data_cleaned <- data[!rows_to_remove, ]

  # Remove columns containing "combined_"
  cols_to_remove <- grepl("combined_", colnames(data_cleaned))
  data_cleaned <- data_cleaned[, !cols_to_remove, with = FALSE]

  # Write the cleaned data to the output directory
  output_file_path <- file.path(output_dir, basename(file_path))
  fwrite(data_cleaned, output_file_path, quote = FALSE, sep = "\t")
}

# List of files to process
files_to_process <- list.files(input_dir, full.names = TRUE)

# Process each file
for (file_path in files_to_process) {
  process_file(file_path, output_dir)
}

cat("Processing completed. Cleaned files are saved in:", output_dir, "\n")
