# Script Description:
# 
# This script automatically sources all R function scripts located in the specified 
# directory ('R_function_scripts'). It checks if the directory exists and contains 
# any R files. If R files are found, they are sourced one by one to make the 
# functions available in the current R environment.
# 
# Steps:
# 1. Defines the path to the folder containing the R function scripts.
# 2. Checks if the folder exists.
# 3. Lists all the R files in the folder (with full file paths).
# 4. Sources each R file if they are found.
# 5. Prints a message if no R files are found or if the directory does not exist.

# Define the path to the folder containing R scripts
R_folder <- paste0(getwd(), "/R_function_scripts")
# Check if the directory exists
if (dir.exists(R_folder)) {
  # List all R files in the folder, including full file paths
  list_of_files <- list.files(R_folder, pattern = "\\.R$", full.names = TRUE)
  # Check if there are any R files to source
  if (length(list_of_files) > 0) {
    # Loop through each file and source it
    for (file in list_of_files) {
      source(file)
    }
  } else {
    cat("No R files found in the directory.\n")
  }
} else {
  cat("Directory does not exist: ", R_folder, "\n")
}