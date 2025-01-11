# Function Name: download_file
#
# This function downloads a file from a specified URL and saves it to a designated directory.
# If the directory doesn't exist, it creates the directory. If the file already exists,
# the function skips the download and prints a message indicating the file's existence.
#
# Inputs:
# - download_dir: The directory where the downloaded file will be saved. If the directory 
#   doesn't exist, it will be created.
# - file_name: The name of the file to be saved in the specified directory.
# - url: The URL from which the file will be downloaded.
#
# Output:
# - If the file doesn't exist in the specified directory, the function downloads the file.
# - If the file already exists, the function prints a message stating that the file already exists.
# - The function doesn't return any value.

download_file <- function(download_dir, file_name, url){
  # Create the directory if it doesn't exist
  if (!dir.exists(download_dir)) {
    dir.create(download_dir, recursive = TRUE)
  }
  # Define the destination of the file
  dest_file <- file.path(download_dir, file_name)
  # If the file doesn't exist download it
  if (!file.exists(dest_file)) {
    cat("Downloading", url, "\n")
    download.file(url, dest_file)
  } else {
    cat("File already exists:", dest_file, "\n")
  }
}
