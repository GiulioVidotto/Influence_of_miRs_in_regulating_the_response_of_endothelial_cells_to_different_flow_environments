#!/bin/bash

# link of the file to download
download_link=$1
# Folder where the files will be stored
download_dir=$2

# If not already present, create a new folder where the downloaded chromosome file will be stored (if it is not already present)
if [ ! -d "$download_dir" ]; then
    mkdir -p "$download_dir"
fi

# Get the base name of the file from the link
file_name=$(basename "$download_link")
# Destination path where the file will be stored
dest_file="$download_dir/$file_name"

# Check if the file already is present in the folder
if [ -f "$dest_file" ]; then
    echo "File $file_name already exists, skipping download."
  else
    # Download the file if the link is valid
    echo "Checking link: $download_link"
    if wget --spider "$download_link" 2>/dev/null; then
        # If the link is valid, proceed to download
        echo "Downloading: $download_link"
        wget -P "$download_dir" "$download_link"
    else
        # If the link is not valid, log an error or skip
        echo "Error: $download_link is not reachable or invalid."
    fi
fi

# Attempt to unzip the file if it exists
if [ -f "$dest_file" ]; then
    file_type=$(file -b "$dest_file")

    if echo "$file_type" | grep -q "gzip compressed data"; then
        echo "Unzipping (gzip): $dest_file"
        gunzip "$dest_file"
    elif echo "$file_type" | grep -q "Zip archive data"; then
        echo "Unzipping (zip): $dest_file"
        unzip "$dest_file" -d "$download_dir"
    else
        echo "File $file_name is not a recognized compressed format. Skipping extraction."
    fi
else
    echo "Error: File $dest_file not found after download."
    exit 1
fi

echo "Download and processing completed successfully!"
echo "Download check and process completed!"
