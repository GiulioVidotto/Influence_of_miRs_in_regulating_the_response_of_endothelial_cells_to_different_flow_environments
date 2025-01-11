#!/bin/bash

# link of the file to download
download_link=$1
# Folder where the files will be stored
download_dir="$HOME/Validation_step/bed_and_fasta_file/chr_fasta_file"

# If not already present, create a new folder where the downloaded chromosome file will be stored (if it is not already present)
if [ ! -d "$download_dir" ]; then
    mkdir -p "$download_dir"
fi

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

echo "Download check and process completed!"
