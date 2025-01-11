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

# If the file is zipped, then unzip it
# Unzip the file if it exists and is a valid gzip or zip file
if [ -f "$dest_file" ]; then
    # Check if the file is a gzip file using the 'file' command
    if file "$dest_file" | grep -q "gzip compressed data"; then
        echo "Unzipping (gzip): $dest_file"
        gunzip "$dest_file"
    elif file "$dest_file" | grep -q "Zip archive data"; then
        echo "Unzipping (zip): $dest_file"
        unzip "$dest_file" -d "$download_dir"
    else
        echo "Error: $dest_file is neither a gzip nor a zip file."
    fi
else
    echo "Error: $dest_file not found."
fi

echo "Download check and process completed!"
