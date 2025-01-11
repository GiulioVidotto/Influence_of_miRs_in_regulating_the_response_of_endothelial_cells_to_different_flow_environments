#!/bin/bash

# link of the file to download
download_link=$1
# Folder where the files will be stored
download_dir="$HOME/Validation_step/bed_and_fasta_file/chr_fasta_file"
# Link domain
allowed_domain="hgdownload.soe.ucsc.edu"  # UCSC download domain

# If not already present, create a new folder where the downloaded chromosome file will be stored (if it is not already present)
if [ ! -d "$download_dir" ]; then
    mkdir -p "$download_dir"
fi

# Extract the domain from the provided link
domain=$(echo "$download_link" | sed -E 's|^https?://([^/]+).*|\1|')

# Check if the domain matches the allowed domain
if [[ "$domain" == "$allowed_domain" ]]; then
    # Check if the link is valid using wget --spider
    echo "Checking link: $download_link"
    if wget --spider "$download_link" 2>/dev/null; then
        # If the link is valid, proceed to download
        echo "Downloading: $download_link"
        wget -P "$download_dir" "$download_link"
    else
        # If the link is not valid, log an error or skip
        echo "Error: $download_link is not reachable or invalid."
    fi
else
    # If the domain is not valid, skip the link or log an error
    echo "Skipping: $download_link (not from $allowed_domain)"
fi

echo "Download check and process completed!"
