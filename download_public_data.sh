#!/bin/bash

# Define project data directory
project_data="$(pwd)/project_data"

# __ miRBase file __
# Decide the folder where the data from miRBase (miRNA notiation) will be stored
miRBase_download_folder="$project_data/miR_database/miRBase_database"
# Download the file stored in miRBase
./download_file.sh "https://www.mirbase.org/download/mature.fa" "$miRBase_download_folder"
# Define input and output file paths
mature_file="$miRBase_download_folder/mature.fa"
filtered_output="$miRBase_download_folder/miRBase_human.txt"
# Filter the data to obtain only the IDs of human miRNAs
if [ -f "$mature_file" ]; then
grep -P "product=\"hsa-let.*|product=\"hsa-miR.*" "$mature_file" | grep -oP "hsa-\w{3}-(\w|\d)*" > "$filtered_output"
echo "Filtered human miRNA IDs saved to $filtered_output"
else
  echo "Error: File $mature_file not found. Download may have failed."
fi

# __ miRNA-mRNA interactions files __
# Decide the folder where the miRNA-mRNA interaction data will be stored
targetScan_download_folder="$project_data/miR_database/targetScan_database"
miRDB_download_folder="$project_data/miR_database/miRDB_database"
miRTarBase_download_folder="$project_data/miR_database/miRTarBase_database"
tarBase_download_folder="$project_data/miR_database/tarBase_database"
# Download the files stored in TargetScan
./download_file.sh "https://www.targetscan.org/vert_80/vert_80_data_download/Conserved_Site_Context_Scores.txt.zip" "$targetScan_download_folder"
./download_file.sh "https://www.targetscan.org/vert_80/vert_80_data_download/Nonconserved_Site_Context_Scores.txt.zip" "$targetScan_download_folder"
# Download the file stored in miRDB
./download_file.sh "https://mirdb.org/download/miRDB_v6.0_prediction_result.txt.gz" "$miRDB_download_folder"
# Download the file stored in miRTarBase
./download_file.sh "https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2022/cache/download/9.0/miRTarBase_MTI.xlsx" "$miRTarBase_download_folder"
# Download the file stored in TarBase
./download_file.sh "https://dianalab.e-ce.uth.gr/tarbasev9/data/Homo_sapiens_TarBase-v9.tsv.gz" "$tarBase_download_folder"

# __ Chromosome files __
# Base URL
base_url="https://hgdownload.soe.ucsc.edu/goldenPath/currentGenomes/Homo_sapiens/chromosomes"

# Target directory for downloads
chr_download_folder="$project_data/genome_data/bed_and_fasta_file/chr_fasta_file"

# Loop over chromosome numbers (1 to 22, X, Y, and M for mitochondrial DNA)
for chr in {1..22} X Y M; do
  # Construct the full URL
  url="${base_url}/chr${chr}.fa.gz"
  # Check and download the file if it doesn't already exist
  dest_file="$chr_download_folder/$chr${chr}.fa"
  if [ -f "$dest_file" ]; then
    echo "File chr${chr}.fa already exists, skipping download."
  else
    echo "Downloading chr${chr}...."
    ./download_file.sh "$url" "$chr_download_folder"
  fi
done

echo "All downloads completed!"
