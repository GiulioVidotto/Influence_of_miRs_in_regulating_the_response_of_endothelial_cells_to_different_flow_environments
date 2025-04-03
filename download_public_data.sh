#!/bin/bash

# Function to check if the dataset is already downloaded or has to be downlowed
download_if_not_exists() {
  local url=$1
  local folder=$2
  local filename=$(basename "$url")

  if [ -f "$folder/$filename" ]; then
    echo "The file $filename already exists, skipping download."
  else
    echo "Downloading $filename...."
    if ! ./download_file.sh "$url" "$folder"; then
      echo "Error: Download failed for $url"
      exit 1
    fi
  fi
}

# Define project data directory
project_data="$(pwd)/project_data"

# __ miRBase file __
# Decide the folder where the data from miRBase (miRNA notiation) will be stored
miRBase_download_folder="$project_data/miR_database/miRBase_database"
# Download the file stored in miRBase
download_if_not_exists "https://www.mirbase.org/download/mature.fa" "$miRBase_download_folder"
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
download_if_not_exists "https://www.targetscan.org/vert_80/vert_80_data_download/Conserved_Site_Context_Scores.txt.zip" "$targetScan_download_folder"
download_if_not_exists "https://www.targetscan.org/vert_80/vert_80_data_download/Nonconserved_Site_Context_Scores.txt.zip" "$targetScan_download_folder"
# Download the file stored in miRDB
download_if_not_exists "https://mirdb.org/download/miRDB_v6.0_prediction_result.txt.gz" "$miRDB_download_folder"
# Download the file stored in miRTarBase
download_if_not_exists "https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2022/cache/download/9.0/miRTarBase_MTI.xlsx" "$miRTarBase_download_folder"
# Download the file stored in TarBase
download_if_not_exists "https://dianalab.e-ce.uth.gr/tarbasev9/data/Homo_sapiens_TarBase-v9.tsv.gz" "$tarBase_download_folder"


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
  download_if_not_exists "$url" "$chr_download_folder"
done


# __ Lab data (example datasets - first 10 rows) __
lab_data_folder="$project_data/lab_data"
# mRNA-seq example dataset
dest_file_mRNA_seq_data="$lab_data_folder/mRNA_seq_database"
download_if_not_exists "" "$dest_file_mRNA_seq_data"
# miRNA-seq example dataset
dest_file_miRNA_seq_data="$lab_data_folder/miRNA_database"
download_if_not_exists "" "$dest_file_miRNA_seq_data"
# proteomics example dataset
dest_file_proteomics_data="$lab_data_folder/proteomics_database"
download_if_not_exists "" "$dest_file_proteomics_data"


echo "All downloads completed!"
