#!/bin/bash
# Array with all the chromosome numbers
chr_num_array=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y") # Chromosome M is not present becaus no UTR regions were identify on that chromosome
# For loop that iterates through all the sample numbers
for chr_num in "${chr_num_array[@]}"; do
    # Define the bed_file name
    bed_file="$(pwd)/project_data/genome_data/bed_and_fasta_file/utr_with_poly_A_bed_file/3UTR_poly_A_regions_chr${chr_num}.bed"
    # Define the name of the fasta file of the chromosome
    chr_fasta_file="$(pwd)/project_data/genome_data/bed_and_fasta_file/chr_fasta_file/chr${chr_num}.fa"
    # Define the sequence_output name
    sequence_output="$(pwd)/project_data/genome_data/bed_and_fasta_file/utr_with_poly_A_output/3UTR_with_poly_A_sequences_DNA/3UTR_with_poly_A_sequences_on_chr${chr_num}.fa"
    # Check if the required files and directories exist
    if [ ! -f "$bed_file" ]; then
        echo "Error: BED file $bed_file not found. Skipping chromosome $chr_num."
        continue
    fi
    if [ ! -f "$chr_fasta_file" ]; then
        echo "Error: FASTA file $chr_fasta_file not found. Skipping chromosome $chr_num."
        continue
    fi
    if [ ! -d "$(dirname "$sequence_output")" ]; then
        echo "Creating output directory: $(dirname "$sequence_output")"
        mkdir -p "$(dirname "$sequence_output")"
    fi
    # Run bedtools getfasta
    echo "Processing chromosome $chr_num..."
    bedtools getfasta -name -s -fi ${chr_fasta_file} -bed ${bed_file} -fo ${sequence_output}
done
