#!/bin/bash

# First part of the file name (use 3UTR_with_poly_A_sequences_on_chr)
file_name=$1

# Array with all the chromosome numbers
chr_num_array=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y")
# There is no M chromosome because predictAPA could not predict any 3UTR regions on that chromosome based on the poly A sites.

for chr_num in ${chr_num_array[@]};do
    # Define the name of the file
    file="${file_name}${chr_num}.fa"
    # Define the file location to the DNA sequence file   
    DNA_file="$(pwd)/project_data/genome_data/bed_and_fasta_file/utr_with_poly_A_output/3UTR_with_poly_A_sequences_DNA/${file}"
    echo "Processing ${file}"
    fasta_output="$(pwd)/project_data/genome_data/bed_and_fasta_file/utr_with_poly_A_output/3UTR_with_poly_A_sequences_mRNA/mRNA_${file}"
    # If it already exists, remove the file
    rm -f ${fasta_output} 
    # Iterate over each line of the file
    while IFS= read -r line; do
        if [[ "${line}" == ">"* ]];then
            # Append the header to the fasta file
            echo "${line}" >> ${fasta_output}
        else
            # Convert the sequence from DNA to RNA
            rna=$(echo "${line}" | tr tT uU )
            # Convert the sequence from RNA to mRNA
            mrna=$(echo "${rna}" | tr uacgUACG AUGCAUGC)
            # Obtain the reverse complement
            mrna_reverse_complement=$(echo "${mrna}" | rev)
            # Append the sequence to the new fasta file
            echo "${mrna_reverse_complement}" >> ${fasta_output}
        fi
    done < ${DNA_file}
done

echo "Finished!"
    
