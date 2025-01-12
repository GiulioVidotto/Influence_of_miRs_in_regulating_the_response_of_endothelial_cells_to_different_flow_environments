#!/bin/bash

# Parameters
# Each line of this file represent an interaction (gene_ID    miRNA_ID)
TARGET_QUERY_FILE=$1

# Each line of this file contains the information about a gene and its sequence
GENE_SEQ_FILE=$2

# This is a fasta file with all the information regarding miRNAs
MIRNA_SEQ_FILE=$3

if [[ "${TARGET_QUERY_FILE}" == "up"* ]]; then

  # If the script is analyzing interaction between up regulated miRNAs and their targets
  OUTPUT_FILE="up_thermodynamics_scores_file"
  
  OUTPUT_FILE_RNA_up="up_thermodynamics_scores_file_rna_up"
  
  # Name of the folder where the accessibility profiles will be stored (for each interaction)
  PROFILES_DIR="up_profiles"
  
else
  
  # If the script is analyzing interaction between down regulated miRNAs and their targets
  OUTPUT_FILE="down_thermodynamics_scores_file"
  
  OUTPUT_FILE_RNA_up="down_thermodynamics_scores_file_rna_up"
  
  # Name of the folder where the accessibility profiles will be stored (for each interaction)
  PROFILES_DIR="down_profiles"
  
fi


# Create the directory for the accessibility profiles (-p: if it is already been created do not throw an error)
mkdir -p $PROFILES_DIR

# EXTRACT_SEQUENCE
# Function to extract sequences of miRNAs based on their names/IDs
# The input file is the one with the miRNA or mRNAs sequences. Being a fasta file, the format is the following:
# >miRNA_name or mRNA_name and other_information
# sequence

function extract_sequence {
    # Define the local variables
    local name=$1
    local seq_file=$2
    local output_file=$3
    grep -A 1 -m 1 "^>${name}" $seq_file | grep -v "^--$" > ${output_file}
}


function extract_sequence {
    # Define the local variables
    local name=$1
    local seq_file=$2
    local output_file=$3
    grep -A 1 -m 1 "^>${name}" $seq_file | grep -v "^--$" | awk 'BEGIN { FS="::" } { if ($1 ~ /^>/) { print $1 } else { print } }' > ${output_file}
}


# Define how many lines are present in the file "$TARGET_QUERY_FILE"
total_lines=$(wc -l < "$TARGET_QUERY_FILE")
current_line=1

# Loop through each line of the target-query file
while IFS=$'\t' read -r line; do

    # Message to the user
    echo "Processing interaction ${current_line} of ${total_lines}"
    
    # Set the start time
    start_time=$(date +%s)

    # Extract the gene ID and miRNA ID from the line that it's being analyzed

    # The gene ID is stored in the first field of the line
    gene=$(echo ${line} | awk '{print $2}')

    # The miRNA ID is stored in the second field of the line
    mirna=$(echo ${line} | awk '{print $1}')

    # Set filenames based on the gene and miRNA names
    gene_file="${gene}.fa"
    mirna_file="${mirna}.fa"

    # Extract the sequences for both the gene and the miRNA and store them in the respective fasta files.
    extract_sequence ${gene} ${GENE_SEQ_FILE} ${gene_file}
    extract_sequence ${mirna} ${MIRNA_SEQ_FILE} ${mirna_file}
    
    
    # Check if one of the two files is empty
    if [ ! -s $gene_file ] || [ ! -s $mirna_file ]; then
    
        echo "Sequence for $gene or $mirna not found, skipping."
        
        rm -f ${gene_file} ${mirna_file}
        
    else
    
        # GENE SEQUENCE ANALYSIS
        # Check the length of the gene sequences
        gene_sequence=$(grep -v "^>" ${gene_file})
        gene_length=$(echo ${#gene_sequence})
        
        # MIRNA SEQUENCE ANALYSIS
        # Check the length of the miRNA sequences
        miRNA_sequence=$(grep -v "^>" ${mirna_file})
        miRNA_length=$(echo ${#miRNA_sequence})

        if [ ${gene_length} -ge 2000 ]; then
            
            # Run RNAplfold for the gene sequence
            RNAplfold -W 240 -L 160 -b  -O < ${gene_file}
  
            if [ ${miRNA_length} -lt 200 ]; then
                # Run RNAplfold for the miRNA sequences
                RNAplfold -W ${miRNA_length} -u 30 -b -O < ${mirna_file}
            else
                # Run RNAplfold for the miRNA sequences
                RNAplfold -W 240 -L 160 -u 30 -b -O < ${mirna_file}
            fi
  
            # Move the RNAplfold output files to the profiles directory
            mv ${gene}_* ${PROFILES_DIR}/
            mv ${mirna}_* ${PROFILES_DIR}/
    
            # Run RNAplex
            RNAplex -l 25 -c 30 -t ${gene_file} -q ${mirna_file} -a ${PROFILES_DIR} -b >> ${OUTPUT_FILE}
            
            # Clean the profiles directory
            rm -f ${PROFILES_DIR}/*
        
            # Clean up temporary files
            rm -f ${gene_file} ${mirna_file}

        else
        
            duplex="${miRNA_sequence}&${gene_sequence}"
        
            duplex_file="${gene}_${mirna}_duplex.txt"
            
            echo "${duplex}" > ${duplex_file}
            
            RNAup_output="${gene}_${mirna}_RNAup.out"
            RNAup -b < ${duplex_file} > ${RNAup_output}
            
            echo ">${mirna}" >> ${OUTPUT_FILE_RNA_up}
            echo ">${gene}" >> ${OUTPUT_FILE_RNA_up}
            cat ${RNAup_output} >> ${OUTPUT_FILE_RNA_up}
            
             # Clean up temporary files
            rm -f ${RNAup_output} ${duplex_file}
            rm -f RNA*.out
            rm -f ${gene_file} ${mirna_file}
  
        fi
      
    fi
    
    # Set the end time
    finish_time=$(date +%s)
    
    # Calculate the total time taken
    total_time=$(( finish_time - start_time ))
    
    # Calculate the total time to finish everything
    total_time_to_finish=$(( total_time * (total_lines - current_line) ))
    
    # Checks to understand the time unit
    
    if [[ ${total_time} -lt 60 ]];then
      # Time in seconds
      total_time="${total_time} sec"
      
    elif [[ ${total_time} -ge 60 ]] && [ ${total_time} -lt 3600 ]];then
      # Time in minutes
      total_time=$(( total_time / 60 ))
      total_time="${total_time} min"
      
    elif [[ ${total_time} -ge 3600 ]];then
      # Time in hours
      total_time=$(( total_time / 3600 ))
      total_time="${total_time} h"
      
    fi
    
    if [[ ${total_time_to_finish} -lt 60 ]];then
      # Time in seconds
      total_time_to_finish="${total_time_to_finish} sec"
      
    elif [[ ${total_time_to_finish} -ge 60 ]] && [[ ${total_time_to_finish} -lt 3600 ]];then
      # Time in minutes
      total_time_to_finish=$(( total_time_to_finish / 60 ))
      total_time_to_finish="${total_time_to_finish} min"
      
    elif [[ ${total_time_to_finish} -ge 3600 ]];then
      # Time in hours
      total_time_to_finish=$(( total_time_to_finish / 3600 ))
      total_time_to_finish="${total_time_to_finish} h"
      
    fi
    
    # Print a message for the user
    echo "Interaction ${current_line} of ${total_lines} correctly analysed! [time taken: ${total_time}, about ${total_time_to_finish} remaining]"
    
    # update the index of the current line
    current_line=$((current_line+1))

done < ${TARGET_QUERY_FILE}

# Tell the user that the loop is finished
echo "Finished!"