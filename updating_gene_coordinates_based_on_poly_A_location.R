# SCRIPT TO GET THE COORDINATES OF EVERY GENE BASED ON THE FURTHEST POLY-A SITE

# This script uses as input the file obtained from APAtrap (this file is the concatenation of all the chromosome files obtain from predictAPA.
# We have multiple files because predictAPA takes a long time to run and so we tried to run the process in parallel) to update the stop 
# position of each mRNA with the furthest poly-A site and in the end the user will obtain a bed file with the updated information.

# 1) Use the function update_stop_codon to update the stop position of each mRNA
## First define the name of the file
APA_output_path <- paste0(path.expand("~"), "/project_data/genome_data/APA_output/all_APA_output_chr.txt")

## Apply the function update_stop_coord
tryCatch({
  mRNA_info_file <- update_stop_coord(file_path = APA_output_path) %>%
    dplyr::mutate(length_3UTR = as.numeric(stop) - as.numeric(start),
                  MSE = as.integer(MSE))
}, error = function(e) {
  stop("Error during updating stop coordinates: ", e$message)
})

# 2) Extract the list of RefSeq Transcript IDs and use it as input for miRBase (only the ones that are protein coding trascript, strat with NM)
Refseq_transcript_ID_vector <- unlist(mRNA_info_file %>% 
                                        dplyr::filter(startsWith(Transcript_ID, "NM")) %>% 
                                        dplyr::select(Transcript_ID))

# 3) Get list of chromosomes to check if there are chromosomes
chr_vector <- sort(unique(mRNA_info_file$chr))
if (length(chr_vector) == 0) {
  stop("Error: No chromosomes found in mRNA_info_file.")
}

# 4)Apply the transformation (the output is a list and each element is regarding a specific chromosome)
chr_coordinates_list <- lapply(chr_vector, function(chr_num) {transform_into_bed(chr_num = chr_num, genomic_data = mRNA_info_file)})

# For each element, write a new file that contains the data of the element so that these files can be used on the
# shell to extarct the corresponding sequences of each gene. There will be 25 files, one for each
# specific chromosome (write the new files into the "utr_with_poly_A_bed_file" directory).

output_dir <- paste0(path.expand("~"), "/project_data/genome_data/bed_and_fasta_file/utr_with_poly_A_bed_file")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)  # Create the directory if it doesn't exist
}

for (file in chr_coordinates_list) {
  # Save the file name in a variable (extract the number of the chromosome from the first coordinate)
  chr_num <- unique(sub("\t.*","",file$gene_info))
  file_name <- paste0("3UTR_poly_A_regions_", chr_num)
  # Remove the colnames
  colnames(file) <- NULL
  # Save the file in the "utr_with_poly_A_bed_file" folder
  write.table(file, file.path(paste0(output_dir, "/", file_name, ".bed")), row.names = FALSE, quote=FALSE, sep = "\t")
}

# This files will be used as input to the get_sequence.sh script. This script extracts the DNA sequences of the genes based on the coordinates given as input.
