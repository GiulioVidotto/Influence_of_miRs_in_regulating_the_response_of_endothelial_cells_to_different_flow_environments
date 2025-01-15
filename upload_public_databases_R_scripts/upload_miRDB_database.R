# miRDB v6.0 (release 2019)
# 
# Overview:
# - name of the file -> "miRDB_human.txt" (it contains all the miRNA IDs and transcript IDs present in miRBB. Each interaction is associated with a target score)
# - Data obtained with an algorithm and then experimentally validated
# - DESCRIPTION:
# miRDB is an online database for predicted microRNA targets in animals. MicroRNAs are involved in many diverse biological processes and they may potentially regulate the functions of thousands of genes. One major issue in miRNA studies is the lack of bioinformatics programs to accurately predict miRNA targets. Animal miRNAs have limited sequence complementarity to their gene targets, which makes it challenging to select relevant biological features to build target prediction models with high specificity. We have developed a new miRNA target prediction program based on support vector machines (SVMs) and high-throughput training datasets. By systematically analyzing high-throughput experimental data, we have identified novel features that are important to miRNA target binding and expression downregulation. These new features as well as other known features have been integrated in an SVM machine learning framework for the training of our target prediction model. Our prediction algorithm has been validated by independent experimental data for its improved selectivity on predicting a large number of miRNA downregulated gene targets.
# 
# The columns in the database are the following:
# 1) miRNA_ID (NCBI accessions number for miRNAs)
# 2) target_mRNA_ID (NCBI accessions number for the transcript of the target gene)
# 3) prediction_score (interaction prediction scores) - This column will not be considered

# --- 1. Get Working Directory ---
# Set the working directory to the location of the miRDB data
current_working_directory <- getwd()

# --- 2. Import the Database ---
# Define the path to the miRDB prediction data file
path_miRDB_miRNA_mRNA_database <- paste0(current_working_directory, "project_data/miR_databases/miRDB_database/miRDB_human.txt")

# Check if the file exists before importing
if (!file.exists(path_miRDB_miRNA_mRNA_database)) {
  stop("Error: The file was not found in the specified directory.")
}

# Import the miRDB prediction data (tab-delimited file)
miRDB_miRNA_mRNA_database <- read.table(path_miRDB_miRNA_mRNA_database,
                                        sep = "\t",
                                        na.strings = "")

# --- 3. Operations on the database ---
# Adding names to the columns
colnames(miRDB_miRNA_mRNA_database) <- c("miRNA_ID", "target_mRNA_ID", "prediction_score")

# Check for Null values
sum(is.na(miRDB_miRNA_mRNA_database)) # 0, there are no Null values in miRDB_miRNA_mRNA_database

# Select only the columns of interest filtering out Null values
miRDB_miRNA_mRNA_database <- miRDB_miRNA_mRNA_database %>% 
  dplyr::distinct(miRNA_ID, target_mRNA_ID)

# Dimension of the database
dim(miRDB_miRNA_mRNA_database) # There are 3375741 interactions (miRNA-mRNA)

# --- 4. Extract Target Genes ---
# Extract the information about the NCBI IDs of the transcripts of the target genes
miRDB_database_RefSeq_IDs_transcripts <- miRDB_miRNA_mRNA_database %>% 
  dplyr::select(target_mRNA_ID)

# Write a new csv file containing the information about the IDs of the transcipts
write.csv(miRDB_database_RefSeq_IDs_transcripts, file.path(paste0(current_working_directory,"/project_data/miR_databases/miRDB_database/miRDB_database_RefSeq_IDs_transcripts.csv")), row.names = FALSE)

# The following database (miRDB_gene_mRNA_database) has been created starting from the original miRDB database and using Ensembl Biomart
# (the second column, the NCBI accession number for the transcipts/mRNAs of the target genes, was used to get the name of the target genes).

# Between Step 4 and 5 we used Ensembl Biomart (https://www.ensembl.org/info/data/biomart/index.html) to convert NCBI gene IDs into Ensembl gene IDs. We used the just saved file
# "miRDB_database_RefSeq_IDs_transcripts.csv". This was done to convert RefSeq_transcript IDs into Ensembl gene IDs.
# Save this file inside the "miRDB_database" folder as "mirDB_database_ENSEMBL_gene.csv"

# --- 5. Import the database after modifying it with Biomart --- 
# Define the path to the file
path_miRDB_ENSEMBL_gene <- paste0(current_working_directory,"/project_data/miR_databases/miRDB_database/mirDB_database_ENSEMBL_gene.csv")

# Check if the file exists before importing
if (!file.exists(path_miRDB_ENSEMBL_gene)) {
  stop("Error: The file was not found in the specified directory.")
}
# Import the dataset with gene names and mRNA IDs
miRDB_database_ENSEMBL_gene <- read.table(path_miRDB_ENSEMBL_gene, stringsAsFactors = FALSE, sep = ",", header = TRUE)

# --- 6. Filter the new database --- 
miRDB_database_ENSEMBL_gene <- miRDB_database_ENSEMBL_gene %>%
  dplyr::filter(!startsWith(chromosome_name, "H")) %>% 
  dplyr::distinct(ensembl_gene_id, refseq_mrna) %>% 
  dplyr::rename(Target_ID = ensembl_gene_id,
                target_mRNA_ID = refseq_mrna)

# --- 7. Check for NaN values ---
sum(is.na(miRDB_database_ENSEMBL_gene)) # 69168 Null values, this large number is due to Ensembl Biomart

# --- 8. Operation on the dataset ---
# Select only the columns regarding the gene names and the transcript IDs. Then remove all the NA values 
# (There was no correspondence between the gene name and the transcript ID in Ensembl mart)

# Function to separate strings and take only the first element of the string
# This function will be used to delete the version of the transcript ID of the target genes
split_only_fist_element <- function(column_name, separator) {
  result <- strsplit(column_name, separator)
  only_first_part <- lapply(result, function(x) x[1])
  return(unlist(only_first_part))
}
  
# Getting only the transcript ID without the version and filtering out possibile NA values  
miRDB_database_ENSEMBL_gene <- miRDB_database_ENSEMBL_gene %>% 
  dplyr::filter(!is.na(target_mRNA_ID) & !is.na(Target_ID)) %>% 
  dplyr::distinct(target_mRNA_ID, Target_ID) %>% 
  dplyr::mutate(target_mRNA_ID = split_only_fist_element(target_mRNA_ID, "\\."))

# dimension of the database
dim(miRDB_database_ENSEMBL_gene)[1] # There are 39641 genes 

# --- 9. Merge the databases ---
# Merging the two datasets (the one with the gene name and mRNA IDs and the one with miR IDs and mRNA IDs)
total_number_genes <- nrow(miRDB_miRNA_mRNA_database %>%  dplyr::distinct(target_mRNA_ID))
number_of_missing_genes <- nrow(miRDB_miRNA_mRNA_database %>% 
  left_join(miRDB_database_ENSEMBL_gene, relationship = "many-to-many") %>% 
  dplyr::distinct(target_mRNA_ID, Target_ID) %>% 
  dplyr::filter(is.na(Target_ID))
)

miRDB_merge_database <- miRDB_miRNA_mRNA_database %>% 
  left_join(miRDB_database_ENSEMBL_gene, relationship = "many-to-many") %>% 
  dplyr::distinct(miRNA_ID, Target_ID) %>% 
  dplyr::filter(!is.na(miRNA_ID) & !is.na(Target_ID))

## Comparison between miRDB_merge_database and miRBase_database
output_miRDB_database <- check_ID(miRBase_database, miRDB_merge_database)

# --- 10. Save the Output File ---
# Specify the path for the new folder
folder_path <- paste0(current_working_directory, "/project_data/miR_databases/final_outputs")

# Check if the folder exists, and create it if it doesn't
if (!file.exists(folder_path)) {
  dir.create(folder_path, recursive = TRUE)
  cat("Folder created:", folder_path, "\n")
} else {
  cat("Folder already exists:", folder_path, "\n")
}
# Save the file file in the final_outputs folder
write.csv(output_miRDB_database, file.path(paste0(current_working_directory, "/project_data/miR_databases/final_outputs/output_miRDB_database.csv")), row.names = FALSE)
