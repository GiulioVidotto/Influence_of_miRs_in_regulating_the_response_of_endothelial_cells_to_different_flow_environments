# miRTarBase v9.0 (released 2021)
# 
# Overview:
# - Original file name: "miRTarBase_MTI.xlsx" (contains all the information in miRTarBase)
# - Experimentally validated miRNA-target interactions (MTIs)
# - Description:
#   miRTarBase accumulates over 360,000 MTIs collected from literature using NLP techniques.
#   The interactions are experimentally validated by reporter assays, western blot, microarray, and sequencing.
#   It includes data from 30 species, but this script focuses on human interactions.
#
# Columns in the database:
# 1) miRTarBase ID
# 2) miRNA
# 3) Species (miRNA)
# 4) Target Gene
# 5) Target Gene (Entrez ID)
# 6) Species
# 7) Experiments
# 8) Support Type
# 9) References (PMID)
#
# The database has been filtered as we can see below:

# --- 1. Get Working Directory ---
# As stated in the README file on github, the working directory should have been set to the cloned github repository
current_working_directory <- getwd()

# --- 2. Import the Database ---
# Define the path to the miRTarBase database
path_miRTarBase_database <- paste0(current_working_directory, "/miR_databases/mirTarBase_database/miRTarBase_MTI.xlsx")
# Import the miRTarBase database (this will read the Excel file and store it as a dataframe)
miRTarBase_database <- readxl::read_excel(path_miRTarBase_database,
                                          na = "")
# --- 3. Filter Data for Human miRNA-Target Interactions ---
# Filter the data to retain only interactions between human miRNAs and target genes
miRTarBase_human_database <- miRTarBase_database %>% 
  dplyr::filter(`Species (miRNA)` == "Homo sapiens" & `Species (Target Gene)` == "Homo sapiens" & startsWith(miRNA,"hsa")) %>% 
  dplyr::distinct(miRNA, `Target Gene (Entrez ID)`) %>% 
  dplyr::rename(miRNA_ID = miRNA,
                NCBI_IDs = `Target Gene (Entrez ID)`)

# --- 4. Extract NCBI IDs ---
# Extract the column containing NCBI IDs
miRTarBase_database_NCBI_gene <- miRTarBase_human_database %>% 
  dplyr::select(NCBI_IDs)
# Write the columns of NCBI gene IDs in a cvs file and export it to use it in Ensembl Biomart (the exported file is called "miRTarBase_database_NCBI_gene")
write.csv(miRTarBase_database_NCBI_gene, file.path(paste0(working_directory_path, "/miR_databases/mirTarBase_database/miRTarBase_database_NCBI_gene")), row.names = FALSE)

# Between Step 4 and 5 we used Ensembl Biomart (https://www.ensembl.org/info/data/biomart/index.html) to convert NCBI gene IDs into Ensembl gene IDs. We used the just saved file
# "miRTarBase_database_NCBI_gene", which contains the list of NCBI gene IDs, as input in Ensembl Biomart. Among the choosen attributes, there were Gene.stable.ID (Ensembl gene ID),
# NCBI.gene..formerly.Entrezgene..ID (NCBI gene IDs) and Chromosome.scaffold.name to consider only the canonical chromosomes and no alternative versions.

# --- 5. Import Ensembl Gene IDs ---
# After having used Ensembl Biomart and converted the NCBI IDs into Ensembl IDs, import the file with the new information (called "mirTarBase_database_ENSEMBL_gene.txt")
path_miRTarBase_ENSEMBL_gene <- paste0(working_directory_path, "/miR_databases/mirTarBase_database/mirTarBase_database_ENSEMBL_gene.txt")
miRTarBase_database_ENSEMBL_gene <- read.table(path_miRTarBase_ENSEMBL_gene, stringsAsFactors = FALSE, sep = ",", header = TRUE)

# --- 6. Filter Ensembl Database for Standard Chromosomes ---
# Ensure only genes corresponding to standard chromosomes are kept
miRTarBase_database_ENSEMBL_gene <- miRTarBase_database_ENSEMBL_gene %>%
  dplyr::filter(!startsWith(Chromosome.scaffold.name, "H")) %>% 
  dplyr::distinct(Gene.stable.ID, NCBI.gene..formerly.Entrezgene..ID) %>% 
  dplyr::rename(NCBI_IDs = NCBI.gene..formerly.Entrezgene..ID,
                Target_ID = Gene.stable.ID)

# --- 7. Check for Null Values ---
# Check if there are any missing values in both databases
sum(is.na(miRTarBase_human_database)) # 0, there are no null values in "miRTarBase_human_database"
sum(is.na(miRTarBase_database_ENSEMBL_gene)) # 0, there are no null values in "miRTarBase_database_ENSEMBL_gene"

# --- 8. Dimensions of Databases ---
# Check the dimensions of the databases
dim(miRTarBase_human_database) # There are 381084 interactions (miRNA-gene)
dim(miRTarBase_database_ENSEMBL_gene) # There are 14981 Ensembl/NCBI IDs

# --- 9. Add Ensembl IDs to Human Database ---
# Merge the two datasets on NCBI IDs to add Ensembl gene IDs to the human database
miRTarBase_human_database <- left_join(miRTarBase_human_database, miRTarBase_database_ENSEMBL_gene, relationship = "many-to-many") %>% 
  # The option "many-to-many" relationships is used because sometimes a single ESEMBL ID corresponds to multiple NCBI IDs and vice versa.
  dplyr::filter(!is.na(Target_ID)) # Not always there is a correspondence between NCBI IDs and ENSEMBL IDs

# --- 10. Compare miRTarBase Database with miRBase Database ---
# Compare miRTarBase_human_database with miRBase database (this function needs to be defined elsewhere)
output_miRTarBase_database <- check_ID(miRBase_database, miRTarBase_human_database)

# --- 11. Save Final Output ---
# Save the final comparison output to the "final_outputs" folder
write.csv(output_miRTarBase_database, file.path(paste0(working_directory_path, "/miR_databases/final_outputs/output_miRTarBase_database.csv")), row.names = FALSE)
