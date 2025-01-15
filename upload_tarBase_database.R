# TarBase database v9.0 (released 2023)
#
# Overview:
# - name of the original file <- "Homo_sapiens_TarBase-v9.tsv" (it contains all the information about the
#   miRNA-gene interactions present in tarBase)
# - experimentally evaluated data
# - Description:
#   The content of TarBase-v9.0 is exclusively experimentally supported and consists of interactions identified 
#   via high- or low-yield methods. Raw datasets produced by state-of-the-art interactome mapping methods (e.g.,
#   HITS-CLIP, PAR-CLIP, CLASH) were uniformly analyzed to maintain consistently high-quality standards, while
#   miRNA targets verified via low-yield methods were amassed by means of manual article curation.
#
# The columns in the database are the following:
# 1) species
# 2) mirna_name
# 3) mirna_id
# 4) gene_name
# 5) gene_id
# 6) gene_location
# 7) transcript_name
# 8) transcript_id
# 9) chromosome
# 10) start
# 11) end
# 12) strand
# 13) experimental_method
# 14) regulation
# 15) tissue
# 16) cell_line
# 17) article_pubmed_id
# 18) confidence
# 19) interaction_group
# 20) cell_type
# 21) microt_score
# 22) comment
#
# The database has been filtered as we can see below:

# --- 1. Set Working Directory ---
# Set the working directory to the location of the GitHub cloned repository.
current_working_directory <- getwd()

# --- 2. Import the TarBase Database ---
# Define the path to the TarBase database file.
path_tarBase_database <- paste0(working_directory_path, "/miR_databases/tarBase_database/Homo_sapiens_TarBase-v9.tsv")
tarBase_database <- read.delim(path_tarBase_database,
                               stringsAsFactors = FALSE,
                               header = TRUE)

# --- 3. Filter Columns of Interest ---
# Take the columns of interest (the IDs of the miRNA and of the genes), only unique rows and rename the two column of interest.
tarBase_human_database <- tarBase_database %>% 
  dplyr::distinct(mirna_name, gene_id) %>% 
  dplyr::rename(miRNA_ID = mirna_name,
                Target_ID = gene_id)

# --- 4. Remove NA Values ---
# Remove possibile NA values
tarBase_human_database <- tarBase_human_database %>% 
  dplyr::filter(!is.na(Target_ID))

# --- 5. Compare with miRBase Database ---
## Comparison between miRDB_merge_database and miRBase_database
output_tarBase_database <- check_ID(miRBase_database, tarBase_human_database)

# --- 6. Save the Output File ---
# Save the file file in the final_outputs folder
write.csv(output_tarBase_database, file.path(paste0(working_directory_path,"/miR_databases/final_outputs/output_tarBase_database")), row.names = FALSE)
