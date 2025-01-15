# TargetScan database (released 2021)
#
# - Name of the original files <- "Conserved_Site_Context_Scores.txt" and "Nonconserved_Site_Context_Scores.txt".
#   They both contains information about miRNA-gene interactions in converved and nonconverved sites.
# - Description:
#   TargetScan predicts biological targets of miRNAs by searching for the presence of conserved 8mer, 7mer,
#   and 6mer sites that match the seed region of each miRNA (Lewis et al., 2005). As options, predictions
#   with only poorly conserved sites and predictions with nonconserved miRNAs are also provided. 
# - The columns in the original database are the following:
#  1) Gene.ID
#  2) Gene.Symbol
#  3) Transcript.ID
#  4) Gene.Tax.ID
#  5) miRNA
#  6) Site.Type
#  7) UTR_start
#  8) UTR.end
#  9) context...score
#  10) context...score.percentile
#  11) weighted.context...score
#  12) weighted.context...score.percentile
#  13) Predicted.relative.KD

# --- 1. Set Working Directory ---
# As stated in the README file on GitHub, the working directory should be set 
# to the location of the cloned GitHub repositor
current_working_directory <- getwd()

# --- 2. Import Conserved and Nonconserved Context Scores ---
# Define paths to the TargetScan database files
path_targetScan_conserved_site <- paste0(current_working_directory, "/project_data/miR_databases/targetScan_database/Conserved_Site_Context_Scores.txt")
# Check if the file exists before importing
if (!file.exists(path_targetScan_conserved_site)) {
  stop("Error: The file was not found in the specified directory.")
}
path_targetScan_non_conserved_site <- paste0(current_working_directory, "/project_data/miR_databases/targetScan_database/Nonconserved_Site_Context_Scores.txt")
# Check if the file exists before importing
if (!file.exists(path_targetScan_non_conserved_site)) {
  stop("Error: The file was not found in the specified directory.")
}

# Import the conserved site context scores
targetScan_conserved_site <- read.delim(path_targetScan_conserved_site, stringsAsFactors = FALSE)

# Import the nonconserved site context scores
targetScan_non_conserved_site <- read.delim(path_targetScan_non_conserved_site, stringsAsFactors = FALSE)

# --- 3. Filter Human Interactions ---
# Keep only interactions involving Homo sapiens genes and miRNAs (only distinc interactions)
targetScan_non_conserved_site_human <- targetScan_non_conserved_site %>% 
  dplyr::filter(startsWith(miRNA, "hsa-")) %>% 
  dplyr::distinct(Gene.ID, miRNA, .keep_all = TRUE) 

# Before combining the two files, this operation needs to be done also on the other database
targetScan_conserved_site_human <- targetScan_conserved_site %>% 
  dplyr::filter(startsWith(miRNA, "hsa-")) %>% 
  dplyr::distinct(Gene.ID, miRNA, .keep_all = TRUE)

# --- 4. Combine Conserved and Nonconserved Data ---
# Combine the two files into a single one called "targetScan_database"
targetScan_database <- rbind(targetScan_conserved_site_human, targetScan_non_conserved_site_human)

# Remove the version from the ENSEMBL gene IDs. Then take the new distinct interactions and rename the two columns of interest.
targetScan_database <- targetScan_database %>%
  dplyr::mutate(Gene.ID = sub("\\.\\d+","",Gene.ID)) %>% 
  dplyr::distinct(miRNA, Gene.ID) %>%
  dplyr::rename(miRNA_ID = miRNA,
                Target_ID = Gene.ID)

# --- 5. Count miRNA-Gene Interactions ---
# Calculate how many times each gene is targeted by miRNAs
how_many_times_list <- lapply(unique(targetScan_non_conserved_site_2$Gene.ID), function(gene_ID) {
  how_many_times <- targetScan_non_conserved_site_2 %>% 
    dplyr::filter(Gene.ID == gene_ID) %>% 
    group_by(miRNA) %>% 
    summarise(times_that_it_binds = n())
  return(how_many_times)
})
how_many_times_df <- do.call(rbind, how_many_times_list)

# --- 6. Compare TargetScan Database with miRBase Database ---
# Compare the TargetScan database against miRBase database for ID consistency
targetScan_database_checked <- check_ID(miRBase_database, targetScan_database)

# --- 7. Separate Perfect Matches and Missing Information ---
# Extract perfect matches and rename the columns appropriately
# 1) The perfect matches
targetScan_database_perfect_match <- targetScan_database_checked %>% 
  dplyr::filter(!(endsWith(checked_miRNA_ID, ":'-5p' missing in miRBase") |
                    endsWith(checked_miRNA_ID, ":'-3p' missing in miRBase") |
                    endsWith(checked_miRNA_ID, ":'-3p' or '-5p' specification needed") |
                    endsWith(checked_miRNA_ID, ":'-5p' specification needed") |
                    endsWith(checked_miRNA_ID, ":'-3p' specification needed") |
                    endsWith(checked_miRNA_ID, ": Not found"))) %>% 
  dplyr::distinct(checked_miRNA_ID, Target_ID) %>% 
  dplyr::rename(miRNA_ID = checked_miRNA_ID)

# 2) matches not found because of the specification (3p or 5p)
targetScan_database_missing <- targetScan_database_checked%>% 
  dplyr::filter(endsWith(checked_miRNA_ID, ":'-3p' or '-5p' specification needed") |
                  endsWith(checked_miRNA_ID, ":'-5p' specification needed") |
                  endsWith(checked_miRNA_ID, ":'-3p' specification needed")
  )

# --- 8. Retrieve Additional Information for Missing Matches ---
# Other than the interaction database there is also another database (called "miR_Family_Info.txt") in targetScan regarding
# more information about the miRNAs (including the sequences and their miRBase Accession numbers). This database has been
# downloaded because, after the comparison with miRBase database, some miRNA IDs (7570) present in targetScan_database are
# not specified as they should. To extract additional information to be able to confirm the "specifications" of these problematic
# miRNA, the column about the miRBase Accession Numbers has been used.
# All the following operation are done to add information to the targetScan_database_missing database
path_targetscan_MAN_database <- paste0(current_working_directory, "/project_data/miR_databases/targetScan_database/miR_Family_Info.txt")
# Check if the file exists before importing
if (!file.exists(path_targetscan_MAN_database)) {
  stop("Error: The file was not found in the specified directory.")
}
# Import the dataset
targetscan_MAN_database <- read.delim(path_targetscan_MAN_database, stringsAsFactors = FALSE, header = TRUE)

# Get only the miRNA in Homo sapiens and the two coulmns of interest (MiRBase.Accession and MiRBase.ID)
targetscan_MAN_database <- targetscan_MAN_database %>% # MAN = Mirbase Accession Number
  dplyr::filter(startsWith(MiRBase.ID, "hsa-")) %>% 
  dplyr::distinct(MiRBase.Accession, MiRBase.ID) %>% 
  dplyr::rename(miRNA_ID = MiRBase.ID,
                miRBase_Accession = MiRBase.Accession)

#Join the two different databases (targetScan_database_missing and targetscan_MAN_database):
# Left join between targetScan_dataset_checked and sequence_targetscan
targetscan_miRNA_MAN_database <- left_join(targetScan_database_missing, targetscan_MAN_database)
# Left join between targetscan_miRNA_MAN_database and miRBase_Name_ID_database
targetscan_miRBase_database <- left_join(targetscan_miRNA_MAN_database, miRBase_Name_ID_database)
# miRNAs that we were able to retrive with the missing information
targetscan_miRBase_database %>% 
  dplyr::distinct(checked_miRNA_ID, miRBase_miRNA_ID)
# Output of the whole process (remove the rows where the Accession number is not available, it means that those are not present in miRBase)
targetscan_final_database <- targetscan_miRBase_database %>% 
  dplyr::filter(!is.na(miRBase_miRNA_ID)) %>% 
  dplyr::distinct(miRBase_miRNA_ID, Target_ID) %>% 
  dplyr::rename(miRNA_ID = miRBase_miRNA_ID)

# --- 9. Combine Perfect Matches with Resolved Matches ---
# Concatenate the dataset of perfect matches and resolved matches
targetscan_output_database <- rbind(targetScan_database_perfect_match, targetscan_final_database)

# Recheck the final database for consistency with miRBase
output_targetScan_database <- check_ID(miRBase_database, targetscan_output_database)

# --- 10. Save the Final Output ---
# Specify the path for the new folder
folder_path <- paste0(current_working_directory, "/project_data/miR_databases/final_outputs")

# Check if the folder exists, and create it if it doesn't
if (!file.exists(folder_path)) {
  dir.create(folder_path, recursive = TRUE)
  cat("Folder created:", folder_path, "\n")
} else {
  cat("Folder already exists:", folder_path, "\n")
}
# Save the final output file in the specified folder
write.csv(output_targetScan_database, file.path(paste0(current_working_directory, "/project_data/miR_databases/final_outputs/output_targetScan_database.csv")), row.names = FALSE)
