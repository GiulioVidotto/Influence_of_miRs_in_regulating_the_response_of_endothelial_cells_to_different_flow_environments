# Proteomics Dataset Overview
#
# - Name of the file -> "Proteomics_data.txt" (contains protein expression levels across different conditions)
#
# - DESCRIPTION:
# This proteomics dataset contains quantitative information on protein expression levels.

# --- 1. Get Working Directory ---
# Set the working directory
current_working_directory <- getwd()

# --- 2. Import the Database ---
# Define the path to the proteomics data file
original_proteomics_data_path <- paste0(current_working_directory, "project_data/lab_data/proteomics_database/20191219_HumphriesJ_prot_01.xlsx")

# Check if the file exists before importing
if (!file.exists(original_proteomics_data_path)) {
  stop("Error: The file was not found in the specified directory.")
}

# Import the proteomics data in R
original_proteomics_data <- read_xlsx(original_proteomics_data_path)

# --- 3. Analysis on the missing ENSEMBL IDs stored in the dataset ---
# The ENSEMBLE IDs stored in the proteomics dataset are going to be used to merge the information in the proteomics data with the one in the transcriptomics data.
# Important: Some entires are missing the ENSEMBL IDs but there is the Entrez IDs for most of them. For this reason, we decided to obtain the corresponding ENSEMBL IDs
#            starting from the NCBI IDs.

# Extract the NCBI IDs for the genes that do not have a corresponding ENSEMBL ID
NCBI_IDs_with_missing_ENSEMBL_IDs <- original_proteomics_data$`Entrez Gene ID`[is.na(original_proteomics_data$`Ensembl Gene ID`)]

# Export the table just created ("NCBI_IDs_with_missing_ENSEMBL_IDs") and use it in ENSEMBL BioMart (link: https://www.ensembl.org/info/data/biomart/index.html)
# Use ENSEMBL BioMart to obtain the correspondent ENSEMBL IDs from the provided NCBI IDs.
# Save the file obtained from ENSEMBL BioMart as "NCBI_ID_to_ENSEBL_ID.txt" in the "proteomics_database" folder.

# Load the data with the ENSEMBL IDs and the respective NCBI IDs (table obtained from ENSEMBL BioMart). While loading the dataset, change the names of the columns
# and make sure that the NCBI IDs are considered as characters. 
missing_ENSEMBL_ID_db_path <- paste0(current_working_directory, "project_data/lab_data/proteomics_database/NCBI_ID_to_ENSEBL_ID.txt")
missing_ENSEMBL_ID_db <- read.table(missing_ENSEMBL_ID_db_path, sep = "\t", header = TRUE) %>% 
  dplyr::rename(Target_ID = Gene.stable.ID,
                NCBI_ID = NCBI.gene..formerly.Entrezgene..ID) %>% 
  dplyr::mutate(NCBI_ID = as.character(NCBI_ID))

# --- 3. Merge all the datasets together ---
# 3) Merge the dataset called "missing_ENSEMBL_ID_db" with the dataset called "original_proteomics_data". 
# We noticed that using ENSEMBL BioMArt allowed us to recover almost all the missing ENSEMBL IDs.
original_proteomics_data_with_missing_ENSEMBL <- original_proteomics_data %>% 
  dplyr::filter(is.na(original_proteomics_data$`Ensembl Gene ID`)) %>% 
  dplyr::rename(NCBI_ID = `Entrez Gene ID`) %>% 
  dplyr::select(-`Ensembl Gene ID`) %>% 
  left_join(missing_ENSEMBL_ID_db, by = "NCBI_ID", relationship = "many-to-many") %>% # It can happen that a single NCBI ID maps to multiple ENSEMBL IDs 
  dplyr::distinct(Target_ID, .keep_all = TRUE) %>% 
  dplyr::relocate(Target_ID, .before = "Gene Symbol")

# --- 4. Remove the rows still without ENSEMBL IDs ---
# The rows still without ENSEMBL IDs are also the ones that do not have other information about the genes
# (like the gene name and also the chromosome etc..). We have decided to exclude these rows from further analysis and
# rename some of the columns.
original_proteomics_data <- original_proteomics_data %>% 
  dplyr::filter(!is.na(`Ensembl Gene ID`)) %>% 
  dplyr::rename(NCBI_ID = `Entrez Gene ID`,
                Target_ID = `Ensembl Gene ID`) %>% 
  rbind(original_proteomics_data_with_missing_ENSEMBL)

# --- 5. Select only the columns of interest for this analysis ---
# In this step we selected only the columns about the protein IDS, ENSEMBL IDs, abundance ratio and adjusted p-values. 
statistics_proteomics_data <- original_proteomics_data %>% 
  dplyr::select(Protein_ID = Accession,
                Target_ID,
                protein_OSS_vs_LSS_Abundance_Ratio = "Abundance Ratio: (Oscillatory) / (Laminar)",
                protein_OSS_vs_LSS_adj_pvalue =  "Abundance Ratio Adj. P-Value: (Oscillatory) / (Laminar)",
                protein_ESS_vs_LSS_Abundance_Ratio = "Abundance Ratio: (Elevated) / (Laminar)",
                protein_ESS_vs_LSS_adj_pvalue = "Abundance Ratio Adj. P-Value: (Elevated) / (Laminar)") %>% 
  separate_rows(Target_ID, sep ="; ") %>%  # Some of the rows have multiple Ensembl IDs divided by ";"
  dplyr::mutate(Target_ID = gsub("\\..*", "", Target_ID),
                Protein_ID = gsub("-\\d*", "", Protein_ID)) %>%  # Remove the version of the ENSEMBL IDS
  dplyr::filter(!rowSums(
    is.na(
      dplyr::select(.,
                    protein_OSS_vs_LSS_Abundance_Ratio,
                    protein_OSS_vs_LSS_adj_pvalue,
                    protein_ESS_vs_LSS_Abundance_Ratio,
                    protein_ESS_vs_LSS_adj_pvalue)
    )
  ) == 4) %>%  # Filtered out all the rows where all 4 columns related to fold changes and p-values in all 4 columns are NA.
  dplyr::distinct(Protein_ID,
                  Target_ID,
                  protein_OSS_vs_LSS_Abundance_Ratio,
                  protein_OSS_vs_LSS_adj_pvalue,
                  protein_ESS_vs_LSS_Abundance_Ratio,
                  protein_ESS_vs_LSS_adj_pvalue)