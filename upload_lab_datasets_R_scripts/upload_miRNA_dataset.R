# miRNA Dataset Overview
#
# - Name of the file -> "miR_sequencing_dataset.xlsx" (contains miRNA expression levels across different conditions)
#
# - DESCRIPTION:
# This script loads and cleans a miRNA sequencing dataset, renames and filters relevant columns, and extracts significantly altered
# miRNAs (log2FC > 1, FDR < 0.05) with at least 15 total counts for the OSS vs LSS and ESS vs LSS comparisons.

# --- 1. Get Working Directory ---
# Set the working directory
current_working_directory <- getwd()

# --- 2. Import the Database ---
# Define the path to the miRNAs data file
miR_sequencing_database_path <- paste0(current_working_directory, "project_data/lab_data/miRNA_seq_example_dataset.xlsx")

# Check if the file exists before importing
if (!file.exists(miR_sequencing_database_path)) {
  stop("Error: The file was not found in the specified directory.")
}

# Import the proteomics data in R
miR_sequencing_database <- readxl::read_excel(miR_sequencing_database_path, na = "")

# --- 3. Change the column names ---
miR_sequencing_database <- miR_sequencing_database %>%
  dplyr::select(-matches("unpaired|UNpaired"),
                -c(pvalue...35,
                   padj...36,
                   pvalue...41,
                   padj...42,
                   pvalue...47,
                   padj...48)
  ) %>% 
  dplyr::rename(miRNA_ID = Name,
                patient_1_OSS_raw = SW1_S1...10,
                patient_1_LSS_raw = SW2_S2...11,
                patient_1_ESS_raw = SW3_S3...12,
                patient_2_OSS_raw = SW4_S4...13,
                patient_2_LSS_raw = SW5_S5...14,
                patient_2_ESS_raw = SW6_S6...15,
                patient_3_OSS_raw = SW7_S7...16,
                patient_3_LSS_raw = SW8_S8...17,
                patient_3_ESS_raw = SW9_S9...18,
                patient_1_OSS = SW1_S1...19,
                patient_1_LSS = SW2_S2...20,
                patient_1_ESS = SW3_S3...21,
                patient_2_OSS = SW4_S4...22,
                patient_2_LSS = SW5_S5...23,
                patient_2_ESS = SW6_S6...24,
                patient_3_OSS = SW7_S7...25,
                patient_3_LSS = SW8_S8...26,
                patient_3_ESS = SW9_S9...27,
                ESS_vs_LSS_paired_log2FoldChange = ESS_vs_LSS_paired_log2FC,
                ESS_vs_LSS_paired_pvalue = pvalue...32,
                ESS_vs_LSS_paired_padj = padj...33,
                LSS_vs_OSS_paired_log2FoldChange = LSS_vs_OSS_paired_log2FC,
                OSS_vs_LSS_paired_pvalue = pvalue...38,
                OSS_vs_LSS_paired_padj = padj...39,
                OSS_vs_ESS_paired_pvalue = pvalue...44,
                OSS_vs_ESS_paired_padj = padj...45
  ) %>% 
  dplyr::mutate(
    patient_col_sum = rowSums(dplyr::select(., matches("patient"))),
    LSS_vs_OSS_paired_log2FoldChange = LSS_vs_OSS_paired_log2FoldChange * -1 # Change the foldchange values for this contrast (because it is the opposite)
  ) %>% 
  dplyr::rename(OSS_vs_LSS_paired_log2FoldChange = LSS_vs_OSS_paired_log2FoldChange)

# --- 4. Obtain the altered miRNAs in the OSS vs LSS contrast --- 
# We selected only the altered miRNAs with a minimum of 15 counts across all samples, considering the rest as not biologically significant.
OSS_vs_LSS_miR_sequencing_database <- miR_sequencing_database %>% 
  dplyr::mutate(ID = gsub("_\\d*", "", ID)) %>% 
  dplyr::distinct(ID, .keep_all = TRUE) %>% 
  dplyr::filter(abs(OSS_vs_LSS_paired_log2FoldChange) > 1 &
                OSS_vs_LSS_paired_padj < 0.05 &
                patient_col_sum > 15)

# --- 4. Obtain the altered miRNAs in the ESS vs LSS contrast --- 
# We selected only the altered miRNAs with a minimum of 15 counts across all samples, considering the rest as not biologically significant.
ESS_vs_LSS_miR_sequencing_database <- miR_sequencing_database %>% 
  dplyr::mutate(ID = gsub("_\\d*", "", ID)) %>% 
  dplyr::distinct(ID, .keep_all = TRUE) %>% 
  dplyr::filter(abs(ESS_vs_LSS_paired_log2FoldChange) > 1 &
                ESS_vs_LSS_paired_padj < 0.05 &
                patient_col_sum > 15)
