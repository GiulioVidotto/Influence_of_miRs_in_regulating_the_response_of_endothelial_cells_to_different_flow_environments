# Function: get_classifier_data
#
# This function processes data for the classification of miRNA-target interactions. It extracts relevant features,
# transforms the data into a long format, and adds binary classification variables based on target groups. The output is 
# a data frame ready for training a machine learning model for classification.
#
# Input:
#   - data_list: A list containing various data components required for the function. The list must include the following:
#     - "contrast": A string representing the contrast between two conditions (e.g., "OSS_vs_LSS").
#     - "altered_miRs_type": A table with miRNA alterations data.
#     - "merged_data": A data frame containing gene information with miRNA-target interactions.
#     - "background_group_1": A data frame with background information for group 1.
#     - "background_group_2": A data frame with background information for group 2.
#     - "test_group_1": A data frame with target genes for test group 1.
#     - "test_group_2": A data frame with target genes for test group 2.
#     - "test_group_3": A data frame with target genes for test group 3.
#   - altered_miR_table: A table with altered miRNA information.
#   - expression_database: A database containing expression data for genes.
#   - Kb_database: A database containing affinity scores and interactions between miRNAs and genes.
#
# Output:
#   - A data frame in long format containing the following columns:
#     - miRNA_ID: ID of the miRNA.
#     - Target_ID: ID of the target gene.
#     - total_num_of_miRNA: Total number of miRNAs targeting the gene.
#     - y, y_1, y_2, y_3: Binary variables indicating whether a gene belongs to test groups 1, 2, or 3 (0 or 1).
#     - Normalised miRNA counts, Normalised mRNA counts, Affinity scores: Features extracted and added for classification.

get_classifier_data <- function(data_list,
                                altered_miR_table,
                                expression_database,
                                Kb_database) {
  if(is.list(data_list)) {
    if(all(names(data_list) %in% c("contrast",
                                   "altered_miRs_type",
                                   "merged_data",
                                   "background_group_1",
                                   "background_group_2",
                                   "test_group_1",
                                   "test_group_2",
                                   "test_group_3" ))) {
      # 1) LOAD THE ENTIRE DATABASE AND EXTRACT THE INFORMATION ABOUT THE CONTRAST AND THE CONDITIONS
      # 1.1) Get the element with all the genes (without considering the different divisions)
      all_gene_data <- data_list$merged_data %>% 
        dplyr::filter(targeted_by_miRNA != 0) %>% # get only the row where there is a gene targeted by a miRNA
        separate_rows(all_miRNA_ID, sep = ", ") %>% 
        dplyr::rename(miRNA_ID = all_miRNA_ID) %>% 
        dplyr::distinct(Target_ID, miRNA_ID, total_num_of_miRNA)
      # 1.2) Contrast
      contrast <- data_list$contrast
      # 1.3) Altered miRs type 
      altered_miRs_type <- data_list$altered_miRs_type
      # 1.4) From the contrast, extract the two conditions
      condition_1 <- gsub("_vs_.*","", contrast)
      condition_2 <- gsub(".*_vs_","", contrast)
      # Extract the right condition
      if (any(condition_1 == "OSS", condition_1 == "ESS") & condition_2 == "LSS") {
        # 2) ADD INFORMATION (MISSING FEATURES)
        ## The classifier will use 4 different features for the classification:
        ### - Normalised miRNA counts
        ### - Normalised mRNA counts
        ### - Affinity scores
        ### - Number of miRNAs targeting a single mRNA
        ### The feature called "Number of miRNAs targeting a single mRNA" is already present in the data so there is no need to add it.
        ### For this reason only the others need to be added:
        # 2.1) Normalised miRNA counts
        first_feature <- get_first_feature(altered_miR_table,
                                           altered_miRs_type,
                                           condition_1,
                                           condition_2)
        # 2.2) Normalised mRNA counts
        second_feature <- get_second_feature(expression_database,
                                             condition_1)
        # 2.3) Affinity scores
        thrid_feature <- get_third_feature(Kb_database,
                                           altered_miRs_type)
        # 3) ADD THE INFORMATION TO THE TABLE
        # 3.1) Add the first feature
        all_gene_data <- all_gene_data %>% 
          left_join(first_feature, by = "miRNA_ID")
        # 3.2) Add the second feature
        all_gene_data <- all_gene_data %>% 
          left_join(second_feature, by = "Target_ID")
        # 3.3) Add the third feature
        all_gene_data <- all_gene_data %>% 
          left_join(thrid_feature, by = c("Target_ID" = "Target_ID",
                                          "miRNA_ID" = "miRNA_ID")) %>% 
          dplyr::filter(!is.na(Kd)) # remove rows where there is an interaction but we do not have information about the affinity scores
        # 4) TRANSFORM THE DATA INTO LONG FORMAT
        long_format_data <- all_gene_data %>%
          pivot_longer(
            cols = starts_with("miRNA_patient") | starts_with("mRNA_patient"),
            names_to = c(".value", "patient", "condition"),
            names_pattern = "(miRNA|mRNA)_(patient\\.\\d+)_(.*)"
          ) %>%
          dplyr::select(-c(patient,condition))
        # 5) ADD THE BINARY VARIABLE TO CLASSIFY
        # The binary variable depends on the test group selected (remember that this script is considering only genes that are targeted by miRNAs)
        # Only those genes are selected because in long_format_data are present only those genes that are targeted by miRNAs (see point 1.1)
        long_format_data <- long_format_data %>% 
          # Add the binary variable without the test group specification
          dplyr::mutate(y = ifelse(Target_ID %in% unique(data_list[["test_group_1"]]$Target_ID)  |
                                     Target_ID %in% unique(data_list[["test_group_2"]]$Target_ID)  |
                                     Target_ID %in% unique(data_list[["test_group_3"]]$Target_ID),
                                   1,
                                   0),
                        # Binary variable only for test group 1
                        y_1 = ifelse(Target_ID %in% unique(data_list[["test_group_1"]]$Target_ID), 1, 0),
                        # Binary variable only for test group 2
                        y_2 = ifelse(Target_ID %in% unique(data_list[["test_group_2"]]$Target_ID), 1, 0),
                        # Binary variable only for test group 3
                        y_3 = ifelse(Target_ID %in% unique(data_list[["test_group_3"]]$Target_ID), 1, 0))
        # 6) RETURN THE DATA
        return(long_format_data)
      } else {
        msg <- "Error: the first condition must be the disease condition (ex. \"OSS\ or \"ESS\") and the second one the baseline condition (ex. \"LSS\")"
        stop(msg)
      }
    } else {
      msg <- "Error: The object \"data_list\" must contain the following elements: \"merged_data\", \"altered_miRs_type\", \"background_group_1\", \"background_group_2\", \"test_group_1\", \"test_group_2\", \"test_group_3\""
      stop(msg)
    }
  } else {
    msg <- "Error: the object \"data_list\" must be a list"
    stop(msg)
  }
}


# Function: get_first_feature
#
# This function extracts miRNA concentration data for a given disease condition from the provided altered miRNA table. 
# The function filters the miRNAs based on their regulation status (either UP or DOWN) and their fold change and adjusted p-value 
# from a differential expression analysis. The function returns a table with miRNA IDs and their associated concentrations 
# in the disease condition for three patients.
#
# Input:
#   - altered_miR_table: A table with miRNA data including the fold change and p-value from differential expression analysis.
#   - altered_miRs_type: A string that indicates whether the miRNAs are "UP" or "DOWN" regulated.
#   - condition_1: The disease condition (e.g., "OSS").
#   - condition_2: The baseline condition (e.g., "LSS").
#
# Output:
#   - A table containing miRNA IDs and the corresponding concentrations for each of the three patients in the disease condition.

get_first_feature <- function(altered_miR_table, altered_miRs_type, condition_1, condition_2) {
  # Obtain the concentration of miRNAs for each patient in the disease condition
  if (altered_miRs_type == "UP") {
    miRNA_counts_db <- altered_miR_table %>% 
      # To make the size of the resulting table smaller, the funcion is taking into account only UP regulated miRNA
      dplyr::filter(!!rlang::sym(paste0(condition_1, "_vs_", condition_2, "_paired_log2FoldChange")) > 1 &
                      !!rlang::sym(paste0(condition_1, "_vs_", condition_2, "_paired_padj")) < 0.05) %>% 
      # Take only the values in the disease condition
      dplyr::select(miRNA_ID,
                    !!rlang::sym(paste0("miRNA_patient.1_", condition_1)) := !!rlang::sym(paste0("patient_1_", condition_1)),
                    !!rlang::sym(paste0("miRNA_patient.2_", condition_1)) := !!rlang::sym(paste0("patient_2_", condition_1)),
                    !!rlang::sym(paste0("miRNA_patient.3_", condition_1)) := !!rlang::sym(paste0("patient_3_", condition_1))) %>%
      # Make sure to have only distinct miRNAs names (in rare cases there are two miRNAs with the same name but encoded on differente genes. We decided to consider them as one)
      dplyr::distinct(miRNA_ID, .keep_all = TRUE)
  } else if (altered_miRs_type == "DOWN") {
    miRNA_counts_db <- altered_miR_table %>% 
      # To make the size of the resulting table smaller, the funcion is taking into account only DOWN regulated miRNA
      dplyr::filter(!!rlang::sym(paste0(condition_1, "_vs_", condition_2, "_paired_log2FoldChange")) < -1 &
                      !!rlang::sym(paste0(condition_1, "_vs_", condition_2, "_paired_padj")) < 0.05) %>% 
      # Take only the values in the disease condition
      dplyr::select(miRNA_ID,
                    !!rlang::sym(paste0("miRNA_patient.1_", condition_1)) := !!rlang::sym(paste0("patient_1_", condition_1)),
                    !!rlang::sym(paste0("miRNA_patient.2_", condition_1)) := !!rlang::sym(paste0("patient_2_", condition_1)),
                    !!rlang::sym(paste0("miRNA_patient.3_", condition_1)) := !!rlang::sym(paste0("patient_3_", condition_1))) %>% 
      # Make sure to have only distinct miRNAs names (in rare cases there are two miRNAs with the same name but encoded on differente genes. We decided to consider them as one)
      dplyr::distinct(miRNA_ID, .keep_all = TRUE)
  }
  return(miRNA_counts_db)
}


# Function: get_second_feature
#
# This function extracts mRNA concentration data for a given disease condition from the provided expression database. 
# The function returns a table with mRNA IDs and their associated concentrations in the disease condition for three patients.
#
# Input:
#   - expression_database: A table with mRNA data, which includes expression levels for different conditions and patients.
#   - condition_1: The disease condition (e.g., "OSS").
#
# Output:
#   - A table containing mRNA IDs and the corresponding concentrations for each of the three patients in the disease condition.

get_second_feature <- function(expression_database, condition_1) {
  # Obtain the concentration of mRNAs for each patient in the disease condition
  mRNA_counts_db <- expression_database %>% 
    dplyr::select(Target_ID,
                  !!rlang::sym(paste0("mRNA_patient.1_", condition_1)) := !!rlang::sym(paste0("patient_1_", condition_1)),
                  !!rlang::sym(paste0("mRNA_patient.2_", condition_1)) := !!rlang::sym(paste0("patient_2_", condition_1)),
                  !!rlang::sym(paste0("mRNA_patient.3_", condition_1)) := !!rlang::sym(paste0("patient_3_", condition_1)))
  return(mRNA_counts_db)
  
}


# Function: get_third_feature
#
# This function retrieves the affinity score information for miRNA-target gene interactions from the provided Kb_database. 
# The function filters and renames columns to retain only the relevant data about gene-miRNA interactions and the corresponding affinity score.
#
# Input:
#   - Kb_database: A list of data frames, each containing information about miRNA-target gene interactions, including affinity scores.
#   - altered_miRs_type: A string that specifies whether the miRNAs being considered are "UP" or "DOWN" regulated.
#
# Output:
#   - A data frame containing the filtered miRNA-gene interaction data with the corresponding affinity scores (excluding the MFE value).

get_third_feature <- function(Kb_database, altered_miRs_type) {
  # Select the MSE for the interactions between genes and the right altered miRNAs
  list_name <- names(Kb_database)
  match_num <- grep(tolower(altered_miRs_type), list_name)
  Kb_table <- Kb_database[[match_num]]
  # Extract only the columns regarding the interaction info (gene and miRNA IDs) and the one about the affinity score
  Kb_table <- Kb_table %>% 
    dplyr::select(-MFE_value) %>% 
    dplyr::rename(Target_ID = gene_ID)
  return(Kb_table)
}
