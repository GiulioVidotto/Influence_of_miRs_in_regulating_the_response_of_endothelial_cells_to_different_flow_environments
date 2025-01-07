# Function Name: get_altered_miRNAs
# Purpose: Retrieves a vector or data frame of altered miRNAs based on a specific contrast and user-defined thresholds.
# Usage: This function is used in the "get_model_data_function.Rmd" and "script_for_statistical_tests.Rmd" scripts.
#
# Description:
# This function identifies altered miRNAs (UP-regulated or DOWN-regulated) from a provided miRNA dataset
# based on a specific contrast between two conditions. It applies user-defined thresholds for adjusted p-values 
# and log2 fold changes, and returns either a vector or a data frame of the results depending on the specified data type.
#
# Inputs:
# - altered_miR_database: Data frame containing miRNA data and associated statistics.
# - altered_miRNAs_type: Type of miRNAs to retrieve ("UP", "DOWN").
# - condition_1: First condition in the contrast (e.g., "LSS", "OSS", or "ESS").
# - condition_2: Second condition in the contrast (e.g., "LSS", "OSS", or "ESS").
# - data_type: Format of the output ("dataframe" or "vector").
# - miRNA_adj_p_value_threshold: Numeric threshold for adjusted p-values (significance level).
# - fold_change_threshold: Numeric threshold for log2 fold changes.
#
# Outputs:
# - A data frame or vector containing all altered miRNAs that meet the specified criteria in the given contrast.

get_altered_miRNAs <- function(altered_miR_database,
                               altered_miRNAs_type,
                               condition_1,
                               condition_2,
                               data_type,
                               miRNA_adj_p_value_threshold,
                               fold_change_threshold) {
  
  # Validate data_type input
  if (data_type %in% c("dataframe", "vector")) {
    # Validate altered_miRNAs_type input
    if (altered_miRNAs_type %in% c("DOWN", "UP")) {
      # Ensure conditions are valid and different
      if (condition_1 %in% c("LSS", "OSS", "ESS") & condition_2 %in% c("LSS", "OSS", "ESS")) {
        if (condition_1 != condition_2) {
          # Check that p-value threshold is numeric
          if (is.numeric(miRNA_adj_p_value_threshold)) {
            # Ensure fold change threshold is non-negative
            if (fold_change_threshold >= 0) {
              # Filter the database for distinct miRNAs and relevant columns
              altered_miR_database <- altered_miR_database %>%
                dplyr::distinct(miRNA_ID, .keep_all = TRUE) %>% 
                dplyr::select(
                  miRNA_ID,
                  matches(paste0("patient_[123]_", condition_1)),
                  matches(paste0("patient_[123]_", condition_2)),
                  miRNA_paired_padj = matches(paste0(condition_1, "_vs_", condition_2, "_paired_padj")),
                  miRNA_paired_log2FoldChange = matches(paste0(condition_1, "_vs_", condition_2, "_paired_log2FoldChange"))
                ) %>% 
                dplyr::filter(
                  miRNA_paired_padj < miRNA_adj_p_value_threshold
                )
              ## This part of the code works based on the choice between UP and DOWN regulated
              # Handle UP-regulated miRNAs
              if (altered_miRNAs_type == "UP") {
                ### Take all the miRNAs with a log2 fold change greater than a specific threshold on the fold changes
                if (data_type == "dataframe") {
                  altered_miR_database <- altered_miR_database %>% 
                    dplyr::filter(
                      miRNA_paired_log2FoldChange > fold_change_threshold
                    )
                  return(altered_miR_database)
                } else if (data_type == "vector") {
                  altered_miRNA_vector <-  unlist(
                    altered_miR_database %>% 
                      dplyr::filter(
                        miRNA_paired_log2FoldChange > fold_change_threshold
                      ) %>% 
                      dplyr::select(miRNA_ID)
                  )
                  return(altered_miRNA_vector)
                }
                # Handle DOWN-regulated miRNAs  
              } else if (altered_miRNAs_type == "DOWN") {
                ### Take all the miRNAs with a log2 fold change lower than a specific threshold on the fold changes
                if (data_type == "dataframe") {
                  altered_miR_database <- altered_miR_database %>% 
                    dplyr::filter(
                      miRNA_paired_log2FoldChange < -fold_change_threshold
                    )
                  return(altered_miR_database)
                } else if (data_type == "vector") {
                  altered_miRNA_vector <-  unlist(
                    altered_miR_database %>%
                      dplyr::filter(
                        miRNA_paired_log2FoldChange < -fold_change_threshold
                      ) %>% 
                      dplyr::select(miRNA_ID)
                  )
                  return(altered_miRNA_vector)
                }
              }
            } else {
              # Error message if the fold change threshold is not numeric, a single value, or is negative
              msg <- "The user should input the absolute value of the threshold on the fold changes"
              stop(msg)
            }
          } else {
            # Error message if the given p-value threshold is not numeric
            msg <- "Error: The threshold on the p-value of miRNAs must be a number"
            stop(msg)
          }
        } else {
          # Error message if the two conditions are the same (a contrast must be specified)
          msg <- "Error: the condition must be different (must be a contrast)"
          stop(msg)
        }
      } else {
        # Error message if the given conditions are not valid (must be "LSS", "OSS", or "ESS")
        if (!condition_1 %in% c("LSS", "OSS", "ESS")) {
          msg <- "Error: Condition 1 must be one of the following: \"LSS\", \"OSS\" or \"ESS\"."
          stop(msg)
        } else if (!condition_2 %in% c("LSS", "OSS", "ESS")) {
          msg <- "Error: Condition 2 must be one of the following: \"LSS\", \"OSS\" or \"ESS\"."
          stop(msg)
        }
      }
    } else {
      # Error message if the altered_miRNAs_type is invalid (must be "UP", "DOWN", or "ALL")
      msg <- "Error: altered_miRNAs_type not correct. Values allowed: UP or DOWN"
      stop(msg)
    }
  } else {
    # Error message if the data_type is invalid (must be "dataframe" or "vector")
    msg <- "Error: data_type is not correct. Only \"dataframe\" or \"vector\" are allowed."
    stop(msg)
  }
}



