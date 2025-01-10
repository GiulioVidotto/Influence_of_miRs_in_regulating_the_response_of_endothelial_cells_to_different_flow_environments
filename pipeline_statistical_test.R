# Function Name: pipeline_statistical_test (used in the "pipeline_statistical_test.Rmd")
#
# This function performs a series of statistical tests on altered miRNAs, mRNA, 
# and proteomics data. It validates inputs, extracts altered miRNAs, identifies 
# common targets, and generates statistical plots based on the selected test.
#
# Inputs:
# - altered_miR_database: The name of the miRNA sequencing database. Must be one of 
#   "OSS_vs_LSS_miR_sequencing_database", "OSS_vs_ESS_miR_sequencing_database", or 
#   "ESS_vs_LSS_miR_sequencing_database".
# - miR_mRNA_interaction: Dataframe or table containing interactions between 
#   miRNAs and mRNAs, including target gene information.
# - mRNA_proteomics_data: Dataframe or table with mRNA expression and proteomics 
#   data, including fold changes and adjusted p-values.
# - altered_miRNAs_type: Specifies the type of altered miRNAs. Allowed values are 
#   "UP" or "DOWN".
# - common_target_threshold: Minimum number of targets required for a miRNA to 
#   be considered in the analysis.
# - adj_p_value_threshold: The adjusted p-value threshold for selecting 
#   significant miRNAs.
# - miRNAs_threshold: Minimum number of miRNAs targeting a gene required for 
#   inclusion in the analysis.
# - background_group_num: Indicates which background group to use for the 
#   statistical comparison. Must be "1" or "2".
# - fold_change_threshold: Minimum fold change required for miRNAs or genes to 
#   be considered significant.
# - about_test: Specifies the type of test to perform. Allowed values are 
#   "KFERQ motif test" or "regulated genes test".
#
# Output:
# - A series of plots showing the results of the statistical tests, including 
#   mosaic plots comparing background and test groups.

pipeline_statistical_test <- function(altered_miR_database,
                                      miR_mRNA_interaction,
                                      mRNA_proteomics_data,
                                      altered_miRNAs_type,
                                      common_target_threshold,
                                      adj_p_value_threshold,
                                      miRNAs_threshold,
                                      background_group_num, 
                                      fold_change_threshold,
                                      about_test,
                                      save_intermediate) {
  # Check if the inputs are correct
  possibile_databases_names <- c("OSS_vs_LSS_miR_sequencing_database",
                                 "OSS_vs_ESS_miR_sequencing_database",
                                 "ESS_vs_LSS_miR_sequencing_database")
  # Change the altered_miRNAs_type from lower case to upper case 
  # (so the user can also write it in lower case and it would not make a difference)
  altered_miRNAs_type <- toupper(altered_miRNAs_type)
  if (deparse(substitute(altered_miR_database)) %in% possibile_databases_names) {
    if (altered_miRNAs_type %in% c("DOWN", "UP", "ALL")) {
      if (common_target_threshold > 0) {
        if (adj_p_value_threshold >= 0 & adj_p_value_threshold <= 1) {
          if (miRNAs_threshold >= 0) {
            if (as.character(background_group_num) %in% c("1","2")) {
              if (fold_change_threshold >= 0) {
                if (about_test %in% c("KFERQ motif test", "regulated genes test")) {
                  if (is.logical(save_intermediate)) {
                    # 1) ___FIRST PART - GET THE INFORMATION ABOUT THE MIRNAs
                    ## Extract the contrast from the name of the database of the miRNAs
                    contrast <- gsub("_miR_sequencing_database", "", deparse(substitute(altered_miR_database)))
                    # From the contast, extract the two conditions
                    condition_1 <- gsub("_vs_.*","", contrast)
                    condition_2 <- gsub(".*_vs_","", contrast)
                    # Extract miRNA information
                    common_target_db <- extract_miRNA_info(altered_miR_database = altered_miR_database,
                                                           altered_miRNAs_type = altered_miRNAs_type,
                                                           adj_p_value_threshold = adj_p_value_threshold,
                                                           fold_change_threshold = fold_change_threshold,
                                                           condition_1 = condition_1,
                                                           condition_2 = condition_2,
                                                           common_target_threshold = common_target_threshold)
                    # 2) ___ SECOND PART - EXECUTE THE STATISTICAL TEST AND GET A GRAPH AS OUTPUT ___ 
                    # Test and create graphs for different statistical tests (chi-square test or fisher test based on the used data)
                    # The analysed groups are:
                    # 1) Background group 1 - genes where there is no change in proteomics and also no change in mRNA expression
                    # 2) Background group 2 - genes where the change in proteomics scales with change in mRNA expression up or down
                    # 3) test group 1 - genes that don't change at mRNA, but do change in proteomics 
                    # 2) test group 2 - genes that don't change at the proteomics level, but do change in mRNA
                    # 3) test group 3 - genes where mRNA changes in one direction, but proteins change in opposite 
                    get_graphs(common_interaction_db = common_target_db,
                               altered_miRNAs_type = altered_miRNAs_type,
                               mRNA_proteomics_data = mRNA_proteomics_data,
                               contrast = contrast,
                               miRNAs_threshold = miRNAs_threshold,
                               background_group_num = background_group_num,
                               about_test = about_test,
                               save_intermediate = save_intermediate,
                               adj_p_value_threshold = adj_p_value_threshold
                    )
                  } else {
                    msg <- "Error: save_intermediate must be logical"
                    stop(msg)
                  }
                } else {
                  msg <- "Error: about_test incorrect. Only \"regulated genes test\" and \"KFERQ motif test\" can be used"
                  stop(msg)
                }
              } else {
                msg <- "The user should input the absolute value of the threshold on the fold changes"
                stop(msg)
              }
            } else {
              msg <- "Error: The background group number can be only 1 or 2"
              stop(msg)
            }
          } else {
            msg <- "Error: The threshold on the number of considered miRNAs must be greater than 0"
            stop(msg)
          }
        } else {
          msg <- "Error: The threshold on the p-value must be between 0 and 1"
          stop(msg)
        }
      } else {
        msg <- "Error: The threshold on the number of considered target genes for each miRNAs must be a number greater than 0"
        stop(msg)
      }
    } else {
      msg <- "Error: altered_miRNAs_type not correct. Values allowed: UP, DOWN or ALL"
      stop(msg)
    }
  } else {
    # If the name of the "altered_miR_database" is incorrect
    msg <- "Error: altered_miR_database not correct. Select among:\n
            OSS_vs_LSS_miR_sequencing_database,\n
            OSS_vs_ESS_miR_sequencing_database,\n
            ESS_vs_LSS_miR_sequencing_database"
    stop(msg)
  }
}


# Function Name: extract_miRNA_info
#
# This function extracts altered miRNAs based on user-defined conditions and 
# identifies their common target genes. It validates inputs, performs filtering, 
# and returns a database of interactions for the selected miRNAs.
#
# Inputs:
# - altered_miR_database: The name of the miRNA database to analyze.
# - altered_miRNAs_type: Specifies the type of altered miRNAs to analyze. 
#   Allowed values are "UP", "DOWN", or "ALL".
# - adj_p_value_threshold: Threshold for the adjusted p-value to filter significant miRNAs. 
#   Must be between 0 and 1.
# - fold_change_threshold: Minimum fold change required for miRNAs to be 
#   considered significant. Must be non-negative.
# - condition_1: The first condition for the contrast analysis. Must be one of 
#   "LSS", "OSS", or "ESS".
# - condition_2: The second condition for the contrast analysis. Must be one of 
#   "LSS", "OSS", or "ESS". Cannot be the same as condition_1.
# - common_target_threshold: Minimum number of target genes required for a 
#   miRNA to be included in the analysis. Must be greater than 0.
#
# Output:
# - A dataframe (common_target_db) containing interactions between the selected 
#   altered miRNAs and their common target genes. Each row provides details 
#   about a specific interaction.

extract_miRNA_info <- function(altered_miR_database,
                               altered_miRNAs_type,
                               adj_p_value_threshold,
                               fold_change_threshold,
                               condition_1,
                               condition_2,
                               common_target_threshold) {
  if (altered_miRNAs_type %in% c("DOWN", "UP", "ALL")) {
    if (adj_p_value_threshold >= 0 & adj_p_value_threshold <= 1) {
      if (fold_change_threshold >= 0) {
        if (common_target_threshold > 0) {
          if (condition_1 %in% c("LSS", "OSS", "ESS") & condition_2 %in% c("LSS", "OSS", "ESS")) {
            if (condition_1 != condition_2) {
              # 1) ___ FIRST PART - EXTRACT THE NAMES OF THE ALTERED MIRNAs ___ 
              ## Obtain the vector of alter miRNAs based on the contrast that the user wants to analyse
              altered_miRNA_vector <- get_altered_miRNAs(altered_miR_database = altered_miR_database,
                                                         altered_miRNAs_type = altered_miRNAs_type,
                                                         condition_1 = condition_1,
                                                         condition_2 = condition_2,
                                                         data_type = "vector",
                                                         miRNA_adj_p_value_threshold = adj_p_value_threshold,
                                                         fold_change_threshold = fold_change_threshold)
              # Remove the names form the miRNA vector
              names(altered_miRNA_vector) <- NULL
              # 2) ___ SECOND PART - FIND THE COMMON TARGETS BASED ON THE NAMES OF THE ALTERED MIRNAs ___ 
              ## Get the dataframe where in each row there is the information about a specific interaction
              ## The gene targets present in this database are the ones targeted by one or more miRNAs
              common_target_db <- find_common_target(miR_mRNA_interaction = miR_mRNA_interaction,
                                                     miRNA_vector = altered_miRNA_vector,
                                                     n_target = common_target_threshold)
              return(common_target_db)
            } else {
              msg <- "Error: the condition must be different (must be a contrast)"
              stop(msg)
            }
          } else {
            if (!condition_1 %in% c("LSS", "OSS", "ESS")) {
              msg <- "Error: Condition 1 must be one of the following: \"LSS\",  \"OSS\" or \"ESS\"."
              stop(msg)
            } else if (!condition_2 %in% c("LSS", "OSS", "ESS")) {
              msg <- "Error: Condition 2 must be one of the following: \"LSS\",  \"OSS\" or \"ESS\"."
              stop(msg)
            }
          }
        } else {
          msg <- "Error: The threshold on the number of considered target genes for each miRNAs must be a number greater than 0"
          stop(msg)
        }
      } else {
        msg <- "The user should input the absolute value of the threshold on the fold changes"
        stop(msg)
      }
    } else {
      msg <- "Error: The threshold on the p-value must be between 0 and 1"
      stop(msg)
    }
  } else {
    msg <- "Error: altered_miRNAs_type not correct. Values allowed: UP, DOWN or ALL"
    stop(msg)
  }
}


get_graphs <- function(common_interaction_db,
                       altered_miRNAs_type,
                       mRNA_proteomics_data,
                       contrast,
                       miRNAs_threshold,
                       background_group_num,
                       about_test,
                       save_intermediate,
                       adj_p_value_threshold) {
  
  # 1) Create the temporary databases
  first_temporary_database <- create_first_temporary_database(common_interaction_db, miRNAs_threshold)
  second_temporary_database <- create_second_temporary_database(mRNA_proteomics_data, contrast, adj_p_value_threshold)
  # 2) Merge the data
  merged_data <- merge_databases(first_temporary_database, second_temporary_database)
  # 3) Divide the data into groups
  groups <- divide_data_into_groups(merged_data, adj_p_value_threshold)
  # 4) Save intermediate data if requested
  if (save_intermediate) {
    intermediate_data <- list(contrast = contrast,
                              altered_miRNAs_type = altered_miRNAs_type,
                              merged_data = merged_data,
                              background_group_1 = groups$background_group_1,
                              background_group_2 = groups$background_group_2,
                              test_group_1 = groups$test_group_1,
                              test_group_2 = groups$test_group_2,
                              test_group_3 = groups$test_group_3)
    save_intermediate_data(intermediate_data)
  }
  # 5) Select the background group
  background_group <- ifelse(background_group_num == "1", groups$background_group_1, groups$background_group_2)
  # 6) Perform statistical tests and plotting
  test_groups_list <- list(test_group_1 = groups$test_group_1,
                           test_group_2 = groups$test_group_2,
                           test_group_3 = groups$test_group_3)
  results_list <- perform_statistical_tests(test_groups_list, background_group, altered_miRNAs_type, miRNAs_threshold, about_test)
  create_plots(results_list)
  # 7) Extract p-values
  p_values_vector <- extract_p_values(results_list)
  return(p_values_vector)
}


# Function Name: create_first_temporary_database
#
# This function creates the first temporary database by filtering out genes that are targeted by 
# fewer than 'miRNAs_threshold' number of miRNAs. It retains only distinct target genes for further analysis.
#
# Inputs:
# - common_interaction_db: A dataframe containing the interaction data of miRNAs and target genes.
# - miRNAs_threshold: The minimum number of miRNAs required to target a gene for inclusion.
#
# Output:
# - A dataframe (first_temporary_database) containing the target genes and associated miRNAs 
#   that meet the threshold criteria.

create_first_temporary_database <- function(common_interaction_db, miRNAs_threshold) {
  first_temporary_database <- common_interaction_db %>%
    dplyr::select(Target_ID = common_target_ENSEMBL_ID,
                  all_miRNA_ID,
                  total_num_of_miRNA) %>%
    dplyr::filter(total_num_of_miRNA >= miRNAs_threshold) %>%
    dplyr::distinct(Target_ID, .keep_all = TRUE)
  return(first_temporary_database)
}


# Function Name: create_second_temporary_database
#
# This function generates the second temporary database by processing the 'mRNA_proteomics_data'.
# It filters out rows with missing values or specific invalid protein abundance values 
# and adds a new column for the log2-transformed protein abundance ratios.
#
# Inputs:
# - mRNA_proteomics_data: A dataframe containing mRNA, protein abundance, and associated metadata.
#
# Output:
# - A dataframe (second_temporary_database) containing the filtered data with log2-transformed 
#   protein abundance ratios for each gene.

create_second_temporary_database <- function(mRNA_proteomics_data, contrast, adj_p_value_threshold) {
  second_temporary_database <- mRNA_proteomics_data %>%
    dplyr::select(Target_ID,
                  gene_name,
                  Protein_ID,
                  protein_Abundance_log2_Ratio = matches(paste0("protein_", contrast, "_Abundance_Ratio")),
                  mRNA_log2_FoldChange = matches(paste0(contrast, "_paired_log2FoldChange")),
                  protein_adj_pvalue = matches(paste0("protein_", contrast, "_adj_pvalue")),
                  mRNA_adj_pvalue = matches(paste0(contrast, "_paired_padj")),
                  KFERQ_motif) %>%
    dplyr::filter(!is.na(protein_Abundance_log2_Ratio) & !is.na(protein_adj_pvalue)) %>%
    dplyr::filter(protein_Abundance_log2_Ratio != 100 & protein_Abundance_log2_Ratio != 0.01) %>%
    dplyr::mutate(protein_Abundance_log2_Ratio = log2(protein_Abundance_log2_Ratio))
  return(second_temporary_database)
}


# Function Name: merge_databases
#
# This function merges two datasets: the first temporary database and the second temporary database.
# It performs a left join by the 'Target_ID' and adds new columns to indicate the number of miRNAs targeting 
# each gene, as well as whether the gene is targeted by any miRNA. Missing values are handled by replacing 
# them with default values.
#
# Inputs:
# - first_temporary_database: A dataframe containing filtered miRNA-target gene interactions.
# - second_temporary_database: A dataframe containing filtered mRNA data with log2-transformed 
#   protein abundance ratios.
#
# Output:
# - A merged dataframe (merged_data) with additional columns for miRNA targeting information and 
#   missing value handling.

merge_databases <- function(first_temporary_database, second_temporary_database) {
  merged_data <- second_temporary_database %>%
    left_join(first_temporary_database, by = "Target_ID") %>%
    dplyr::mutate(total_num_of_miRNA = ifelse(is.na(total_num_of_miRNA), 0, total_num_of_miRNA),
                  targeted_by_miRNA = ifelse(total_num_of_miRNA > 0, 1, 0),
                  all_miRNA_ID = ifelse(is.na(all_miRNA_ID), "No miRNA targeting", all_miRNA_ID))
  return(merged_data)
}


# Function Name: divide_data_into_groups
#
# This function divides the merged data into multiple groups based on adjusted p-values for protein and mRNA data.
# The groups are divided into background and test groups based on specific thresholds and conditions.
#
# Inputs:
# - merged_data: A dataframe containing the merged data from miRNA and protein/mRNA datasets.
# - adj_p_value_threshold: The threshold for the adjusted p-value to define significant data points.
#
# Output:
# - A list of dataframes containing the divided groups: background_group_1, background_group_2, 
#   test_group_1, test_group_2, and test_group_3.

divide_data_into_groups <- function(merged_data, adj_p_value_threshold) {
  # Background group 1
  background_group_1 <- merged_data %>%
    dplyr::filter(protein_adj_pvalue > adj_p_value_threshold & mRNA_adj_pvalue > adj_p_value_threshold)
  # Background group 2
  background_group_2 <- merged_data %>%
    dplyr::filter(protein_adj_pvalue <= adj_p_value_threshold & mRNA_adj_pvalue <= adj_p_value_threshold &
                    (protein_Abundance_log2_Ratio > 0 & mRNA_log2_FoldChange > 0 | protein_Abundance_log2_Ratio < 0 & mRNA_log2_FoldChange < 0))
  # Test group 1
  test_group_1 <- merged_data %>%
    dplyr::filter(protein_adj_pvalue <= adj_p_value_threshold & mRNA_adj_pvalue > adj_p_value_threshold)
  # Test group 2
  test_group_2 <- merged_data %>%
    dplyr::filter(protein_adj_pvalue > adj_p_value_threshold & mRNA_adj_pvalue <= adj_p_value_threshold)
  # Test group 3
  test_group_3 <- merged_data %>%
    dplyr::filter(protein_adj_pvalue <= adj_p_value_threshold & mRNA_adj_pvalue <= adj_p_value_threshold) %>%
    dplyr::filter((protein_Abundance_log2_Ratio > 0 & mRNA_log2_FoldChange < 0) | 
                    (protein_Abundance_log2_Ratio < 0 & mRNA_log2_FoldChange > 0))
  return(list(background_group_1 = background_group_1, 
              background_group_2 = background_group_2, 
              test_group_1 = test_group_1, 
              test_group_2 = test_group_2, 
              test_group_3 = test_group_3))
}


# Function Name: save_intermediate_data
#
# This function saves the intermediate data to an RData file for later use. 
# It provides a message confirming that the data has been saved.
#
# Inputs:
# - intermediate_data: The data to be saved, typically containing the intermediate results or processed data.
#
# Output:
# - The function saves the intermediate data to the file "intermediate_data.RData".

save_intermediate_data <- function(intermediate_data) {
  save(intermediate_data, file = "intermediate_data.RData")
  message("Intermediate data saved to 'intermediate_data.RData'.")
}


# Function Name: perform_statistical_tests
#
# This function performs statistical tests on the test groups against a background group to identify 
# overrepresented genes. The function iterates through different test conditions and calls the 
# 'find_over_represented_genes' function to obtain the results.
#
# Inputs:
# - test_groups_list: A list containing the test groups to be analyzed.
# - background_group: A dataframe containing the background group for comparison.
# - altered_miRNAs_type: The type of altered miRNAs to be considered (e.g., "UP", "DOWN").
# - miRNAs_threshold: The threshold for the number of miRNAs to target a gene.
# - about_test: A string providing additional information about the test.

perform_statistical_tests <- function(test_groups_list, background_group, altered_miRNAs_type, miRNAs_threshold, about_test) {
  results_list <- lapply(seq(1,3), function(t_num) {
    find_over_represented_genes(t_num = t_num,
                                test_groups_list = test_groups_list,
                                background_group = background_group,
                                bc_num = bc_num,
                                altered_miRNAs_type = altered_miRNAs_type,
                                miRNAs_threshold = miRNAs_threshold,
                                about_test = about_test)
  })
  
  return(results_list)
}


# Function Name: create_plots
#
# This function creates a set of plots based on the results of statistical tests. It arranges the individual plots 
# horizontally and returns a combined statistical plot.
#
# Inputs:
# - results_list: A list containing the results of the statistical tests, including data for plotting.
#
# Output:
# - A combined plot (statistical_plot) generated by arranging individual plots horizontally.

create_plots <- function(results_list) {
  plotlist <- lapply(seq(1,3), function(t_num) results_list[[t_num]][[paste0("test_", t_num, "_plot")]])
  statistical_plot <- ggpubr::ggarrange(plotlist = plotlist, nrow = 1)
  grid.newpage()
  grid.draw(statistical_plot)
  return(statistical_plot)
}


# Function Name: extract_p_values
#
# This function extracts the p-values from the results list of statistical tests. It returns a vector 
# containing all p-values associated with the test groups.
#
# Inputs:
# - results_list: A list containing the results of statistical tests, including p-value data for each group.
#
# Output:
# - A vector (p_values_vector) containing all extracted p-values.

extract_p_values <- function(results_list) {
  p_values <- lapply(seq(1,3), function(t_num) {
    test_p_value <- results_list[[t_num]][[paste0("test_", t_num, "_p_value")]]
    names(test_p_value) <- paste0("test_", t_num, "_p_value")
    return(test_p_value)
  })
  p_values_vector <- unlist(p_values)
  return(p_values_vector)
}


# Function Name: find_over_represented_genes
#
# This function identifies over-represented genes in the test groups compared to the background group using statistical tests. 
# It generates contingency tables and performs chi-square or Fisher tests to assess the association between the binary gene feature 
# (e.g., being targeted by miRNAs or having the KFERQ motif) and the group type (background vs. test group).
# The function also creates a mosaic plot to visually represent the results of the statistical test.
#
# Inputs:
# - t_num: Number indicating the test group to analyze (1, 2, or 3).
# - test_groups_list: List of dataframes for each test group.
# - background_group: Dataframe representing the background group for comparison.
# - bc_num: Number of the background group ("1" or "2").
# - altered_miRNAs_type: Type of altered miRNAs ("UP" or "DOWN").
# - miRNAs_threshold: Minimum number of miRNAs required to target a gene.
# - about_test: Type of test to perform ("regulated genes test" or "KFERQ motif test").
#
# Output:
# - A list containing:
#   1. p-value from the statistical test (p-value of the association between the gene feature and the group type).
#   2. A mosaic plot representing the results of the contingency table and statistical test.

find_over_represented_genes <- function(t_num,
                                        test_groups_list,
                                        background_group,
                                        bc_num,
                                        altered_miRNAs_type,
                                        miRNAs_threshold,
                                        about_test) {
  # 1) CHOOSE THE NAME OF THE BINARY COLUMN BASED ON THE PARAMETER ABOUT THE MEANING OF THE TEST
  if (about_test == "regulated genes test") {
    binary_col <- "targeted_by_miRNA"
  } else if (about_test == "KFERQ motif test") {
    binary_col <- "KFERQ_motif" 
  } else {
    msg <- "Error: about_test incorrect. Only \"regulated genes test\" and \"KFERQ motif test\" can be used."
    stop(msg)
  }
  # 2) GET THE BACKGROUND AND TEST GROUPS
  # gene background group
  background_db <- background_group %>%
    dplyr::select(!!rlang::sym(binary_col)) %>%
    dplyr::mutate(group = paste0("background\ngroup ", bc_num, "\n"),
                  !!rlang::sym(binary_col) := factor(ifelse(!!rlang::sym(binary_col) == "1", "Yes", "No"), levels = c("Yes", "No")))
  # gene test group
  test_db <- test_groups_list[[t_num]] %>%
    dplyr::select(!!rlang::sym(binary_col)) %>%
    dplyr::mutate(group = paste0("test\ngroup ", t_num, "\n"),
                  !!rlang::sym(binary_col) := factor(ifelse(!!rlang::sym(binary_col) == "1", "Yes", "No"), levels = c("Yes", "No")))
  # 3) MERGE THE TWO TABLES
  # Create a dataframe with data from the background and test group
  data_for_test <- rbind(background_db, test_db)
  # 4) IF THE MEANING OF THE TEST IS "regulated genes test" CHANGE THE LABELS OF THE PLOT. IF NOT, CONTINUE
  if (about_test == "regulated genes test") {
    # Based on the number of chosen miRNAs the colnames will be slightly different
    if (miRNAs_threshold == 1) {
      colnames(data_for_test) <- c(paste0("Genes targeted by\n",
                                          #miRNAs_threshold,
                                          #" ",
                                          altered_miRNAs_type,
                                          "-regulated miRNAs\n"),
                                   "group type\n\n")
    } else if(miRNAs_threshold > 1) {
      colnames(data_for_test) <- c(paste0("Genes targeted by at least ",
                                          #miRNAs_threshold,
                                          #" ",
                                          altered_miRNAs_type,
                                          "-regulated miRNAs\n"),
                                   "group type\n")
      
    } 
  } else if(about_test == "KFERQ motif test") {
    colnames(data_for_test) <- c(binary_col, "group type\n\n")
  } else {
    msg <- "ERROR: the number of miRNAs must be a positive number"
    stop(msg)
  }
  # 5) CREATE THE CONTINGENCY TABLE
  ct_table <- xtabs(~ ., data = data_for_test)
  # Print the observed value table
  message(paste0("Observed value table for test bc", bc_num, " vs ", "t_", t_num, ":"))
  print(ct_table)
  # 6) PERFORM THE STATISTICAL TEST
  statistical_test <- perform_statistical_test(ct_table = ct_table,
                                               bc_num = bc_num,
                                               t_num = t_num)
  # Extract p-value from the test result
  p_value <- statistical_test$p.value
  # 7) CREATE THE MOSAIC PLOTS
  plot_expr <- create_mosaic_plot(contingency_table = ct_table,
                                  background_group = background_group,
                                  bc_num = bc_num,
                                  test_group = test_group,
                                  t_num = t_num,
                                  statistical_test = statistical_test
  )
  # Add plot and p-value as a list
  result_list <- list(p_value = p_value, plot_expr = plot_expr)
  # Add the names to the list
  names(result_list) <- c(paste0("test_", t_num, "_p_value"),
                          paste0("test_", t_num, "_plot")
  )
  # Return the results
  return(result_list)
}


# Function Name: perform_statistical_test
#
# Description:
# This function performs a statistical test (either Chi-square or Fisher test) on a contingency table.
# Initially, it attempts to perform a Chi-square test to determine if there is an association between two categorical variables.
# If the Chi-square test raises warnings (likely due to small expected counts in the contingency table), the function automatically
# switches to performing a Fisher test, which is more reliable for small sample sizes and when expected frequencies are low.
#
# Inputs:
# - ct_table: A contingency table (created using xtabs() or similar methods), representing the observed counts of the categorical variables.
# - bc_num: Number of the background group ("1" or "2").
# - t_num: Number indicating the test group to analyze (1, 2, or 3).
#
# Output:
# - A statistical test result (either from the Chi-square test or Fisher test), including p-values, test statistics, and other relevant metrics.

perform_statistical_test <- function(ct_table,
                                     bc_num,
                                     t_num) {
  # 1) TRY TO PERFORM THE CHI-SQUARE TEST
  message(paste0("Performing Chi-square test for bg_", bc_num, " vs ", "t_", t_num, "..."))
  # Check if the chi-square test is giving some warnings/error and if so handle them
  warning_output <- FALSE
  tryCatch({
    statistical_test <- chisq.test(ct_table, correct = TRUE)
    print("Chi-square test performed")
  }, warning = function(w) {
    warning_output <- TRUE
  })
  # 2) IF THERE ARE SOME WARNINGS, PERFORM THE FISHER TEST
  if (isTRUE(warning_output)) {
    message(paste0("Not able to perform Chi-square test for bg_", bc_num, " vs ", "t_", t_num, "..."))
    # If the chi-square is raising a warning for the low observed values than perform fisher test
    statistical_test <- fisher.test(ct_table)
    message(paste0("Fisher test performed for bg_", bc_num, " vs ", "t_", t_num, "..."))
  }
  # Separate the messages of the different tests
  message("\n")
  return(statistical_test)
}


# Function Name: create_mosaic_plot
#
# Description:
# This function generates a mosaic plot visualizing the relationship between two categorical variables based on a contingency table.
# The plot shows the proportions of different combinations of the variables and highlights their association.
# The function also displays the p-value from a statistical test (Chi-square or Fisher) in the plot title.
#
# Libraries:
# - vcd: For creating mosaic plots.
# - grid: For advanced graphical functions and customization.
#
# Inputs:
# - contingency_table: A contingency table representing the observed counts of the categorical variables.
# - background_group: The background group identifier (e.g., "1" or "2").
# - bc_num: The background group number, used in the title and labeling.
# - test_group: The test group identifier (e.g., "1", "2", or "3").
# - t_num: The test group number, used in the title and labeling.
# - statistical_test: The result of a statistical test (either Chi-square or Fisher test), which contains the p-value to be displayed in the plot title.
#
# Output:
# - A mosaic plot showing the relationship between the two categorical variables, with labels indicating the observed counts and the p-value in the title.

create_mosaic_plot <- function(contingency_table, background_group, bc_num, test_group, t_num, statistical_test) {
  library(vcd)
  library(grid)
  # Create the mosaic plot
  plot_expr <- grid.grabExpr({
    mosaic(
      contingency_table,
      pop = FALSE,
      legend = TRUE,
      main = paste0(
        "P-value: ",
        signif(as.numeric(statistical_test$p.value), digits = 3),
        "\n\n\n\n"
      ),
      gp_main = gpar(fontsize = 40, col = "black", fontface = "bold"), # Customize title font size and style
      gp_labels = gpar(fontsize = 30, col = "black"), # Customize font size and color of side labels
      gp_varnames = gpar(fontsize = 30, col = "black", fontface = "bold"),# Customize legend font size
      margin = unit(c(5,5,5,5), "lines")
    )
    # Customize cell labels
    labeling_cells(
      text = contingency_table,
      gp_text = gpar(fontsize = 30, col = "black", fontface = "bold"), # Change cell font size and color
      margin = 0
    )(contingency_table) 
    # Return to the previous viewport
    popViewport()
  }, wrap.grobs = TRUE)  # Adjust the width and height for better fitting
  # Return the plot object
  return(plot_expr)
}