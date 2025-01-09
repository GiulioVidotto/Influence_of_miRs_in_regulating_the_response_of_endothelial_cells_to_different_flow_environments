# Function: find_target
#
# This function retrieves the top miRNA targets from the interaction dataset. 
# It filters interactions based on their database reliability (at least two databases) 
# and then calculates the number of targets for each miRNA across different database counts (4, 3, or 2 databases).
# The function further extracts and returns the specified number of top targets for each miRNA based on the n_target parameter.
#
# Input:
#   - miR_mRNA_interaction: A data frame containing miRNA-mRNA interaction data with columns 
#     `miRNA_ID`, `mRNA_ID`, and `number_of_database`, which indicates the number of databases the interaction is present in.
#   - miRNA_vector: A vector of miRNA IDs for which the target genes should be retrieved.
#   - n_target: The number of top targets to retrieve for each miRNA from the interaction data.
#
# Output:
#   - A list of miRNA-target interactions, where each list element corresponds to a specific miRNA ID and contains 
#     a list of target genes based on the specified number of targets (n_target).

find_common_target <- function (miR_mRNA_interaction, miRNA_vector, n_target) {
  
  # Take the miRNAs present in the interaction database
  miRNA_vector <- miRNA_vector[miRNA_vector %in% miR_mRNA_interaction$miRNA_ID]
  
  # Find the common targets
  common_target_list <- find_target(miR_mRNA_interaction = miR_mRNA_interaction,
                                    miRNA_vector = miRNA_vector,
                                    n_target = n_target)
  
  # Re-organise the data
  list_top_target_re_organized <- lapply(common_target_list, function(element) {re_organise_data(element)})
  
  # Obtain a list with the genes and the miRNAs that are targeting the same gene
  common_target <- get_common_target(list_top_target_re_organized)
  
  # Transform the list into a data frame
  df_common_target <- do.call(rbind, common_target)
  
  # Return the results as a data frame
  return(df_common_target)
  
}


# Function: find_common_target
#
# This function identifies the common target genes for a set of miRNAs from the provided interaction dataset. 
# It first filters the miRNA IDs from the query that are present in the interaction data. 
# Then, it finds the targets for each miRNA and re-organizes the data. 
# Finally, it identifies and returns the common target genes that are shared by the miRNAs in the input vector.
#
# Input:
#   - miR_mRNA_interaction: A data frame containing miRNA-mRNA interaction data with columns `miRNA_ID`, `mRNA_ID`, etc.
#   - miRNA_vector: A vector of miRNA IDs for which common target genes are to be identified.
#   - n_target: The number of top targets to retrieve for each miRNA from the interaction data.
#
# Output:
#   - A data frame containing the common target genes for the miRNAs provided in the input `miRNA_vector`. 
#     Each row represents a common target gene and the miRNAs that target it.

find_target <- function(miR_mRNA_interaction, miRNA_vector, n_target) {
  
  # Filtered out all the interactions present in only one database (they are not reliable)
  miR_mRNA_interaction <- miR_mRNA_interaction %>% 
    dplyr::filter(number_of_database > 1)
  
  # Create a dataset where there are the numbers of targets for each miRNA in different numbers of databases (in all 4, at least 3, at least 2)
  selected_miR_mRNA_interaction <- miR_mRNA_interaction %>% 
    dplyr::filter(miRNA_ID %in% miRNA_vector) %>%
    dplyr::group_by(miRNA_ID) %>% 
    dplyr::summarize(
      number_of_targets_in_4_databases = sum(number_of_database == 4),
      number_of_targets_in_3_databases = sum(number_of_database == 3),
      number_of_targets_in_2_databases = sum(number_of_database == 2)
    )
  
  # Print to the used this information (so the user can see the number of interactions for each group of databases)
  #print(selected_miR_mRNA_interaction)
  
  # Create a list with all the common targets and miRNAs that are targeting those genes
  list_miRNA_gene <- lapply(miRNA_vector, function(miRNA) {extract_targets(selected_miR_mRNA_interaction = selected_miR_mRNA_interaction,
                                                                           miR_mRNA_interaction = miR_mRNA_interaction,
                                                                           miRNA = miRNA,
                                                                           n_target = n_target
  )
  })
  
  return(list_miRNA_gene)
  
}


# Function: extract_targets
#
# This function retrieves the top target genes for a specific miRNA based on the number of databases that the interactions
# are found in (4, 3, or 2 databases). It determines how many target genes exist for the miRNA in each group of databases,
# and then it retrieves the top target genes according to the number of available databases.
#
# Input:
#   - selected_miR_mRNA_interaction: A data frame containing the number of targets for each miRNA in various databases (4, 3, 2 databases).
#     This data frame includes the `miRNA_ID` and columns representing the number of targets in 4, 3, and 2 databases.
#   - miR_mRNA_interaction: A data frame containing the full miRNA-mRNA interaction data.
#   - miRNA: The specific miRNA ID for which the target genes should be extracted.
#   - n_target: The number of top target genes to retrieve based on the available interactions in the databases.
#
# Output:
#   - A list of the top target genes for the specified miRNA. The targets are selected based on the number of databases 
#     they appear in, starting with 4 databases, then 3, and finally 2. The function attempts to gather the number of 
#     targets specified by `n_target`.

extract_targets <- function(selected_miR_mRNA_interaction, miR_mRNA_interaction, miRNA, n_target) {
  
  # Extract the row related to the specific miRNA
  one_miR_mRNA_interaction <- selected_miR_mRNA_interaction %>% 
    dplyr::filter(miRNA_ID == miRNA) %>% 
    dplyr::select(matches("number_of_targets_in_[1-9]_databases"))
  
  # Extract the number of target genes based on the numer of databases that they are found in
  n_4_databases <- one_miR_mRNA_interaction$number_of_targets_in_4_databases
  
  n_3_databases <- one_miR_mRNA_interaction$number_of_targets_in_3_databases
  
  n_2_databases <- one_miR_mRNA_interaction$number_of_targets_in_2_databases
  
  # Based on how many targets genes there are in each group of databases apply the function "target_gene" with different parameters
  if (n_4_databases >= n_target) {
    
    target_gene <- find_top_target(miR_mRNA_interaction = miR_mRNA_interaction,
                                   miRNA = miRNA,
                                   n_database = 4
    )
    
    return(target_gene)
    
  } else if ((n_4_databases + n_3_databases) >= n_target) {
    
    target_gene <- find_top_target(miR_mRNA_interaction = miR_mRNA_interaction,
                                   miRNA = miRNA,
                                   n_database = 3
    )
    
    return(target_gene)
    
  } else if ((n_4_databases + n_3_databases + n_2_databases) >= n_target) {
    
    target_gene <- find_top_target(miR_mRNA_interaction = miR_mRNA_interaction,
                                   miRNA = miRNA,
                                   n_database = 2
    )
    
    return(target_gene)
    
  }
  
}


# Function: find_top_target
#
# This function retrieves the target genes for a specific miRNA based on the number of databases in which the interaction
# is found. The target genes are filtered by the minimum number of databases specified by the `n_database` parameter.
# 
# Input:
#   - miR_mRNA_interaction: A data frame containing the miRNA-mRNA interaction data, with columns such as `miRNA_ID` 
#     and `number_of_database` (indicating the number of databases in which the interaction was found).
#   - miRNA: The specific miRNA ID for which the target genes should be retrieved.
#   - n_database: The minimum number of databases that the interaction must be found in for the target gene to be included.
#
# Output:
#   - A data frame containing the target genes for the specified miRNA, where the interactions are found in at least `n_database` databases.

find_top_target <- function(miR_mRNA_interaction, miRNA, n_database) {
  
  # Find the targets genes based on the number of databases and the name of the miRNA
  gene <- miR_mRNA_interaction %>% 
    dplyr::filter(miRNA_ID == miRNA & number_of_database >= n_database)
  
  # Return the results
  return(gene)
  
}


# Function: re_organise_data
#
# This function reorganizes the input data by creating a new column that concatenates the `Target_ID` and `gene_name`
# separated by a colon. The function also renames the resulting column to reflect the specific miRNA ID in the input element.
#
# Input:
#   - element: A data frame containing columns `Target_ID` and `gene_name`, as well as a column for `miRNA_ID`. 
#     The function expects that the element corresponds to a single miRNA, and each row contains interactions with target genes.
#
# Output:
#   - A data frame with a new column `Target_ID_gene_name`, where the `Target_ID` and `gene_name` are concatenated by a colon.
#     The function renames this column to reflect the miRNA ID for that element.

re_organise_data <-function(element) {
  
  # for each element get only the target_ID and the gene name separated by ":"
  element_result <- element %>% 
    dplyr::mutate(Target_ID_gene_name = paste0(Target_ID, ":", gene_name)) %>% 
    dplyr::select(Target_ID_gene_name)
  
  
  # Add column names based on the miRNA ID and on the number of databases were that interaction was found
  if (length(unique(element$miRNA_ID)) == 1) {
    
    colnames(element_result) <- paste0(unique(element$miRNA_ID), ":")
    
  } else {
    
    mgs <- "There should be only one miRNA for list element"
    stop(msg)
    
  }
  
  # Return the output
  return(element_result)
  
}


# Function: get_common_target
#
# This function identifies the common target genes that are targeted by multiple miRNAs. It takes as input 
# a list of reorganized target data (miRNAs and their target genes) and returns a list of data frames with 
# information about the common target genes, their corresponding ENSEMBL ID, gene name, and the miRNAs that target them.
#
# Input:
#   - list_top_target_re_organized: A list of data frames, where each element contains a `Target_ID_gene_name` column,
#     which concatenates `Target_ID` and `gene_name` for miRNAs and their target genes.
#
# Output:
#   - A list of data frames where each data frame contains information about one common target gene. For each common gene, 
#     the data frame includes:
#     - `common_target_ENSEMBL_ID`: The ENSEMBL ID of the gene.
#     - `common_target_name`: The name of the gene.
#     - `miRNA_ID`: A vector of miRNAs targeting the gene.
#     - `all_miRNA_ID`: A comma-separated string of all miRNAs targeting the gene.
#     - `total_num_of_miRNA`: The total number of miRNAs targeting the gene.

get_common_target <- function(list_top_target_re_organized) {
  
  # Create a vector with all the genes
  gene_vector <- unlist(list_top_target_re_organized)
  
  # Find the genes that are repeated at least 2 times (the common target genes)
  common_target <- unique(gene_vector[duplicated(gene_vector)])
  
  # For each duplicated gene, find which are the miRNAs that are targeting that gene
  final_result <- lapply(unique(gene_vector), function(gene) {
    
    # Create a vector of true and false depending on where a specific gene is found in the dataframe, where each column is a list of genes targeted by a specific miRNA
    TRUE_FALSE_vector <- grepl(gene, gene_vector)
    
    # Based on the true and false vector, get the names of the miRNAs
    miRNA_info_vector <- names(gene_vector[TRUE_FALSE_vector])
    
    # This vector contains: miRNA_name
    unique_miRNA_vector <- unique(sub(":\\d*$", "", miRNA_info_vector))
    
    # Get the ENSEMBL ID and the gene name
    gene_info <- unlist(strsplit(gene, ":"))
    ENSEMBL_ID <- gene_info[1]
    gene_name <- gene_info[2]
    
    # Build a new dataframe with the name of the genes and the names of the miRNAs that are targeting that gene
    result <- data.frame(common_target_ENSEMBL_ID = ENSEMBL_ID,
                         common_target_name = gene_name,
                         miRNA_ID = unique_miRNA_vector,
                         all_miRNA_ID = paste0(unique_miRNA_vector, collapse = ", "),
                         total_num_of_miRNA = length(unique_miRNA_vector),
                         stringsAsFactors = FALSE)
    
    # Return the output
    return(result)
    
  })
  
  # Return the output
  return(final_result)
  
}