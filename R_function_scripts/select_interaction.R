# Function: select_interaction
#
# This function extracts interactions that are common across different databases based on user-specified settings. 
# The function also utilizes a helper function `split_into_columns` to separate miRNA IDs from target gene IDs.
#
# Input:
#   - data: A data frame (e.g., from an upset plot) containing interaction information across multiple databases.
#   - number: The number of databases to consider for identifying common interactions.
#   - comparison: A string that defines the type of comparison to make (e.g., "equal", "greater or equal", "greater", "lower or equal", "lower").
#
# Output:
#   - A data frame with miRNA IDs, target gene IDs, and the corresponding database information for the common interactions.

# Function to get the desired interaction from the upset plot
select_interaction <- function(data, number, comparison) {
  # Get the names of the interactions
  interaction_vector <- rownames(data)
  # Number of total databases to consider
  n_databases <- length(colnames(data)[colnames(data) != c("source", "row_sum")])
  if(number > n_databases) {
    stop("Error: the number selected is higher than the considered databases") 
  } else {
    if (comparison == "equal") {
      # Get the common interactions among the selected databases (equal)
      common_interaction <- interaction_vector[data$row_sum == number]
      source <- data$source[data$row_sum == number]
      number_of_database <- data$row_sum[data$row_sum == number]
      info_database <- paste0(common_interaction, "_", source, "_", number_of_database)
    } else if(comparison == "greater or equal") {
      # Get the common interactions among the selected databases (greater or equal)
      common_interaction <- interaction_vector[data$row_sum >= number]
      source <- data$source[data$row_sum >= number]
      number_of_database <- data$row_sum[data$row_sum >= number]
      info_database <- paste0(common_interaction, "_", source, "_", number_of_database)
    } else if(comparison == "greater") {
      # Get the common interactions among the selected databases (greater)
      common_interaction <- interaction_vector[data$row_sum > number]
      source <- data$source[data$row_sum > number]
      number_of_database <- data$row_sum[data$row_sum > number]
      info_database <- paste0(common_interaction, "_", source, "_", number_of_database)
    } else if(comparison == "lower or equal") {
      # Get the common interactions among the selected databases (lower)
      common_interaction <- interaction_vector[data$row_sum < number]
      source <- data$source[data$row_sum < number]
      number_of_database <- data$row_sum[data$row_sum < number]
      info_database <- paste0(common_interaction, "_", source, "_", number_of_database)
    } else if(comparison == "lower") {
      # Get the common interactions among the selected databases (lower or equal)
      common_interaction <- interaction_vector[data$row_sum <= number]
      source <- data$source[data$row_sum <= number]
      number_of_database <- data$row_sum[data$row_sum <= number]
      info_database <- paste0(common_interaction, "_", source, "_", number_of_database)
    } else {
      msg <- paste0(comparison, " not valid\n", "Select between: equal, greater or equal, greater, lower or equal, lower")
      stop(msg)
    }
    # Once the interactions are obtained, split them between miRNA IDs and Target gene IDs
    info_database_split <- strsplit(info_database, "_")
    # Get the single elements from the list previously obtained (common_interaction_split)
    miRNA_ID <- unlist(lapply(info_database_split, split_into_columns, 1))
    Target_ID <- unlist(lapply(info_database_split, split_into_columns, 2))
    source_database <- unlist(lapply(info_database_split, split_into_columns, 3))
    number_of_database <- unlist(lapply(info_database_split, split_into_columns, 4))
    # Create a dataframe with the miRNA IDs and the target gene IDs
    results <- data.frame(miRNA_ID, Target_ID, source_database, number_of_database, stringsAsFactors = FALSE)
    # Return the results
    return(results)
  }
}

# helper function
# Function to divide the miRNA IDs from the target gene IDs
split_into_columns <- function(list, index) {
  return(list[index])
}