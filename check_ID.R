# Function: check_ID (Quality control on miRNA notation)
#
# This function checks if the miRNA IDs in a specific database are updated to the latest version 
# of miRBase. If the notation differs, the function allows adding specifications to the miRNA_ID to 
# identify what part of the ID is missing (e.g., whether there is no specification for 5p or 3p).
# 
# The function first attempts to find perfect matches, and if some matches are not exact, it prompts 
# the user to decide whether to proceed with approximate matching. Approximate matching involves 
# adding additional information about the ID, specifying the differences compared to miRBase.
# 
# After running, the function provides a summary of how many miRNAs have a different nomenclature 
# compared to miRBase, specifying the differences (e.g., missing 5p/3p or other parts).
# 
# The function contains multiple sub-functions:
# - find_match: To find the matches between miRNA IDs
# - search_ID_column: To identify the column corresponding to the miRNA IDs
# - exact_match: To find the exact matches between the IDs
# - extract_ID: To extract miRNA IDs, as the match function returns indices
# - approx_match: To find approximate matches when exact ones are unavailable
# - obtain_updated_ID: To update approximate matches by adding missing details (e.g., 5p/3p)
# 
# Input: A database with standardized nomenclature (e.g., miRBase) and a query database with 
#        non-standardized nomenclature
# Output:
# - The function returns the modified query database, where the miRNA IDs are updated according 
#   to miRBase nomenclature. If approximate matches are found and updated, the IDs will be 
#   annotated with additional information (e.g., '-5p', '-3p', or specification of missing parts).
# - A summary is printed to the console, providing a breakdown of how many miRNA IDs had differences 
#   in their nomenclature compared to miRBase, and specifies what the differences are (e.g., missing 
#   5p/3p or other parts).
# - The summary is in the form of a named vector, with each name corresponding to a specific type of 
#   discrepancy (e.g., missing 5p, missing 3p, etc.), and each value indicating the number of IDs that 
#   exhibit that particular discrepancy.

check_ID <- function(std_db, query_db) {
  # Use the find_match function to compare miRNA IDs between the standard and query databases
  checked_ID <- find_match(std_db, query_db)
  # Create a vector that counts different categories of miRNA ID matches/mismatches
  vector_output <- c(
    sum(endsWith(checked_ID, ":'-5p' missing in miRBase")),
    sum(endsWith(checked_ID, ":'-3p' missing in miRBase")),
    sum(endsWith(checked_ID, ":'-3p' or '-5p' specification needed")),
    sum(endsWith(checked_ID, ":'-5p' specification needed")),
    sum(endsWith(checked_ID, ":'-3p' specification needed")),
    sum(endsWith(checked_ID, ": Not found")))
  # Assign descriptive names to the categories in the list_output vector
  names(vector_output) <- c(
    "'-5p' missing in miRBase",
    "'-3p' missing in miRBase",
    "'-3p' or '-5p' specification needed",
    "'-5p' specification needed",
    "'-3p' specification needed",
    "Not found")
  # Give to the user a summary of what has been found in the database
  print(vector_output)
  # Add the checked miRNA IDs to the query database for further analysis
  query_db <- query_db %>% 
    mutate(checked_miRNA_ID = checked_ID)
  # Return the updated query database
  return(query_db)
}


# Function: find_match
#
# Description:
# This function compares miRNA IDs from a standardized reference database (such as miRBase) with a query
# database to identify exact and approximate matches. The function first ensures that both ID lists are
# vectors and that all IDs are converted to lowercase. It also removes any version numbers from the miRNA IDs
# in the query database (if present). The function starts by checking for exact matches between the query and
# standard databases. If any IDs don't match exactly, the function prompts the user to decide whether to proceed
# with approximate matching. If the user agrees, the function performs approximate matching using a regular 
# expression pattern to find partial matches. Finally, the function returns a vector of the exact and/or approximate
# matches, providing a summary of the matching process to the user.
#
# Input:
# - std_db: A database containing standardized miRNA IDs (e.g., from miRBase).
# - query_db: A database with miRNA IDs that may not be standardized.
#
# Output:
# - A vector containing the miRNA IDs from the query database, either updated with exact or approximate matches
#   based on the comparison with the standard database (miRBase).

find_match <- function(std_db, query_db) {
  # Make sure that both ID lists are vectors and all letters are lowercase
  std_db_ID <- search_ID_column(std_db)
  query_db_ID <- search_ID_column(query_db)
  # If present, remove the version of the miRNA IDs
  query_db_ID <- if (any(grepl("\\.[1-9]+$", query_db_ID))) {
    sub("\\.[1-9]+$", "", query_db_ID)
  } else {
    query_db_ID
  }
  # Check exact matches
  exact_matches <- exact_match(query_db_ID, std_db_ID)
  # Check the number of NA (matches that are not exact matches)
  n_NA <- sum(is.na(exact_matches))
  n_extact_matches <- sum(!is.na(exact_matches))
  # Ask the user to continue with the approximate matching in case there are some IDs for which the exact match was not found
  if (n_NA > 0) {
    # Message to display to the user  
    message_to_display <- paste0(
      "Number of perfect matches found: ",
      n_extact_matches,
      "\n",
      "Number of IDs for which the perfect match was not found: ",
      n_NA,
      "\n",
      "Do you want to apply approximate matches? "
    )
    # Ask the user
    ANSWER <- switch(menu(c("Yes", "No"), title = message_to_display), "Yes", "No")
    # The answer could be "Yes" or "No". If another operation is written, the user will be asked again to write the answer
    if (ANSWER == "Yes") {
      # Create an index
      index <- seq(1, length(exact_matches))
      # set the progress bar (pb)
      pb <- txtProgressBar(min = 1,
                           max = length(index), 
                           initial = 1,
                           style = 3,
                           width = 50,
                           char = "=") 
      # Create the pattern to match based on the standard IDs
      pattern_to_match <- paste0(std_db_ID, ".*")
      # Function to get a list of approximate matches and exact matches
      approx_exact_matches <- unlist(map(index,
                                         approx_match,
                                         exact_matches = exact_matches,
                                         query_db_ID = query_db_ID,
                                         std_db_ID = std_db_ID,
                                         pattern_to_match = pattern_to_match,
                                         pb = pb
      )
      )
      # Close the progress bar
      close(pb)
      # Return the vector of approximate and exact matches
      return(approx_exact_matches)
    } else {
      return(exact_matches)
    }
  } else {
    # If everything is already updated and the number of NA values is equal to 0 return the list of exact matches
    message(paste0("All the IDs (", dim(query_db)[1], ") are already updated"))
    return(exact_matches)
  }
}


# Function: search_ID_column
# 
# Description:
# This function searches for the column in a given database (data frame) that contains miRNA IDs
# following a specific pattern. The pattern is defined as three letters, followed by a hyphen,
# another three letters, and a final alphanumeric sequence (e.g., `abc-def-123`). The function
# performs this search case-insensitively. It returns the values in the column that matches the pattern,
# providing a subset of the original database that corresponds to the miRNA IDs.
#
# Input:
# - database: A data frame containing various columns, one of which is expected to have miRNA IDs
#   with a specific pattern (3-letter - 3-letter - alphanumeric).
#
# Output:
# - A vector containing the miRNA IDs from the column that matches the defined pattern.
#   The function returns only the values of the column where the pattern is detected.

search_ID_column <- function(database) {
  # Searching for the ID column where the values has a specific pattern (3 letter - 3 letter - numbers/letters)
  # N.B the search is not sensible to upper and lower cases
  ID_column_True_False <- apply(database, 2, function (x) {stringr::str_detect(x, "^\\w{3}-\\w{3}-.*$")})
  ID_column <- as.vector(database[,ID_column_True_False])
  return(ID_column)
}


# Function: exact_match
#
# Description:
# This function performs an exact match check between two vectors of miRNA IDs: one from the query database
# and the other from the standard database (e.g., miRBase). The function uses the `match` function to find 
# the indices of the query database IDs that correspond to IDs in the standard database. It then converts 
# these indices into the actual miRNA IDs from the standard database using the `extract_ID` function. 
# The result is a vector containing the exact matches between the query and standard miRNA IDs.
#
# Input:
# - query_db_ID: A vector of miRNA IDs from the query database that need to be checked against the standard database.
# - std_db_ID: A vector of miRNA IDs from the standardized reference database (e.g., miRBase).
#
# Output:
# - A vector of the miRNA IDs from the standard database that correspond to exact matches found in the query database.

exact_match <- function(query_db_ID, std_db_ID) {
  # Inspect possible matches between the two vectors
  # N.B. The function match will return the index of the IDs in the standard vector that are found in the ID vector that needs to be updated)
  exact_matches <- match(query_db_ID, std_db_ID)
  # Transform each index into the actual ID
  exact_mir_ID <- unlist(map(exact_matches, extract_ID, std_db_ID = std_db_ID))
  # Return the output of the function exact_match
  return(extact_mir_ID)
}


# Function: extract_ID
#
# This function takes an index and extracts the corresponding miRNA ID from the 
# provided standard miRNA ID database. The index is used to retrieve the exact 
# match from the standard database vector.
#
# Input:
#   - index: The index corresponding to the miRNA ID to extract from the standard database (std_db_ID).
#   - std_db_ID: A vector of miRNA IDs from the standard miRNA database.
#
# Output:
#   - Returns the miRNA ID from the standard database corresponding to the provided index.

extract_ID <- function(index, std_db_ID) {
  elem <- std_db_ID[index]
  return(elem)
}


# Function: approx_match
#
# This function attempts to find approximate matches for miRNA IDs that do not have an exact match.
# If the ID at the given index has a perfect match, it returns that exact match. If the match is not 
# exact, the function attempts to find an approximate match based on a pattern and updates the ID accordingly.
#
# Input:
#   - index: The index of the current miRNA ID being checked.
#   - exact_matches: A vector containing the exact matches found previously.
#   - query_db_ID: The list of miRNA IDs from the query database that need to be matched.
#   - std_db_ID: The list of miRNA IDs from the standard database (e.g., miRBase).
#   - pattern_to_match: A regular expression pattern used to attempt matching with standard miRNA IDs.
#   - pb: A progress bar object to track the progress of the matching process.
#
# Output:
#   - Returns either the exact match if found, or an approximate match (updated miRNA ID) if no exact match is available.

approx_match <- function(index, exact_matches, query_db_ID, std_db_ID, pattern_to_match, pb) {
  
  # Check if the value is actually NULL, is should be NULL for the match to be a non perfect match
  if (!is.na(exact_matches[index])) {
    exact_or_updated_miRNA <- exact_matches[index]
    return(exact_or_updated_miRNA )
    # If the if-condition is FALSE it means that the match is not perfect and so approximate match can be tried
  } else {
    elem <- unlist(obtain_updated_ID(index = index,
                                     query_db_ID = query_db_ID,
                                     std_db_ID = std_db_ID,
                                     pattern_to_match = pattern_to_match,
                                     exact_matches = exact_matches,
                                     pb = pb
    )
    )
    return(elem)
  }
}

# Function: obtain_updated_ID
#
# This function attempts to find approximate matches for miRNA IDs that do not have an exact match in the standard reference database.
# If the ID at the given index has a perfect match in the reference database, it returns that exact match. 
# If the match is not exact, the function checks for possible discrepancies such as missing '-3p' or '-5p' annotations and updates the miRNA ID accordingly.
# Specifically, the function handles cases where:
#   - The query ID includes a '-5p' or '-3p' annotation not present in the reference database.
#   - The query ID is shorter and lacks a '-5p' or '-3p' annotation that exists in the reference database.
#   - No approximate match can be found.
#
# Input:
#   - index: The index of the current miRNA ID being checked within the query database.
#   - query_db_ID: A vector containing miRNA IDs from the query database that need to be matched.
#   - std_db_ID: A vector of miRNA IDs from the standard database (e.g., miRBase) to match against.
#   - exact_matches: A vector containing the exact matches found earlier (NA for unmatched IDs).
#   - pattern_to_match: A regular expression pattern used to attempt matching with standard miRNA IDs.
#   - pb: A progress bar object to track the progress of the matching process.
#
# Output:
#   - Returns either the exact match if found, or an updated miRNA ID with additional information (e.g., missing '-3p' or '-5p') if no exact match is found.

obtain_updated_ID <- function(index, query_db_ID, std_db_ID, exact_matches, pattern_to_match, pb) {
  # Check if the values if actually NA for this index
  if (is.na(exact_matches[index])) {
    # Based on the Na index, get the element in the query ID vector
    element_to_find <- query_db_ID[index]
    # Try to search if the element in the query database is longer than the one in miRBase (no 3p or 5p specification in miRBase)
    if (any(stringr::str_detect(element_to_find, pattern_to_match))) {
      end_specification <- ifelse(endsWith(element_to_find, "-5p"),
                                  ":'-5p' missing in miRBase", 
                                  ":'-3p' missing in miRBase")
      element_not_specified <- paste0(element_to_find, end_specification)
      setTxtProgressBar(pb,index)
      return(element_not_specified)
      # Try the search if the element in the query database is shorter (no 3p or 5p specification in the query database)
    } else if (any(grepl(paste0(element_to_find, "-"), pattern_to_match))) {
      true_false_matrix <- grepl(paste0(element_to_find, "-"), pattern_to_match)
      elements_std_ID <- std_db_ID[true_false_matrix]
      # Check if no elements are found
      if (length(elements_std_ID) == 0) {
        return(paste0(element_to_find, ": No matching elements in standard database"))
      }
      end_specification <- ifelse(length(elements_std_ID) == 2,
                                  ":'-3p' or '-5p' specification needed",
                                  ifelse(endsWith(elements_std_ID, "-5p"),
                                         ":'-5p' specification needed",
                                         ":'-3p' specification needed"))
      element_not_specified <- paste0(element_to_find, end_specification)
      setTxtProgressBar(pb,index)
      return(element_not_specified)
      # Otherwise no element is found
    } else {
      element_not_found <- paste0(element_to_find, ": Not found")
      setTxtProgressBar(pb,index)
      return(element_not_found)
    }
  } else {
    msg <- paste("the value at index", index, "is not Na")
    stop(msg, call. = FALSE)
  }
}