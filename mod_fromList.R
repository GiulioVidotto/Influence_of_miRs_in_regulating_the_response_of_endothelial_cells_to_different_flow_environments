# Function: mod_fromList
#
# This function takes a list of vectors and returns a binary data frame indicating 
# the presence or absence of each unique element across the input lists.
# Each row corresponds to a unique element (e.g., a gene, miRNA, etc.), and each 
# column corresponds to one of the input lists. The cell values are `1` if the 
# element is present in the respective list, and `0` if it is absent.
#
# The function works by:
#   1. Extracting all unique elements across the input lists.
#   2. Matching each element against every list and converting the result to binary 
#      (1 for presence, 0 for absence).
#   3. Reshaping the result into a data frame, where the row names are the unique elements 
#      and the column names are the names of the input lists.
#   4. Removing rows where all values are `0`, i.e., elements that are not present in any list.
#
# Input:
#   - input: A list where each element is a vector containing elements (e.g., genes, miRNAs).
#
# Output:
#   - A data frame with binary values (1 or 0) indicating the presence or absence of the 
#     unique elements in each list. Row names are the unique elements, and column names are 
#     the names of the input lists.

mod_fromList <- function (input) 
{
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F), stringsAsFactors = FALSE)
  row.names(data) <- elements
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  return(data)
}