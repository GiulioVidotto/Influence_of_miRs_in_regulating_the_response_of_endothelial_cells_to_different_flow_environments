# Function: get_MFE
#
# This function processes a file line to extract the gene ID, miRNA ID, and the Minimum Free Energy (MFE) value.
# It identifies lines that represent a gene ID, a miRNA ID, or an MFE value and extracts the relevant information.
# 
# The function operates as follows:
# 1. If the line contains a gene ID (starting with ">ENSG"), it extracts and stores the gene ID.
# 2. If the line contains a miRNA ID (starting with ">hsa-"), it extracts and stores the miRNA ID.
# 3. If the line contains an MFE value in a specific format, it extracts the value and returns it along with the gene ID and miRNA ID in a data frame.
# 
# Input:
#   - file_line: A string containing a line from a file. The line may contain:
#     - A gene ID starting with ">ENSG".
#     - A miRNA ID starting with ">hsa-".
#     - A string representing an MFE value in the format of "(-10.5, 0.5)" or similar.
#
# Output:
#   - If the line contains an MFE value, a data frame with the following columns is returned:
#     - gene_ID: The gene identifier (ENSG format).
#     - miRNA_ID: The miRNA identifier (hsa- format).
#     - MFE_value: The Minimum Free Energy (MFE) value extracted from the line.
#   - If the line contains only gene or miRNA information, no output is returned but the corresponding ID is stored for later use.

get_MFE <- function(file_line) {
  if (startsWith(file_line, ">ENSG")) {
    gene_ID <<- gsub(">ENSG","ENSG", file_line)
    NULL
  } else if (startsWith(file_line, ">hsa-")) {
    miRNA_ID <<- gsub(">hsa-", "hsa-", file_line)
    NULL
  } else {
    MFE_value <- gsub("[(\\.]*&[)\\.]*\\s*\\d*,\\d*\\s*:\\s*\\d*,\\d*\\s*\\(|\\s+.*", "", file_line)
    result <- data.frame(gene_ID = gene_ID,
                         miRNA_ID = miRNA_ID,
                         MFE_value = MFE_value,
                         stringsAsFactors = FALSE)
    return(result)
  }
}