# Function Name: update_stop_coord
#
# Description:
# This function processes the output from APAtrap to extract and clean data regarding transcript IDs,
# positions of predicted APA (Alternative Polyadenylation) sites, and coordinates of 3' UTR regions.
# It reformats the data by removing invalid entries, redefining column names, and extracting relevant
# genomic coordinates (chromosome, strand, start, and stop positions).
#
# Inputs:
# - file_path: Path to a file containing the APAtrap output. The file should be a comma-separated text
#              file with headers, containing predicted APA site information.
#
# Outputs:
# - A cleaned data frame containing the following columns:
#   - Transcript_ID: Unique identifier for the transcript.
#   - MSE: Mean squared error for the APA prediction.
#   - chr: Chromosome on which the transcript is located.
#   - strand: Strand information for the transcript ("+" or "-").
#   - start: Start position of the 3' UTR.
#   - stop: Stop position of the 3' UTR.

update_stop_coord <- function(file_path) {
  # Check that the file exists
  if (!file.exists(APA_output_path)) {
    stop("Error: APA output file does not exist at the specified path.")
  } else {
    # Upload the data
    poly_A_site_db <- read.table(file_path, sep = ",", stringsAsFactors = FALSE, header = TRUE)
  }
  # Checks on the columns
  required_columns <- c("Transcript_ID", "MSE", "UTR_coordinates")
  if (!all(required_columns %in% colnames(poly_A_site_db))) {
    stop("Error: Missing required columns in the input file. Ensure it contains 'Transcript_ID', 'MSE', and 'UTR_coordinates'.")
  }
  # Check on the file
  if (nrow(poly_A_site_db) == 0) {
    stop("Error: The input file is empty.")
  }
  # Take only the infomation about the transcript_ID, position of the predicted APA sites and the coordinates of the 3'UTR region
  poly_A_site_db <- poly_A_site_db %>% 
    dplyr::filter(Gene != "Gene") # There are some lines where there is "Gene" as gene field (this is due to using cat and these are the headers of all the predictAPA output files)
  # Re-define the colnames
  colnames(poly_A_site_db) <- c("Transcript_ID", "MSE", "UTR_coordinates")
  # Extract the coordinates of the furthest poly-A site
  poly_A_site_db <- poly_A_site_db %>%
    dplyr::mutate(
      chr = gsub(".*\\|NA\\||\\|[+-]", "", Transcript_ID),
      strand = gsub(".*\\|chr.*\\|", "", Transcript_ID),
      Transcript_ID = gsub("\\|.*", "", Transcript_ID),
      start = gsub("chr.*:|-.*", "", UTR_coordinates),
      stop = gsub(".*-", "", UTR_coordinates)) %>% 
    dplyr::filter(!is.na(chr)) %>% 
    dplyr::select(-UTR_coordinates)
  # Return the results
  return(poly_A_site_db)
}
