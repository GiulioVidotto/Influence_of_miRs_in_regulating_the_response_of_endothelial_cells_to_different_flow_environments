# Function Name: transform_into_bed
#
# This function converts genomic data for a specific chromosome into a BED file-like format.
# It filters the input data for the specified chromosome, extracts relevant columns (start, stop,
# gene ID, score, and strand), and constructs a new data frame with the gene information and
# coordinates formatted as required for BED files.
#
# Inputs:
# - chr_num: The chromosome number to filter the genomic data. This can be numeric (e.g., "1") 
#            or a string (e.g., "X").
# - genomic_data: A data frame containing genomic information. It must include the following 
#                 columns:
#   - chr: Chromosome identifier.
#   - start: Start coordinate of the gene.
#   - stop: Stop coordinate of the gene.
#   - gene_ID: Unique identifier for the gene.
#   - score: Optional score for the region (numeric).
#   - strand: The strand information ("+" or "-").
#
# Outputs:
# - A data frame with one column (`gene_info`), where each row is formatted as:
#   chr_num <tab> start <tab> stop <tab> gene_ID <tab> score <tab> strand.

transform_into_bed <- function(chr_num, genomic_data) { 
  # Check for required columns
  required_columns <- c("chr", "start", "stop", "gene_ID", "score", "strand")
  missing_columns <- setdiff(required_columns, colnames(genomic_data))
  # If there are missing columns -> error
  if (length(missing_columns) > 0) {
    stop(paste("Error: Missing required columns in genomic_data:", 
               paste(missing_columns, collapse = ", ")))
  }
  # Filter only the coordinates of a specific chromosome
  chr_rows <- genomic_data %>% 
    dplyr::filter(chr == chr_num)
  # Take the strat and stop of each transcript and tne gene ID that is found at those coordinates
  start <- chr_rows$start
  stop <- chr_rows$stop
  gene_ID <- chr_rows$gene_ID
  score <- chr_rows$score
  strand <- chr_rows$strand
  # Create a new line with this format: chr_num:start-stop
  bed_format_info <- paste0(chr_num, "\t", start, "\t", stop, "\t", gene_ID, "\t", score, "\t", strand)
  # Create a new dataframe with the gene ID and the coordinates of that gene
  result <- data.frame(gene_info = bed_format_info,
                       stringsAsFactors = FALSE)
  # Return the results
  return(result)
}