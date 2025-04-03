# Function: modified_plotSparsity
#
# This function modifies the original "plotSparsity" function by adding more details to enhance 
# the visualization of count sparsity in RNA-seq data. The function identifies problematic genes 
# where all counts are concentrated in a single sample and allows users to save a list of these 
# genes for further filtering.
#
# The function first calculates the sum of counts per gene and the proportion of counts coming 
# from the most abundant sample. Genes with all counts originating from one sample are labeled 
# as "Problematic genes." These genes are highlighted in the plot, and the user is given the option 
# to save a list of these genes for further analysis.
#
# After running, the function provides a summary of how many genes are considered problematic 
# and offers an option to store their IDs for filtering or removal.
#
# The function contains multiple sub-functions:
# - Compute gene-level statistics: Calculates the total sum of counts per gene and the proportion 
#   of counts coming from the most abundant sample.
# - Classify genes into categories: "Problematic genes" (counts concentrated in a single sample) 
#   and "Non-problematic genes" (distributed counts).
# - Generate an enhanced visualization: Uses log-scale transformation for better readability, 
#   and modifies color, shape, and size to highlight problematic genes.
# - Optionally save a list of problematic genes: If enabled, saves a vector of noisy gene IDs 
#   as "noisy_gene_vector.RData".
#
# Input: 
# - x: A DESeqDataSet object or a count matrix where rows represent genes and columns represent samples.
# - normalized (default = TRUE): Whether to use normalized counts (TRUE) or raw counts (FALSE).
# - save_intermediate (default = TRUE/FALSE): Whether to save a list of problematic genes for further analysis.
#
# Output:
# - A modified scatter plot that visualizes sparsity in RNA-seq data:
#   - X-axis: Sum of counts per gene (log scale).
#   - Y-axis: Proportion of counts from the most abundant sample.
#   - Problematic genes are colored red for easy identification.
# - A summary message specifying the number of problematic genes found.
# - A saved file ("noisy_gene_vector.RData") containing the IDs of problematic genes (if applicable).


# Modify the original "plotSparsity" function by addig more details
modified_plotSparsity <- function (x,
                                   normalized = TRUE,
                                   save_intermediate = c(TRUE, FALSE)) {
  if (is(x, "DESeqDataSet")) {
    x <- counts(x, normalized = normalized)
  }
  
  # Compute the values that will be plotted
  rs <- rowSums(x)
  rmx <- apply(x, 1, max)
  
  # Add these values to a data frame 
  df <- data.frame(
    sum_counts = rs[rs > 0], # sum of counts per gene
    max_div_sum = (rmx / rs)[rs > 0] # max count divided by the sum
  )
  
  # Convert the rownames into a column
  df$geneID <- rownames(df)
  
  df$combined <- ifelse(df$max_div_sum == 1, "Problematic genes", "Non problematic genes")
  
  # Create the plot
  sparity_plot <- ggplot(data = df, aes(x = sum_counts, y = max_div_sum)) +
    geom_point(aes(
      color = combined,  # Use the new combined column for color
      shape = combined,  # Use the new combined column for shape
      size = ifelse(max_div_sum == 1, 1, 0.25)  # Adjust size (bigger for problematic genes)
    )) +
    scale_color_manual(values = c("Problematic genes" = "red", 
                                  "Non problematic genes" = "black")) +  # Color mapping
    scale_shape_manual(values = c("Problematic genes" = 15, 
                                  "Non problematic genes" = 16)) +  # Shape mapping
    scale_size_identity() +  # Use identity mapping for size
    scale_x_log10() +  # Log scale for x-axis
    ylim(0, 1) +  # Y-axis limits
    labs(
      title = "Concentration of counts over total sum of counts",
      x = "Sum of counts per gene",
      y = "Max count / sum"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Center and enlarge title
      axis.title = element_text(size = 30),   # Increase axis title size
      axis.text = element_text(size = 30),    # Increase axis text size
      legend.position = "right",  # Adjust legend position
      legend.title = element_blank(),  # Remove legend title
      legend.text = element_text(size = 30)   # Increase legend text size
    ) +
    guides(
      size = guide_legend(override.aes = list(shape = 16, color = "black"))  # Adjust size legend
    )
  
  # Get the IDs of the genes with a max_div_sum equals to 1
  noisy_gene_vector <-  as.vector(df$geneID[df$max_div_sum == 1])
  
  # Starting from the count table, get the row counts of the nosy genes
  gene_count_table <- counts(paired_dds)
  noisy_gene_count_table <- gene_count_table[rownames(gene_count_table) %in% noisy_gene_vector,]
  
  # If there are noisy genes, save the vector of noisy genes (to be remove if the user wants to)
  if(save_intermediate == TRUE & length(noisy_gene_vector) > 0) {
    
    message("Noisy genes present in the data. The noisy gene vector will be saved to 'noisy_gene_vector.RData'")
    save(noisy_gene_vector, file = "noisy_gene_vector.RData")
    
    # Return the plot
    return(sparity_plot)
    
  } else if (save_intermediate == TRUE & length(noisy_gene_vector) == 0) {
    
    message("No noisy genes present in the data. The noisy gene vector won't be saved")
    
    # Return the plot
    return(sparity_plot)
    
  } else if (save_intermediate == FALSE & length(noisy_gene_vector) == 0) {
    
    message("No noisy genes present in the data.")
    
    # Return the plot
    return(sparity_plot)
    
  }
  
}