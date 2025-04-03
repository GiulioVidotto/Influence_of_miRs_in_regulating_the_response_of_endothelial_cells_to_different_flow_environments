# mRNA-seq Dataset Overview
#
# - Name of the file -> "mRNA_expression_over_0.xlsx" (contains transcript expression levels across different conditions)
#
# - DESCRIPTION:
# This mRNA-seq dataset contains quantitative information on transcript expression levels.
# While the other two datasets obtained from the lab (proteomics and miRNAs datasets) were already processed
# we processed the mRNA-seq data starting from the row counts of the transcripts following the Deseq2 pipeline.

# --- 1. Get Working Directory ---
# Set the working directory
current_working_directory <- getwd()

# --- 2. Import the Database ---
# Define the path to the mRNA-seq data file
original_proteomics_data_path <- paste0(current_working_directory, "project_data/lab_data/mRNA_seq_database/mRNA_expression_over_0.xlsx")

# Check if the file exists before importing
if (!file.exists(original_proteomics_data_path)) {
  stop("Error: The file was not found in the specified directory.")
}

# Import the proteomics data in R
mRNA_expression_database <- readxl::read_excel(path_mRNA_expression_database,
                                               na = "")

# --- 3. Operations on the dataset before performing the DESeq2 analysis ---
# Add the rownames to the dataset
mRNA_expression_database <- as.data.frame(mRNA_expression_database)
rownames(mRNA_expression_database) <- gsub("\\.\\d*", "", mRNA_expression_database$geneID)

# Extract the row counts from the matrix (gene IDs as rownames)
count_matrix <- mRNA_expression_database %>% 
  dplyr::select(patient_1_OSS = SW1_S37_rawCounts_AL4.19_OSS,
                patient_1_LSS = SW2_S38_rawCounts_AL4.10_LSS,
                patient_1_ESS = SW3_S39_rawCounts_AL4.1_ESS,
                patient_2_OSS = SW4_S40_rawCounts_BP5.19_OSS,
                patient_2_LSS = SW5_S41_rawCounts_BP5.13_LSS,
                patient_2_ESS = SW6_S42_rawCounts_BP5.1_ESS,
                patient_3_OSS = SW7_S43_rawCounts_BP6.15_OSS,
                patient_3_LSS = SW8_S44_rawCounts_BP6.8_LSS,
                patient_3_ESS = SW9_S45_rawCounts_BP6.1_ESS
  )

# --- 4. DESeq2 analysis ---

# --- 4.1 Check the order of the columns and rows ---
# Create a dataframe with the column data information
coldata <- data.frame(condition = factor(gsub("patient_[0-9]+_", "", colnames(count_matrix))),
                      donor = factor(gsub("patient_|_[A-Z]*", "", colnames(count_matrix)))
)

# Add rownmaes to this dataset
rownames(coldata) <- colnames(count_matrix)

# Check that the order of the rows of the column data and of the columns of the count matrix is the same
all(rownames(coldata) == colnames(count_matrix)) # TRUE

# --- 4.1 Construct the DESeqDataSet object ---

# Create the DESeqDataSet object
# Paired (condition and donor)
paired_dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                     colData = coldata,
                                     design = ~ condition + donor)

# --- 4.2 Pre-filtering ---
smallestGroupSize <- 1
paired_keep <- rowSums(counts(paired_dds) >= 15) >= smallestGroupSize
paired_dds <- paired_dds[paired_keep,]

paired_dds$condition <- relevel(paired_dds$condition, ref = "LSS")

# --- 4.3 Plot the counts before and after normalization ---
# Counts before normalization
# Extract the counts which are not normalized yet
not_normalised_counts <- counts(paired_dds)

# Apply a log transformation to these counts
log_not_normalised_counts <- log2(not_normalised_counts + 1)

# Convert to a long-format data frame
log_not_normalised_counts <- as.data.frame(log_not_normalised_counts) 
  pivot_longer(cols = c(patient_1_OSS, patient_1_LSS, patient_1_ESS, patient_2_OSS, patient_2_LSS, patient_2_ESS, patient_3_OSS, patient_3_LSS, patient_3_ESS),
               names_to = "patient_condition",
               values_to = "Values")

# Plot the counts before normalization
ggplot(log_not_normalised_counts, aes(y = Values, x = as.factor(patient_condition), fill = as.factor(patient_condition))) +
  geom_boxplot(outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Samples", y = "Log2(Normalized Counts)") +
  scale_fill_manual(values = c("patient_1_ESS" = "#B3E0A3",
                               "patient_1_OSS" = "#B3C6FF",
                               "patient_1_LSS" = "#FFB3BA",
                               "patient_2_ESS" = "#B3E0A3",
                               "patient_2_OSS" = "#B3C6FF",
                               "patient_2_LSS" = "#FFB3BA",
                               "patient_3_ESS" = "#B3E0A3",
                               "patient_3_OSS" = "#B3C6FF",
                               "patient_3_LSS" = "#FFB3BA"
  ),
  labels = c("patient_1_ESS" = "Patient 1 - ESS",
             "patient_1_OSS" = "Patient 1 - OSS",
             "patient_1_LSS" = "Patient 1 - LSS",
             "patient_2_ESS" = "Patient 2 - ESS",
             "patient_2_OSS" = "Patient 2 - OSS",
             "patient_2_LSS" = "Patient 2 - LSS",
             "patient_3_ESS" = "Patient 3 - ESS",
             "patient_3_OSS" = "Patient 3 - OSS",
             "patient_3_LSS" = "Patient 3 - LSS")
  ) +
  scale_x_discrete(labels = c("patient_1_ESS" = "Patient 1 - ESS",
                              "patient_1_OSS" = "Patient 1 - OSS",
                              "patient_1_LSS" = "Patient 1 - LSS",
                              "patient_2_ESS" = "Patient 2 - ESS",
                              "patient_2_OSS" = "Patient 2 - OSS",
                              "patient_2_LSS" = "Patient 2 - LSS",
                              "patient_3_ESS" = "Patient 3 - ESS",
                              "patient_3_OSS" = "Patient 3 - OSS",
                              "patient_3_LSS" = "Patient 3 - LSS"
  )) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Center and enlarge title
    axis.title = element_text(size = 30),   # Increase axis title size
    axis.text = element_text(size = 40),    # Increase axis text size
    legend.position = "none"
    #legend.title = element_text(size = 30), # Increase legend title size
    #legend.text = element_text(size = 30)   # Increase legend text size
  ) +
  guides(fill = guide_legend(title = NULL))

# Calculate the normalized counts by calculating the estimate size factors
paired_dds <- estimateSizeFactors(paired_dds)
normalized_counts <- counts(paired_dds, normalized = TRUE)

# Apply a log transformation to these counts and convert them to a long-format data frame
log_normalized_counts <- as.data.frame(log2(normalized_counts + 1)) %>% 
  pivot_longer(cols = c(patient_1_OSS, patient_1_LSS, patient_1_ESS, patient_2_OSS, patient_2_LSS, patient_2_ESS, patient_3_OSS, patient_3_LSS, patient_3_ESS),
               names_to = "patient_condition",
               values_to = "Values")

# Plot the counts after normalization
ggplot(log_normalized_counts, aes(y = Values, x = as.factor(patient_condition), fill = as.factor(patient_condition))) +
  geom_boxplot(outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Samples", y = "Log2(Normalized Counts)") +
  scale_fill_manual(values = c("patient_1_ESS" = "#B3E0A3",
                               "patient_1_OSS" = "#B3C6FF",
                               "patient_1_LSS" = "#FFB3BA",
                               "patient_2_ESS" = "#B3E0A3",
                               "patient_2_OSS" = "#B3C6FF",
                               "patient_2_LSS" = "#FFB3BA",
                               "patient_3_ESS" = "#B3E0A3",
                               "patient_3_OSS" = "#B3C6FF",
                               "patient_3_LSS" = "#FFB3BA"
  ),
  labels = c("patient_1_ESS" = "Patient 1 - ESS",
             "patient_1_OSS" = "Patient 1 - OSS",
             "patient_1_LSS" = "Patient 1 - LSS",
             "patient_2_ESS" = "Patient 2 - ESS",
             "patient_2_OSS" = "Patient 2 - OSS",
             "patient_2_LSS" = "Patient 2 - LSS",
             "patient_3_ESS" = "Patient 3 - ESS",
             "patient_3_OSS" = "Patient 3 - OSS",
             "patient_3_LSS" = "Patient 3 - LSS")
  ) +
  scale_x_discrete(labels = c("patient_1_ESS" = "Patient 1 - ESS",
                              "patient_1_OSS" = "Patient 1 - OSS",
                              "patient_1_LSS" = "Patient 1 - LSS",
                              "patient_2_ESS" = "Patient 2 - ESS",
                              "patient_2_OSS" = "Patient 2 - OSS",
                              "patient_2_LSS" = "Patient 2 - LSS",
                              "patient_3_ESS" = "Patient 3 - ESS",
                              "patient_3_OSS" = "Patient 3 - OSS",
                              "patient_3_LSS" = "Patient 3 - LSS"
  )) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Center and enlarge title
    axis.title = element_text(size = 30),   # Increase axis title size
    axis.text = element_text(size = 40),    # Increase axis text size
    legend.position = "none"
    #legend.title = element_text(size = 30), # Increase legend title size
    #legend.text = element_text(size = 30)   # Increase legend text size
  ) +
  guides(fill = guide_legend(title = NULL))

# --- 4.4 Differential expression analysis (DEA) on paired data ---
# Performed the DEA on the paired data
paired_dds <- DESeq(paired_dds)

# --- 4.5 Quality check: PCA plot ---
# Perform regularized log transformation (rlog)
rld <- rlog(paired_dds, blind = FALSE)

# Extract PCA data with grouping variables
pcaData <- plotPCA(rld, intgroup=c("condition", "donor"), returnData=TRUE)

# Calculate percentage variance explained by each principal component
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Plot the PCA plot to understand the distance between samples
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=donor)) +
  geom_point(size=10) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  #labs(title = "Principal Component Analysis (PCA) biplot") +
  #theme_minimal() +
  theme(
    #plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Center and enlarge title
    axis.title = element_text(size = 30),   # Increase axis title size
    axis.text = element_text(size = 30),    # Increase axis text size
    legend.title = element_text(size = 30), # Increase legend title size
    legend.text = element_text(size = 30)   # Increase legend text size
  )

# --- 4.6 Quality check: Sparsity plot ---
# Plot the sparsity plot with the modified function
modified_plotSparsity(paired_dds, normalized = TRUE, save_intermediate = TRUE)

# Load the noisy gene vector
load("noisy_gene_vector.RData")

# Some genes highlighted in red could be a problem/ noisy. Remove them
paired_dds <- paired_dds[!rownames(paired_dds) %in% noisy_gene_vector,]

# --- 4.7 OSS vs LSS contrast results ---
# Results of the OSS vs LSS contrast
paired_res_OSS_vs_LSS <- results(paired_dds,
                                 contrast = c("condition", "OSS", "LSS"),
                                 cooksCutoff = TRUE,
                                 independentFiltering = TRUE)

# Add the rownames as a column
paired_res_OSS_vs_LSS$geneID <- rownames(paired_res_OSS_vs_LSS)

# Transfrom the table into a tibble
paired_res_OSS_vs_LSS <- as_tibble(paired_res_OSS_vs_LSS)

# Extract only the columns about the fold changes and the p-values with the associated gene IDs
paired_res_OSS_vs_LSS <- paired_res_OSS_vs_LSS %>% 
  dplyr::select(geneID,
                OSS_vs_LSS_paired_log2FoldChange = log2FoldChange,
                OSS_vs_LSS_paired_pvalue = pvalue,
                OSS_vs_LSS_paired_padj = padj,
                OSS_vs_LSS_paired_baseMean = baseMean)

# --- 4.8 Quality check: distribution of the p-values ---
# Plot an histogram of the original p-values
OSS_vs_LSS_paired_p_value_histogram <- ggplot(data = paired_res_OSS_vs_LSS) +
  geom_histogram(aes(x = OSS_vs_LSS_paired_pvalue), binwidth = 0.005) +
  labs(title = "OSS vs LSS - paired p-value",
       x = "p-value",
       y = "frequency") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Center and enlarge title
    axis.title = element_text(size = 14),   # Increase axis title size
    axis.text = element_text(size = 12),    # Increase axis text size
    legend.position = "right",  # Adjust legend position
    legend.title = element_blank(),  # Remove legend title
    legend.text = element_text(size = 10)   # Increase legend text size
  )

show(OSS_vs_LSS_paired_p_value_histogram)

# The shape of the histogram is similar to the "Scenario A: Anti-conservative p-values" in the article "http://varianceexplained.org/statistics/interpreting-pvalue-histogram/"
# This means that there is no evident issue with the p-values and it is right to apply the correction.

# --- 4.9 ESS vs LSS contrast results ---
# Results for ESS vs LSS contrast
paired_res_ESS_vs_LSS <- results(paired_dds,
                                 contrast = c("condition", "ESS", "LSS"),
                                 cooksCutoff = TRUE,
                                 independentFiltering = TRUE)

# Add the rownames as a column
paired_res_ESS_vs_LSS$geneID <- rownames(paired_res_ESS_vs_LSS)

# Transfrom the table into a tibble
paired_res_ESS_vs_LSS <- as_tibble(paired_res_ESS_vs_LSS)

# Extract only the columns about the fold changes and the p-values with the associated gene IDs
paired_res_ESS_vs_LSS <- paired_res_ESS_vs_LSS %>% 
  dplyr::select(geneID,
                ESS_vs_LSS_paired_log2FoldChange = log2FoldChange,
                ESS_vs_LSS_paired_pvalue = pvalue,
                ESS_vs_LSS_paired_padj = padj,
                ESS_vs_LSS_paired_baseMean = baseMean)

# --- 4.10 Quality check: distribution of the p-values ---
# Plot an histogram of the original p-values
ESS_vs_LSS_paired_p_value_histogram <- ggplot(data = paired_res_ESS_vs_LSS) +
  geom_histogram(aes(x = ESS_vs_LSS_paired_pvalue), binwidth = 0.005) +
  labs(title = "ESS vs LSS - paired p-value",
       x = "p-value",
       y = "frequency") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Center and enlarge title
    axis.title = element_text(size = 14),   # Increase axis title size
    axis.text = element_text(size = 12),    # Increase axis text size
    legend.position = "right",  # Adjust legend position
    legend.title = element_blank(),  # Remove legend title
    legend.text = element_text(size = 10)   # Increase legend text size
  )

show(ESS_vs_LSS_paired_p_value_histogram)

# The shape of the histogram is similar to the "Scenario A: Anti-conservative p-values" in the article http://varianceexplained.org/statistics/interpreting-pvalue-histogram/
# This means that there is no evident issue with the p-values and it is right to apply the correction.

# --- 4.11 Creation of a table with all the information from the DESeq2 analysis ---
# Normalised counts for paired dds object
paired_dds_normalized_counts <- as.data.frame(counts(paired_dds, normalized=TRUE)) %>% 
  rownames_to_column(var = "geneID")

# Merge both tables for ESS and LSS. Add the information about the fold changes and p-values to the counts
mRNA_expresion_deseq2_table <- paired_dds_normalized_counts %>% 
  left_join(paired_res_ESS_vs_LSS, by = "geneID") %>% 
  left_join(paired_res_OSS_vs_LSS, by = "geneID")

# There are rows where there are NAs for both adjusted p-values. Remove them
mRNA_expresion_deseq2_table <- mRNA_expresion_deseq2_table %>% 
  dplyr::filter(!(is.na(OSS_vs_LSS_paired_padj) & is.na(ESS_vs_LSS_paired_padj)))
