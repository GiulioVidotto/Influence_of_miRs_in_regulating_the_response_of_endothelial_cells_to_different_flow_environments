# __ PARAMETERS __
config <- list(
  # __ PARAMETERS ABOUT THE FILTERING OF ALTERING MIRNAs __
  altered_miRNAs_type = "DOWN", # ("UP", "DOWN") To analyse UP or DOWN regulated miRNAs (used to find the type of altered miRNAs)
  adj_p_value_threshold = 0.05, # Threshold to apply on the adjusted p-values of the miRNAs (used to find the altered miRNAs)
  fold_change_threshold = 1, # Threshold on the fold change (used to find the altered miRNAs)
  # __ PARAMETERS ABOUT THE FINDING OF GENES TARGETED BY ONE OR MORE MIRNA __
  common_target_threshold = 20, # Threshold on the number of common targets (used in the find_common_target function)
  # __ PARAMETERS ABOUT THE STATISTICAL TEST __
  miRNAs_threshold = 1, # Threshold to apply on the number of miRNAs targeting a gene (use in the get_graphs function)
  background_group_num = 1, # Background group to use for the statistical test
  about_test = "regulated genes test" # meaning of the test ("regulated genes test" or "KFERQ motif test")
)

# __ DATA __
# Data about mRNA and proteomics
mRNA_proteomics_data <- statistics_proteomics_data %>% 
  left_join(expression_database, by = "Target_ID") %>% 
  dplyr::filter(!is.na(gene_name))


# Calculate miRNA and mRNA concentrations
# 1) miRNA concentration
# Create the dataframe with the column data information
mIRNA_count_matrix <- miR_sequencing_database %>% 
  dplyr::mutate(miRNA_length = stop - start) %>% 
  dplyr::distinct(miRNA_ID, .keep_all = TRUE) %>% 
  dplyr::select(miRNA_ID, miRNA_length, contains("raw"))

# Calculate the length of all miRNAs
miRNA_legth <- as.vector(mIRNA_count_matrix$miRNA_length)
names(miRNA_legth) <- mIRNA_count_matrix$miRNA_ID

# Convert mIRNA_count_matrix to a data frame
mIRNA_count_matrix <- as.data.frame(mIRNA_count_matrix)

# Add miRNA IDs as rownames in the count matrix
rownames(mIRNA_count_matrix) <- mIRNA_count_matrix$miRNA_ID

# To obtain the counts run the DESeq2 pipeline
# Select only the raw count columns
mIRNA_count_matrix <- mIRNA_count_matrix %>% 
  dplyr::select(contains("raw"))

# Create column data frame for miRNA conditions and donors
coldata_miRNA <- data.frame(condition = factor(gsub("patient_[0-9]+_", "", colnames(mIRNA_count_matrix))),
                            donor = factor(gsub("patient_|_[A-Z]*", "", colnames(mIRNA_count_matrix)))
)

# Set the column names of coldata_miRNA as rownames
rownames(coldata_miRNA) <- colnames(mIRNA_count_matrix)

# Check that the order of the rows of the column data and of the columns of the count matrix is the same
all(rownames(coldata_miRNA) == colnames(mIRNA_count_matrix)) # TRUE

# Create DESeqDataSet object from count data and column data (metadata)
dds_miRNA <- DESeqDataSetFromMatrix(countData = mIRNA_count_matrix, 
                                    colData = coldata_miRNA, 
                                    design = ~ condition)

# Add miRNA length as an additional column in the DESeqDataSet
mcols(dds_miRNA)$basepairs <- miRNA_legth

# Perform differential expression analysis using DESeq
dds_miRNA <- DESeq(dds_miRNA)

# Calculate the RPKM values for miRNA using the DESeq2 results
miRNA_rpkm_values <- as.data.frame(fpkm(dds_miRNA))

# Add miRNA IDs as a column to the RPKM values dataframe
miRNA_rpkm_values <- miRNA_rpkm_values %>% 
  dplyr::mutate(miRNA_ID = rownames(miRNA_rpkm_values))

# Join miRNA RPKM values with fold change and p-value data for OSS vs LSS comparison
OSS_vs_LSS_miRNA_rpkm_values <- miRNA_rpkm_values %>% 
  left_join(miR_sequencing_database %>% 
              dplyr::distinct(miRNA_ID, .keep_all = TRUE) %>% 
              dplyr::select(miRNA_ID, OSS_vs_LSS_paired_log2FoldChange, OSS_vs_LSS_paired_padj),
            by = "miRNA_ID") %>% 
  dplyr::filter(miRNA_ID %in% OSS_vs_LSS_miR_sequencing_database$miRNA_ID)

# Join miRNA RPKM values with fold change and p-value data for ESS vs LSS comparison
ESS_vs_LSS_miRNA_rpkm_values <- miRNA_rpkm_values %>% 
  left_join(miR_sequencing_database %>% 
              dplyr::distinct(miRNA_ID, .keep_all = TRUE) %>% 
              dplyr::select(miRNA_ID, ESS_vs_LSS_paired_log2FoldChange, ESS_vs_LSS_paired_padj),
            by = "miRNA_ID") %>% 
  dplyr::filter(miRNA_ID %in% ESS_vs_LSS_miR_sequencing_database$miRNA_ID)



# 2) mRNA concentration
# Calculate the mRNA length for each gene
mRNA_length <- as.vector(expression_database$stop - expression_database$start)
names(mRNA_length) <- expression_database$Target_ID

# Add mRNA length as base pairs column to the DESeq2 dataset for mRNA data
mcols(paired_dds)$basepairs <- mRNA_length

# Calculate RPKM values for mRNA using the DESeq2 results
mRNA_rpkm_values <- as.data.frame(fpkm(paired_dds)) 

# Add Target_IDs as a column to the mRNA RPKM values dataframe
mRNA_rpkm_values <- mRNA_rpkm_values%>% 
  dplyr::mutate(Target_ID = rownames(mRNA_rpkm_values))

# Define the database of altered miRNAs in a specific contrast
# Possible altered_miR_database:
# 1) "OSS_vs_LSS_miR_sequencing_database"
# 3) "ESS_vs_LSS_miR_sequencing_database"
# Write this directly inside the function
# Create a list of classifier data by applying the function for both UP and DOWN miRNA alterations
classifier_data_list <- lapply(c("UP", "DOWN"), function(altered_miRNAs_type) {
  # Run the function to get the results
  pipeline_statistical_test(altered_miR_database = NULL
                            miR_mRNA_interaction = NULL,
                            mRNA_proteomics_data = NULL,
                            altered_miRNAs_type = altered_miRNAs_type,
                            adj_p_value_threshold = config$adj_p_value_threshold,
                            fold_change_threshold = config$fold_change_threshold,
                            common_target_threshold = config$common_target_threshold,
                            miRNAs_threshold = config$miRNAs_threshold,
                            background_group_num = config$background_group_num,
                            about_test = NULL,
                            save_intermediate = TRUE)
  # Load the intermediate data that was saved in the previous step
  load('intermediate_data.RData')
  # Depending on the contrast from the intermediate data (OSS_vs_LSS or ESS_vs_LSS), get the corresponding classifier data
  if (intermediate_data$contrast == "OSS_vs_LSS") {
    # Get the classifier data for the OSS vs LSS comparison
    classifier_data <- get_classifier_data(intermediate_data,
                                           OSS_vs_LSS_miR_sequencing_database,
                                           expression_database,
                                           OSS_vs_LSS_MFE_info)
    
  } else if (intermediate_data$contrast == "ESS_vs_LSS") {
    # Get the classifier data for the ESS vs LSS comparison
    classifier_data <- get_classifier_data(intermediate_data,
                                           ESS_vs_LSS_miR_sequencing_database,
                                           expression_database,
                                           ESS_vs_LSS_MFE_info)
    
  }
  # Add the miRNA alteration type (UP or DOWN) to the classifier data for reference
  classifier_data <- classifier_data %>% 
    dplyr::mutate(miRNA_type = altered_miRNAs_type)
  
  return(classifier_data)
  
})

# Load the intermediate data directly from the statistical test script

# The intermediate data is a list of multiple elements:
# - contrast -> it could be "ESS vs LSS" or "OSS vs LSS" (it depends on the contrast used in the statistical tests)
# - merged_data -> table with all he information with all the genes (targeted and non targeted)
# - background_group_1 <- table with genes where there is no change in proteomics and also no change in mRNA expression
# - background_group_2 <- table with genes where genes where the change in proteomics scales with change in mRNA expression up or down
# - test_group_1 <- table with genes where there are genes that don't change at mRNA, but do change in proteomics
# - test_group_2 <- table with genes where there are genes that change in proteomics, but don't change in mRNA
# - test_group_3 <- table with genes where there are genes where mRNA changes in one direction, but proteins change in the opposite

# Load the intermediate data that was saved earlier
load('intermediate_data.RData')

# Combine the list of classifier data into a single data frame
classifier_data_df <- do.call(rbind, classifier_data_list)

# Calculate the net effect score for each target gene
net_effect_score <- classifier_data_df %>%
  dplyr::distinct(Target_ID, miRNA_ID, .keep_all = TRUE) %>% 
  mutate(sign = ifelse(miRNA_type == "UP", 1, -1)) %>% 
  group_by(Target_ID) %>%
  summarise(net_score = sum(sign))

# Merge the calculated net effect score back into the original classifier data frame
classifier_data_df <- classifier_data_df %>% 
  left_join(net_effect_score, by = "Target_ID") %>% 
  dplyr::select(Target_ID, miRNA_ID, y, y_1, y_2, y_3, net_score, miRNA, mRNA, Kd)

# Write the final classifier data frame to a CSV file
write.csv(classifier_data_df, paste0(getwd(), "/binary_classification_scripts/", intermediate_data$contrast, "_classifier_data_final.csv"))
