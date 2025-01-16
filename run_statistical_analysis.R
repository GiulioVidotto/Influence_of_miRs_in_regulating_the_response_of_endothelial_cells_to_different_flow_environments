# Statistical analysis script

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

# Define the database of altered miRNAs in a specific contrast
# Possible altered_miR_database:
# 1) "OSS_vs_LSS_miR_sequencing_database"
# 2) "OSS_vs_ESS_miR_sequencing_database"
# Write this directly inside the function
# Run the function to get the results
pipeline_statistical_test(altered_miR_database = ESS_vs_LSS_miR_sequencing_database,
                          miR_mRNA_interaction = miR_mRNA_interaction,
                          mRNA_proteomics_data = mRNA_proteomics_data,
                          altered_miRNAs_type = config$altered_miRNAs_type,
                          adj_p_value_threshold = config$adj_p_value_threshold,
                          fold_change_threshold = config$fold_change_threshold,
                          common_target_threshold = config$common_target_threshold,
                          miRNAs_threshold = config$miRNAs_threshold,
                          background_group_num = config$background_group_num,
                          about_test = config$about_test,
                          save_intermediate = TRUE)