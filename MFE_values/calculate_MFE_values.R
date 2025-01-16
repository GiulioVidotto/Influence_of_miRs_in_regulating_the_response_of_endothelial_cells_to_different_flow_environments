# Before runnning the script obtain the path of the working directory
current_working_directory <- getwd()

# ___ Up-regulated miRNA and their targets interactions ___

# Short 3'UTR (under 2000 bp) - these values are calculated by RNAup
up_short_3UTR_thermodynamics_file_path <- paste0(current_working_directory, "/project_data/genome_data/MFE_values/filtered_up_thermodynamics_scores_file_rna_up")
up_short_3UTR_thermodynamics_file <- readLines(up_short_3UTR_thermodynamics_file_path)

# Long 3'UTR (over 2000 bp) - these values are calculated by RNAplfold and RNAplex
up_long_3UTR_thermodynamics_file_path <- paste0(current_working_directory, "/project_data/genome_data/MFE_values/up_thermodynamics_scores_file")
up_long_3UTR_thermodynamics_file <- readLines(up_long_3UTR_thermodynamics_file_path)

# After having uploaded the file on R, extract the MFE values and the miRNA and gene name for each interaction

gene_ID <- NULL
miRNA_ID <- NULL

up_short_3UTR_MFE_info_list <- lapply(up_short_3UTR_thermodynamics_file, function(line) {get_MFE(file_line = line)})
up_long_3UTR_MFE_info_list <- lapply(up_long_3UTR_thermodynamics_file, function(line) {get_MFE(file_line = line)})

# Create two different dataframes
up_short_3UTR_MFE_info <- do.call(rbind, up_short_3UTR_MFE_info_list)
up_long_3UTR_MFE_info <- do.call(rbind, up_long_3UTR_MFE_info_list)

# Unite the dataframes in a single one
up_MFE_info <- rbind(up_short_3UTR_MFE_info, up_long_3UTR_MFE_info)

# Calculate the equilibrium constant at a temperature of 37 C
up_MFE_info <- up_MFE_info %>% 
  mutate(Kd = get_Kd(as.integer(MFE_value),37))



# ___ Down-regulated miRNA and their targets interactions ___

# Short 3'UTR (under 2000 bp) - these values are calculated by RNAup
down_short_3UTR_thermodynamics_file_path <- paste0(current_working_directory, "/project_data/genome_data/MFE_values/filtered_down_thermodynamics_scores_file_rna_up")
down_short_3UTR_thermodynamics_file <- readLines(down_short_3UTR_thermodynamics_file_path)

# Long 3'UTR (over 2000 bp) - these values are calculated by RNAplfold and RNAplex
down_long_3UTR_thermodynamics_file_path <- paste0(current_working_directory, "/project_data/genome_data/MFE_values/down_thermodynamics_scores_file")
down_long_3UTR_thermodynamics_file <- readLines(down_long_3UTR_thermodynamics_file_path)

# After having uploaded the file on R, extract the MFE values and the miRNA and gene name for each interaction

gene_ID <- NULL
miRNA_ID <- NULL

down_short_3UTR_MFE_info_list <- lapply(down_short_3UTR_thermodynamics_file, function(line) {get_MFE(file_line = line)})
down_long_3UTR_MFE_info_list <- lapply(down_long_3UTR_thermodynamics_file, function(line) {get_MFE(file_line = line)})

# Create two different dataframes
down_short_3UTR_MFE_info <- do.call(rbind, down_short_3UTR_MFE_info_list)
down_long_3UTR_MFE_info <- do.call(rbind, down_long_3UTR_MFE_info_list)

# Unite the dataframes in a single one
down_MFE_info <- rbind(down_short_3UTR_MFE_info, down_long_3UTR_MFE_info)

# Calculate the equilibrium constant at a temperature of 37 C
down_MFE_info <- down_MFE_info %>% 
  mutate(Kd = get_Kd(as.integer(MFE_value),37))
