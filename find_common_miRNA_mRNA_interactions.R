#The databases used for comparison are:
# - miRTarBase
# - miRDB
# - TarBase
# - Targetscan

## Uplaod the different databases
# Get the working directory path
current_working_directory <- getwd()
# MirTarBase
path_output_miRTarBase_database <- paste0(current_working_directory, "/project_data/miR_databases/final_outputs/output_miRTarBase_database.csv")
output_miRTarBase_database <- read.csv(path_output_miRTarBase_database, stringsAsFactors = FALSE)
# MiRDB
path_output_miRDB_database <- paste0(current_working_directory, "/project_data/miR_databases/final_outputs/output_miRDB_database.csv")
output_miRDB_database <- read.csv(path_output_miRDB_database, stringsAsFactors = FALSE)
# TarBase
path_output_tarBase_database <- paste0(current_working_directory, "/project_data/miR_databases/final_outputs/output_tarBase_database.csv")
output_tarBase_database <- read.csv(path_output_tarBase_database, stringsAsFactors = FALSE)
# TargetScan
path_output_targetScan_database <- paste0(current_working_directory, "/project_data/miR_databases/final_outputs/output_targetScan_database.csv")
output_targetScan_database <- read.csv(path_output_targetScan_database, stringsAsFactors = FALSE)

## Quick analysis on the databases
# Search how many miRNA present in Steve's data are also present in the different databases
n_miRNA_database <- function(database, miR_database) {
  if (!"checked_miRNA_ID" %in% colnames(database)) {
    msg <- "Error: The database must have the checked_miRNA_ID coulmn. this column is found after the check on the IDs"
    stop(msg)
  }
  column <- database$checked_miRNA_ID
  unique_miRNA <- miR_database %>% 
    dplyr::filter(miRNA_ID %in% column) %>% 
    dplyr::distinct(miRNA_ID)
  return(unique_miRNA)
}

# Create the list of databases and set the names for each one of them
list_of_coulmns <- list(output_miRTarBase_database, output_miRDB_database, output_tarBase_database, output_targetScan_database)
names(list_of_coulmns) <- c("miRTarBase", "miRDB", "tarBase", "targetScan")

# Apply the fucntion to each element of the list
# - OSS vs LSS
OSS_vs_LSS_n_miRNA_in_each_database <- lapply(list_of_coulmns, function(database) {
  n_miRNA_database(database = database,
                   miR_database = OSS_vs_LSS_miR_sequencing_database)
})

print(OSS_vs_LSS_n_miRNA_in_each_database)

# - ESS vs LSS
ESS_vs_LSS_n_miRNA_in_each_database <- lapply(list_of_coulmns, function(database) {
  n_miRNA_database(database = database,
                   miR_database = ESS_vs_LSS_miR_sequencing_database)
})

print(ESS_vs_LSS_n_miRNA_in_each_database)

# Check the miRNAs
# Function to find elements unique to one dataset
find_unique_elements <- function(data_list) {
  # Get the names of datasets
  names_list <- names(data_list)
  # Initialize a list to store results
  unique_elements <- list()
  # Loop through each dataset
  for (name in names_list) {
    # Current dataset
    current_data <- data_list[[name]]$miRNA_ID
    # Combine all datasets except the current one
    other_data <- unlist(lapply(data_list[names_list != name], function(df) df$miRNA_ID))
    # Find elements in the current dataset not in the other datasets
    unique_in_current <- setdiff(current_data, other_data)
    if (length(unique_in_current) > 0) {
      unique_elements[[name]] <- unique_in_current
    }
  }
  return(unique_elements)
}

# Find unique elements
unique_elements_OSS_LSS <- find_unique_elements(OSS_vs_LSS_n_miRNA_in_each_database)
unique_elements_ESS_LSS <- find_unique_elements(ESS_vs_LSS_n_miRNA_in_each_database)

## Create the interaction columns for all the databases (these databases are ckecked and they contain the info if there is a mismatch)
# The only one different is targetscan for which we were able to update some of the miRNAs so it is specified after

# creating a new column interactions for the miRtarBase database
interaction_miRTarBase <- output_miRTarBase_database %>% 
  mutate(interactions = paste(checked_miRNA_ID, Target_ID, sep = "_"))
# creating a new column interactions for the miRDB database
interaction_miRDB <- output_miRDB_database %>% 
  mutate(interactions = paste(checked_miRNA_ID, Target_ID, sep = "_"))

# creating a new column interactions for the tarBase database
interaction_tarBase <- output_tarBase_database %>% 
  mutate(interactions = paste(checked_miRNA_ID, Target_ID, sep ="_"))

# Find the interasections among the 4 databases
interaction_set <- list(miRTarBase = as.vector(interaction_miRTarBase %>% dplyr::distinct(checked_miRNA_ID)), 
                        miRDB = as.vector(interaction_miRDB %>% dplyr::distinct(checked_miRNA_ID)),
                        tarBase = as.vector(interaction_tarBase %>% dplyr::distinct(checked_miRNA_ID)),
                        targetScan = as.vector(targetScan_database_checked %>% dplyr::distinct(checked_miRNA_ID))
)
list_of_databases <- names(interaction_set)
databse_statistics_list <- lapply(list_of_databases, function(database_name) {
  database_checked_ID <- unlist(interaction_set[[database_name]])
  resul_df <- data.frame(database_name = database_name,
                         `5p_missing_in_miRBase` = sum(endsWith(database_checked_ID, ":'-5p' missing in miRBase")),
                         `3p_missing_in_miRBase` = sum(endsWith(database_checked_ID, ":'-3p' missing in miRBase")),
                         `3p_or_5p_specification_needed` = sum(endsWith(database_checked_ID, ":'-3p' or '-5p' specification needed")),
                         `5p_specification_needed` = sum(endsWith(database_checked_ID, ":'-5p' specification needed")),
                         `3p_specification_needed` = sum(endsWith(database_checked_ID, ":'-3p' specification needed")),
                         `Not_found` = sum(endsWith(database_checked_ID, ": Not found")),
                         stringsAsFactors = FALSE)
  return(resul_df)
})

databse_statistics_df <- do.call(rbind, databse_statistics_list)
colnames(databse_statistics_df) <- c("database name", "-5p missing in miRBase", "-3p missing in miRBase", "-3p or -5p specification needed", "-5p specification needed", "-3p specification needed", "Not found in miRBase")
databse_statistics_long <- pivot_longer(databse_statistics_df, 
                                        cols = -`database name`, 
                                        names_to = "category", 
                                        values_to = "count")

# Remove the categories were there are zeros
databse_statistics_long <- databse_statistics_long %>% 
  dplyr::filter(count > 0)

# Create the bar plot
ggplot(data = databse_statistics_long, aes(x = category, y = count, fill = category)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ `database name`) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Category",
       y = "Counts",
       title = "Database Statistics by Category"
  ) +
  theme(
    axis.title = element_text(size = 30),   # Increase axis title size
    axis.text = element_text(size = 30),    # Increase axis text size
    legend.position = "right",  # Adjust legend position
    legend.title = element_text(size = 30), 
    legend.text = element_text(size = 30),   # Increase legend text size
    strip.text = element_text(size = 30)
  ) +
  scale_fill_manual(values = c("-5p missing in miRBase" = "lightblue",
                               "-3p missing in miRBase" = "blue",
                               "-3p or -5p specification needed" = "darkblue",
                               "-5p specification needed" = "darkgreen",
                               "-3p specification needed" = "green",
                               "Not found in miRBase" = "darkgray"))

# get the interaction from the public databases but this time only the ones that match
# creating a new column interactions for the miRtarBase database
interaction_miRTarBase_miRBase <- output_miRTarBase_database %>% 
  dplyr::filter(!":" %in% checked_miRNA_ID) %>% 
  mutate(interactions = paste(checked_miRNA_ID, Target_ID, sep = "_"))

# creating a new column interactions for the miRDB database
interaction_miRDB_miRBase <- output_miRDB_database %>% 
  dplyr::filter(!":" %in% checked_miRNA_ID) %>%
  mutate(interactions = paste(checked_miRNA_ID, Target_ID, sep = "_"))

# creating a new column interactions for the tarBase database
interaction_tarBase_miRBase <- output_tarBase_database %>% 
  dplyr::filter(!":" %in% checked_miRNA_ID) %>%
  mutate(interactions = paste(checked_miRNA_ID, Target_ID, sep ="_"))

# creating a new column interactions for the targetScan database
interaction_targetScan_miRBase <- output_targetScan_database %>% 
  dplyr::filter(!":" %in% checked_miRNA_ID) %>%
  mutate(interactions = paste(checked_miRNA_ID, Target_ID, sep ="_"))

# UpSet plot
# Fist create a list with all the interactions
interaction_set <- list(miRTarBase = as.vector(interaction_miRTarBase_miRBase$interactions), 
                        miRDB = as.vector(interaction_miRDB_miRBase$interactions),
                        tarBase = as.vector(interaction_tarBase_miRBase$interactions),
                        targetScan = as.vector(interaction_targetScan_miRBase$interactions)
)

# Visualize the different sets as bars inside a barplot
bar_plot_df <- data.frame(df_names = c("miRDB", "miRTarBase", "tarBase", "targetScan"),
                          dim = c(length(interaction_set$miRTarBase),
                                  length(interaction_set$miRDB),
                                  length(interaction_set$tarBase),
                                  length(interaction_set$targetScan)),
                          stringsAsFactors = FALSE)

ggplot(data = bar_plot_df, aes(x = df_names, y = dim)) +
  geom_bar(width = 0.8, stat = "identity", fill = c("darkblue", "darkgreen", "gold", "red")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 20)) +
  xlab("public repostiory") +
  ylab("Set Size") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title = element_text(size = 30),   # Increase axis title size
    axis.text = element_text(size = 30),    # Increase axis text size
    legend.position = "right",  # Adjust legend position
    legend.title = element_text(size = 30), 
    legend.text = element_text(size = 30),   # Increase legend text size
    strip.text = element_text(size = 30)
  )

# Give a name to the databases
names(interaction_set) <- c("miRTarBase", "miRDB", "tarBase", "targetScan")

# Set up the upset plot
s = 1.2
object <- mod_fromList(input = interaction_set)
upset(object,
      order.by = "freq",
      mb.ratio = c(0.55, 0.45),
      sets = c("miRTarBase", "miRDB", "tarBase", "targetScan"),
      keep.order = TRUE,
      point.size = s*3.5,
      line.size = s*1,
      sets.bar.color = c("darkblue", "darkgreen", "gold", "red"),
      text.scale = c(4, 4, 4, 4, 4, 4)
)

# Show the plot
show(upset_graph)


# Modify the database resulting from the upset plot by adding two new columns (row sum and source)
# Data from the upset plot
upset_plot_data <- upset_graph$New_data
# Get the names of the databases
database_vector <- colnames(upset_plot_data)
# Function to get the source
add_source <- function(row, database_vector) {
  # Transform the binary values of a row in logical values
  logical_row <- as.logical(row)
  # Select the databases
  selected_database <- database_vector[logical_row]
  # Write everything in one string
  result <- paste(selected_database, collapse = ",")
  return(result)
}

# Create and add the two new columns to the database
upset_plot_data_new_columns <- upset_plot_data %>% 
  mutate(source = apply(upset_plot_data, MARGIN = 1, add_source, database_vector = database_vector),
         row_sum = rowSums(upset_plot_data))

## Use the function "select_interaction" to get the different elements in different intersection
# All interactions (among TargetScan, TarBase, miRDB, MiRTarBase)
all_interaction <- select_interaction(upset_plot_data_new_columns, 1, "greater or equal")
# Save the filae with all the interactions
write.csv(all_interaction, file.path(paste0(current_working_directory, "/project_data/miR_databases/final_outputs/all_interaction.csv")), row.names = FALSE)

