# IPA workflow

# STEP 1: CREATE THE TABLES THAT WILL BE USED IN IPA

# mRNA for both contrast

# FIRST TABLE WITH ALL THE GENES
# create the table by selected only specific columns
expression_database <- read_tsv(file = "~/data/Steve_White/results/Giulio/2026/expression_database.tsv") %>%
                        rename(Target_ID = geneID)
mRNA_table_for_IPA <- expression_database %>% 
  dplyr::select(Target_ID,
                OSS_vs_LSS_paired_log2FoldChange,
                OSS_vs_LSS_paired_padj,
                OSS_vs_LSS_paired_baseMean,
                ESS_vs_LSS_paired_log2FoldChange,
                ESS_vs_LSS_paired_padj,
                ESS_vs_LSS_paired_baseMean,
                OSS_vs_ESS_paired_log2FoldChange,
                OSS_vs_ESS_paired_padj,
                OSS_vs_ESS_paired_baseMean
                )


# ----------------------------------------

# SECOND TABLE
# create the table by selected only specific columns
# 
protein_table_for_IPA <- statistics_proteomics_data %>% 
  # Filtered out all the abundance ratios that are equal to 100 or 0.01 (they might be false positive)
  dplyr::mutate(protein_OSS_vs_LSS_adj_pvalue = ifelse(protein_OSS_vs_LSS_Abundance_Ratio == 100 | protein_OSS_vs_LSS_Abundance_Ratio == 0.01, NA, protein_OSS_vs_LSS_adj_pvalue),
                protein_OSS_vs_LSS_Abundance_Ratio = ifelse(protein_OSS_vs_LSS_Abundance_Ratio == 100 | protein_OSS_vs_LSS_Abundance_Ratio == 0.01, NA, protein_OSS_vs_LSS_Abundance_Ratio),
                protein_ESS_vs_LSS_adj_pvalue = ifelse(protein_ESS_vs_LSS_Abundance_Ratio == 100 | protein_ESS_vs_LSS_Abundance_Ratio == 0.01, NA, protein_ESS_vs_LSS_adj_pvalue),
                protein_ESS_vs_LSS_Abundance_Ratio = ifelse(protein_ESS_vs_LSS_Abundance_Ratio == 100 | protein_ESS_vs_LSS_Abundance_Ratio == 0.01, NA, protein_ESS_vs_LSS_Abundance_Ratio)) %>% 
   # Convert the abundance ratio of the proteins with a log2 transformation
   dplyr::mutate(protein_OSS_vs_LSS_Abundance_Ratio = ifelse(!is.na(protein_OSS_vs_LSS_Abundance_Ratio), log2(protein_OSS_vs_LSS_Abundance_Ratio), protein_OSS_vs_LSS_Abundance_Ratio),
                 protein_ESS_vs_LSS_Abundance_Ratio = ifelse(!is.na(protein_ESS_vs_LSS_Abundance_Ratio), log2(protein_ESS_vs_LSS_Abundance_Ratio), protein_ESS_vs_LSS_Abundance_Ratio)) %>% 
   # Select only column of interest
   dplyr::select(Protein_ID,
                 Target_ID,
                 protein_OSS_vs_LSS_Abundance_Ratio,
                 protein_OSS_vs_LSS_adj_pvalue,
                 protein_ESS_vs_LSS_Abundance_Ratio,
                 protein_ESS_vs_LSS_adj_pvalue
                 )

#Both tables will be used as input in IPA to run the pathway analysis

# STEP 2: FROM IPA TO THE VISUALIZATION ON RSTUDIO

#From now onward, the data used in this script are the ones obtained by IPA. 2 different analysis were run on IPA, one analysis per contrast (ESS vs LSS, OSS vs LSS). The whole list of expressed genes was used as input for each analysis. During each single analysis no filter on the cell types was used. Regarding the Expr. Fold Change, cutoff values of -1 and 1 were used while for the Expr FDR, a cut off value of 0.05 was used.

#Question that we can address with IPA:
#- which signaling and metabolic canonical pathways are enriched in the data?

# This can be shown in the canonical pathways tab. Here we can find a box plot (Also present in this script) that tell us which are the most significant pathways associated with our data and also, depending on the z-score, if those pathways are activated or inhibited.
  
# - Which predicted upstream regulators are significant?

# IPA considers any molecule that has a downstream effect to be an upstream regulator. So it considers Transcription factors, miRNAs, Drugs, compounds, etc.
# IPA is looking for upstream regulators that are able to change the expression, transcription, protein-DNA binding, phosphorylation of their targets.
# In this type of analysis, as for the canonical pathway, we are looking at the overlap between what is know in literature regarding targets of upstream regulators
# and my data.
  
# Each row in the table is an upstream regulator with the corresponding p-value or FDR (if we decide to include it). Some upstream regulators have a value
# for Expr. Fold Change while others do not. We have to keep in mind that an upstream regulator can be activated without inducing its expression. For example,
# there could be post-transnational modifications that induce the activation of that upstream regulator. Then this activated molecule has an affect on its targets.
# This takes into account the biological complexity.
  
# For each upstream regulator, we get the Expr. Fold change of the targets and if the information about their regulation coming from literature.
# Ex. A gene is regulated and the Fold Change is high, at the same time we could see that also in literature this gene is considered to be over expressed when regulated by that up-regulator.
  
# There is another option called Mechanistic Network and Causal Network that helps when we want to consider the effect of multiple regulators (Between the targets and the master upstream regulator there could be multiple intermediate regulators).
  
# - Which diseases and functions are over-represented?

# There is a heat map (the dimension of the boxes are given by the p-value of that disease. The bigger, the more significant that disease is. This can be changed with the z-score. Same thing also for the color). If we click on one of the boxes we will find in the end a list of genes with data about their expression coming from our data and the prediction coming from the literature (as in the output of the upstream regulators).
  
# - Regulator effects
  
# We can get a table where each row is a network (ranked by the consistency score). We get the information about the number of nodes present in the network divided by regulator nodes, target nodes and disease/function nodes. The last column is a rate between the relation between regulators and diseases/functions known in literature and the ones present in our data (Ex. 43% means that 43% of the relations between regulators and diseases/functions present in the network are also present in literature. All the others are novel). We can use the path tracer toll to highlight a specific path.
  
# - Networks
  
# To understand which molecules are highly interconnected and we can also see the diseases/functions that are related to these molecules/genes. Each network has a score and it depends on the number of focused molecules (molecules that passed the cutoffs and filters set during the core analysis. The higher the score, the more significant are the molecules that take part in that network). If we click on the molecule some of them are bold which are the ones present in our data or in general that are focus molecules. If we go on overlapping networks we can see if there are genes in common among networks.

#Statistical analysis in IPA
# 1) p-value (right-tailed Fisher's Exact test is used for this calculation)
#   - Null Hypothesis: Any association between the molecules in my dataset and the molecules in a particular disease/pathway is due to chance alone.
#   - Significant p-values under 0.05
#   - Possible correction: Benjamini-Hochberg correction for multiple testing
   
# 2) z-score
#   - Used for prediction about activation or inhibition (it compared the expression patterns of the genes in my data with what is known in literature databases)
#   - z-score >= 2 (activation) and z-score <= -2 (inhibition)
   
# The z-score can be used also for measure the match between two analysis:
#    - z-score >= 2 (match) and z-score <= -2 (anti-match)


## Import the data regarding the canonical pathways and the upstream regulators for the ESS vs LSS confrontation
### N.B. In the file names, the abbreviation "wctf" means "without cell type filter"

## First of all, define a function that allowed to modify the structure of the data to read the files in a correct way. Each file contains as first line "© 2000-2024 QIAGEN. All rights reserved." while the second line is empty and then the third contains the headers for the tabular data. For this reason the first two lines where removed.

# GRS
ipa_dir <- "" #Select the directory
dir.create(ipa_dir, recursive = TRUE)

# GRS - add check and upset plot
### checks

# check effect sizes add as expected
# GRS 
check <- mRNA_table_for_IPA %>% mutate(OSS_vs_LSS_minus_ESS_vs_LSS = 
                                         OSS_vs_LSS_paired_log2FoldChange - ESS_vs_LSS_paired_log2FoldChange)
ggplot(check, aes(x = OSS_vs_ESS_paired_log2FoldChange, y = OSS_vs_LSS_minus_ESS_vs_LSS)) + geom_point() + 
  ggtitle("sum of effect sizes") + theme_bw(base_size = 14)


## upset plot 
upset_data <- list()
for (ctr in c("OSS_vs_LSS", "ESS_vs_LSS", "OSS_vs_ESS")) {
   fckey = paste0(ctr, "_paired_log2FoldChange")
   pkey = paste0(ctr, "_paired_padj")  
   bmkey = paste0(ctr, "_paired_baseMean")  # could add & !!sym(bmkey) > 10
   upset_data[[paste0(ctr, "_up")]] <- filter(check, !!sym(pkey) < 0.05 & !!sym(fckey) > 1 ) %>% 
     pull(Target_ID)
   upset_data[[paste0(ctr, "_down")]] <- filter(check, !!sym(pkey) < 0.05 & !!sym(fckey) < -1 ) %>% pull(Target_ID)
}
upset_data4 <- upset_data[1:4]

upset_fn <- function(upset_data) {
  require(UpSetR)
  s <- 1.2
  upset(fromList(upset_data), nintersects = 32, order.by = "degree", 
        mb.ratio = c(0.55, 0.45), empty.intersections = NULL,
        sets = rev(names(upset_data)), keep.order = TRUE ,
        point.size = s*2.5, line.size = s*1, 
        text.scale = s*c(2.0, 2.0, 2.0, 1.5, 1.8, 1.6))

}
upset_fn(upset_data4)
upset_fn(upset_data)
png(file.path(ipa_dir, "threeway_overlaps_mRNA_upset.png"), w=2400, h=1375, res=300) # 3000, 1500
upset_fn(upset_data)
dev.off()

# GRS - output OSS vs ESS

library(annotables)

ctr <-  "OSS_vs_ESS"
fckey = paste0(ctr, "_paired_log2FoldChange")
pkey = paste0(ctr, "_paired_padj")  
bmkey = paste0(ctr, "_paired_baseMean")  # could add & !!sym(bmkey) > 10


check2 <- left_join(check, grch38, by = c("Target_ID" = "ensgene"))
outdat <- filter(check2, !!sym(pkey) < 0.05 & (!!sym(fckey) > 1 | !!sym(fckey) < -1))  %>% 
  select(Target_ID, symbol, all_of(c(fckey, pkey))) %>% 
  rename("gene_name" = "symbol") %>%
  arrange(desc(!!sym(fckey)))

write_csv(outdat, file = file.path(ipa_dir, "DE_genes_OSS_vs_ESS.csv"))

remove_first_two_lines <- function(file_name) {
  
  # Define the complete file path
  complete_file_path <- paste0(ipa_dir, file_name)
    
  # Remove the first two lines from the file and save the new file without the first two lines in the "data_content" variable
  data_content <- readLines(complete_file_path)[-c(1,2)]
  
  # Overwrite the already existing file, with the one where the first two lines were removed
  writeLines(data_content, paste0(ipa_dir, "/modified_", file_name))
  
  message(paste0(file_name, " modified!"))
  
}

## Modify the data by removing the first two lines

# Define a vector with all the names of the files that need to have the first two lines removed
# file_name_matrix <- matrix(c("mRNA_ESS_vs_LSS_expression_analysis.tab=pathways.bar_plot_table.txt",
#                              "mRNA_OSS_vs_LSS_expression_analysis.tab=pathways.bar_plot_table.txt",
#                              "mRNA_ESS_vs_LSS_expression_analysis.tab=Upstream_Analaysis.Upstream_regulators_table.txt",
#                              "mRNA_OSS_vs_LSS_expression_analysis.tab=Upstream_Analaysis.Upstream_regulators_table.txt",
#                              "protein_ESS_vs_LSS_expression_analysis.tab=pathways.bar_plot_table.txt",
#                              "protein_OSS_vs_LSS_expression_analysis.tab=pathways.bar_plot_table.txt",
#                              "protein_ESS_vs_LSS_expression_analysis.tab=Upstream_Analaysis.Upstream_regulators_table.txt",
#                              "protein_OSS_vs_LSS_expression_analysis.tab=Upstream_Analaysis.Upstream_regulators_table.txt",
#                              "comparison_analysis.tab=pathways.heat_map_table.txt",
#                              "comparison_analysis.tab=Upstream_Anslysis.Upstream_regulators_heat_map_table.txt",
#                              "network_table.txt"
#                              ),
#                            ncol = 1)

file_name_matrix <- matrix(list.files(ipa_dir))

# Remove the first two lines from all the files
apply(file_name_matrix, 1, remove_first_two_lines)


## Import the data regarding the canonical pathways and the upstream regulators for the OSS vs LSS confrontation

# mRNA

# Cananical pathways ("cp") 
mRNA_OSS_vs_LSS_contrast_wctf_cp_path <- file.path(ipa_dir, "modified_mRNA_OSS_vs_LSS_bar_plot_table.txt")
mRNA_OSS_vs_LSS_contrast_wctf_cp <- read.delim(mRNA_OSS_vs_LSS_contrast_wctf_cp_path, header = TRUE)

# Modify the dataframe by adding the name of the contrast to the column X and then change its name to "contrast_name"
mRNA_OSS_vs_LSS_contrast_wctf_cp <- mRNA_OSS_vs_LSS_contrast_wctf_cp %>% 
  dplyr::mutate(X = "OSS_vs_LSS (mRNAs)") %>% 
  dplyr::rename(contrast_name = X)

# Upstream regulators ("ur")
mRNA_OSS_vs_LSS_contrast_wctf_ur_path <- file.path(ipa_dir, "modified_mRNA_OSS_vs_LSS_Upstream_regulators_table.txt")
mRNA_OSS_vs_LSS_contrast_wctf_ur <- read.delim(mRNA_OSS_vs_LSS_contrast_wctf_ur_path, header = TRUE)

# ----------------------------------------

# Proteomics

# # Cananical pathways ("cp") 
# protein_OSS_vs_LSS_contrast_wctf_cp_path <- "~/Project_Giulio_Vidotto/IPA_output/output/modified_protein_OSS_vs_LSS_expression_analysis.tab=pathways.bar_plot_table.txt"
# protein_OSS_vs_LSS_contrast_wctf_cp <- read.delim(protein_OSS_vs_LSS_contrast_wctf_cp_path, header = TRUE)
# 
# # Modify the dataframe by adding the name of the contrast to the column X and then change its name to "contrast_name"
# protein_OSS_vs_LSS_contrast_wctf_cp <- protein_OSS_vs_LSS_contrast_wctf_cp %>% 
#   dplyr::mutate(X = "OSS_vs_LSS (proteins)") %>% 
#   dplyr::rename(contrast_name = X)
# 
# # Upstream regulators ("ur")
# protein_OSS_vs_LSS_contrast_wctf_ur_path <- "~/Project_Giulio_Vidotto/IPA_output/output/modified_protein_OSS_vs_LSS_expression_analysis.tab=Upstream_Analaysis.Upstream_regulators_table.txt"
# protein_OSS_vs_LSS_contrast_wctf_ur <- read.delim(protein_OSS_vs_LSS_contrast_wctf_ur_path, header = TRUE)


## Import the data regarding the canonical pathways and the upstream regulators for the ESS vs LSS confrontation


# mRNA

# Cananical pathways ("cp") 
mRNA_ESS_vs_LSS_contrast_wctf_cp_path <- file.path(ipa_dir,"modified_mRNA_ESS_vs_LSS_bar_plot_table.txt")
mRNA_ESS_vs_LSS_contrast_wctf_cp <- read.delim(mRNA_ESS_vs_LSS_contrast_wctf_cp_path, header = TRUE)

# Modify the dataframe by adding the name of the contrast to the column X and then change its name to "contrast_name"
mRNA_ESS_vs_LSS_contrast_wctf_cp <- mRNA_ESS_vs_LSS_contrast_wctf_cp %>% 
  dplyr::mutate(X = "ESS vs LSS (mRNAs)") %>% 
  dplyr::rename(contrast_name = X)

# Upstream regulators ("ur")
mRNA_ESS_vs_LSS_contrast_wctf_ur_path <- file.path(ipa_dir, "modified_mRNA_ESS_vs_LSS_Upstream_regulators_table.txt")
mRNA_ESS_vs_LSS_contrast_wctf_ur <- read.delim(mRNA_ESS_vs_LSS_contrast_wctf_ur_path, header = TRUE, na = "")

# ----------------------------------------

# Import the data for canonical pathways for Proteomics
# GRS addition

# # Cananical pathways ("cp") 
# protein_ESS_vs_LSS_contrast_wctf_cp_path <- "~/Project_Giulio_Vidotto/IPA_output/output/modified_protein_ESS_vs_LSS_expression_analysis.tab=pathways.bar_plot_table.txt"
# protein_ESS_vs_LSS_contrast_wctf_cp <- read.delim(protein_ESS_vs_LSS_contrast_wctf_cp_path, header = TRUE)
# 
# # Modify the dataframe by adding the name of the contrast to the column X and then change its name to "contrast_name"
# protein_ESS_vs_LSS_contrast_wctf_cp <- protein_ESS_vs_LSS_contrast_wctf_cp %>% 
#   dplyr::mutate(X = "ESS vs LSS (proteins)") %>% 
#   dplyr::rename(contrast_name = X)
# 
# # Upstream regulators ("ur")
# protein_ESS_vs_LSS_contrast_wctf_ur_path <- ""
# protein_ESS_vs_LSS_contrast_wctf_ur <- read.delim(protein_ESS_vs_LSS_contrast_wctf_ur_path, header = TRUE, na = "")
protein_ipa_path <- ""
prot_ESS_vs_LSS_contrast_wctf_cp <- read.delim(file.path(protein_ipa_path, "prot_ESS_vs_LSS_bar_plot_table.txt"), skip = 2, header = TRUE, na = "") %>%
    dplyr::mutate(X = "ESS vs LSS (proteins)") %>% 
    dplyr::rename(contrast_name = X) %>% 
    dplyr::rename(X.log.p.value. = X.log.B.H.p.value.)
prot_OSS_vs_LSS_contrast_wctf_cp <- read.delim(file.path(protein_ipa_path, "prot_OSS_vs_LSS_bar_plot_table.txt"), skip = 2, header = TRUE, na = "") %>%
    dplyr::mutate(X = "OSS vs LSS (proteins)") %>% 
    dplyr::rename(contrast_name = X) %>% 
    dplyr::rename(X.log.p.value. = X.log.B.H.p.value.)
prot_OSS_vs_ESS_contrast_wctf_cp <- read.delim(file.path(protein_ipa_path, "prot_OSS_vs_ESS_bar_plot_table.txt"), skip = 2, header = TRUE, na = "") %>%
    dplyr::mutate(X = "OSS vs ESS (proteins)") %>% 
    dplyr::rename(contrast_name = X) %>% 
    dplyr::rename(X.log.p.value. = X.log.B.H.p.value.)

```


## Import the data regarding the canonical pathways and the upstream regulators for the OSS vs ESS confrontation

#GRS addition

# mRNA

# Cananical pathways ("cp") 
mRNA_OSS_vs_ESS_contrast_wctf_cp_path <- file.path(ipa_dir,"modified_mRNA_OSS_vs_ESS_bar_plot_table.txt")
mRNA_OSS_vs_ESS_contrast_wctf_cp <- read.delim(mRNA_OSS_vs_ESS_contrast_wctf_cp_path, header = TRUE)

# Modify the dataframe by adding the name of the contrast to the column X and then change its name to "contrast_name"
mRNA_OSS_vs_ESS_contrast_wctf_cp <- mRNA_OSS_vs_ESS_contrast_wctf_cp %>% 
  dplyr::mutate(X = "OSS vs ESS (mRNAs)") %>% 
  dplyr::rename(contrast_name = X)

# Upstream regulators ("ur")
mRNA_OSS_vs_ESS_contrast_wctf_ur_path <- file.path(ipa_dir, "modified_mRNA_OSS_vs_ESS_Upstream_regulators_table.txt")
mRNA_OSS_vs_ESS_contrast_wctf_ur <- read.delim(mRNA_OSS_vs_ESS_contrast_wctf_ur_path, header = TRUE, na = "")

# ----------------------------------------


# For each single canonical pathway, a bar plot will be created. The bar plot shows which are the pathways that are significant over a specific threshold for the adjusted p-values

#GRS: added mRNA_OSS_vs_ESS to the below 

source("../function_scripts/IPA_barplot_function.Rmd")

# Create a list with all the dataframes for each confrontation
# IPA_list <- list(mRNA_ESS_vs_LSS = mRNA_ESS_vs_LSS_contrast_wctf_cp,
#                  mRNA_OSS_vs_LSS = mRNA_OSS_vs_LSS_contrast_wctf_cp,
#                  protein_ESS_vs_LSS = protein_ESS_vs_LSS_contrast_wctf_cp,
#                  protein_OSS_vs_LSS = protein_OSS_vs_LSS_contrast_wctf_cp)

IPA_list <- list(mRNA_ESS_vs_LSS = mRNA_ESS_vs_LSS_contrast_wctf_cp,
                 mRNA_OSS_vs_LSS = mRNA_OSS_vs_LSS_contrast_wctf_cp,
                 mRNA_OSS_vs_ESS = mRNA_OSS_vs_ESS_contrast_wctf_cp)

# Plot the graphs for each one of the dataframes (they will show the 10 most significant pathways)
# In these graph, the -log(q-values) are plotted. This means that the longer the bar, the higher the significance is (this means that the q-value is really small)
canonical_pathways_plot_list <- lapply(IPA_list, function(IPA_file) {IPA_barplot(file = IPA_file,
                                                                                 threshold = 0.05,
                                                                                 n = 30)}
  )

i=3
contr <- names(IPA_list)[i]
  show(canonical_pathways_plot_list[i])
  png(file.path(ipa_dir, paste0(contr, "_30_IPA_CP.png")), w=1800, h=1800, res=300)
  show(canonical_pathways_plot_list[i])
  dev.off()
                                      
for (i in 1:3) {
  contr <- names(IPA_list)[i]
  show(canonical_pathways_plot_list[i])
  png(file.path(ipa_dir, paste0(contr, "_IPA_CP.png")), w=1800, h=1000, res=300)
  show(canonical_pathways_plot_list[i])
  dev.off()
}
# IPA_barcode_for_specific_gene_set

# OSS vs LSS
# Only for test group 1
# GRS IPA_barcode_for_specific_gene_set(protein_OSS_vs_LSS_contrast_wctf_cp, "test group 1", threshold = 0.05, n = 30)

# Only for test group 2 and 3
IPA_barcode_for_specific_gene_set(mRNA_OSS_vs_LSS_contrast_wctf_cp, "test group 2", threshold = 0.05, n = 30)

# ESS vs LSS
# Only for test group 1
# GRS IPA_barcode_for_specific_gene_set(protein_ESS_vs_LSS_contrast_wctf_cp, "test group 1", threshold = 0.05, n = 30)

# Only for test group 2 and 3
IPA_barcode_for_specific_gene_set(mRNA_ESS_vs_LSS_contrast_wctf_cp, "test group 2", threshold = 0.05, n = 30)


# GRS
IPA_barcode_for_specific_gene_set(mRNA_OSS_vs_ESS_contrast_wctf_cp, "test group 2", threshold = 0.05, n = 30)

# GRS: bar plots for proteomics 

# GRS: this does not work; open file and click the code chunk
# source("../function_scripts/IPA_barplot_function.Rmd")

# Create a list with all the dataframes for each confrontation
# IPA_list <- list(mRNA_ESS_vs_LSS = mRNA_ESS_vs_LSS_contrast_wctf_cp,
#                  mRNA_OSS_vs_LSS = mRNA_OSS_vs_LSS_contrast_wctf_cp,
#                  protein_ESS_vs_LSS = protein_ESS_vs_LSS_contrast_wctf_cp,
#                  protein_OSS_vs_LSS = protein_OSS_vs_LSS_contrast_wctf_cp)

IPA_list <- list(prot_ESS_vs_LSS = prot_ESS_vs_LSS_contrast_wctf_cp,
                 prot_OSS_vs_LSS = prot_OSS_vs_LSS_contrast_wctf_cp,
                 prot_OSS_vs_ESS = prot_OSS_vs_ESS_contrast_wctf_cp)

# Plot the graphs for each one of the dataframes (they will show the 10 most significant pathways)
# In these graph, the -log(q-values) are plotted. This means that the longer the bar, the higher the significance is (this means that the q-value is really small)
canonical_pathways_plot_list <- lapply(IPA_list, function(IPA_file) {IPA_barplot(file = IPA_file,
                                                                                 threshold = 0.05,
                                                                                 n = 10)}
  )

i=3
contr <- names(IPA_list)[i]
show(canonical_pathways_plot_list[i])
  png(file.path(ipa_dir, paste0(contr, "_30_IPA_CP.png")), w=1800, h=1800, res=300)
  show(canonical_pathways_plot_list[i])
  dev.off()
                                      
for (i in 1:3) {
  contr <- names(IPA_list)[i]
  show(canonical_pathways_plot_list[i])
  png(file.path(ipa_dir, paste0(contr, "_IPA_CP.png")), w=1800, h=1000, res=300)
  show(canonical_pathways_plot_list[i])
  dev.off()
}
# IPA_barcode_for_specific_gene_set

# OSS vs LSS
# Only for test group 1
# GRS IPA_barcode_for_specific_gene_set(protein_OSS_vs_LSS_contrast_wctf_cp, "test group 1", threshold = 0.05, n = 30)

# Only for test group 2 and 3
IPA_barcode_for_specific_gene_set(mRNA_OSS_vs_LSS_contrast_wctf_cp, "test group 2", threshold = 0.05, n = 30)

# ESS vs LSS
# Only for test group 1
# GRS IPA_barcode_for_specific_gene_set(protein_ESS_vs_LSS_contrast_wctf_cp, "test group 1", threshold = 0.05, n = 30)

# Only for test group 2 and 3
IPA_barcode_for_specific_gene_set(mRNA_ESS_vs_LSS_contrast_wctf_cp, "test group 2", threshold = 0.05, n = 30)


# GRS
IPA_barcode_for_specific_gene_set(mRNA_OSS_vs_ESS_contrast_wctf_cp, "test group 2", threshold = 0.05, n = 30)

IPA_barcode_for_specific_gene_set <- function(IPA_file, group_name, threshold, n) {
  
  if (group_name %in% c("test group 1", "test group 2", "test group 3")) {
    
    # Load the intermediate data
    load("intermediate_data.RData")
    
    # Get the gene names from the set specified by the user
    gene_set <- intermediate_data[[gsub(" ", "_", group_name)]] %>% 
      dplyr::filter(targeted_by_miRNA == 1) %>% 
      dplyr::select(gene_name)
    print(gene_set)
    # Transform the gene set into a vector
    gene_set <- unlist(gene_set)
    
    # Get the pathways for the molecules only present in the gene set
    IPA_file <- IPA_file %>% 
      separate_rows(Molecules, sep = ",") %>% 
      dplyr::filter(Molecules %in% gene_set) %>% 
      dplyr::distinct(Ingenuity.Canonical.Pathways, .keep_all = TRUE) %>% 
      dplyr::select(Ingenuity.Canonical.Pathways, X.log.p.value., Ratio, z.score, contrast_name)
    print(IPA_file)
    IPA_barplot(file = IPA_file,
                threshold = 0.05,
                n = 30
                )
    
  } else {
    
    msg <- "The group name must be \"test group 1\", \"test group 2\" or \"test group 3\""
    stop(msg)
    
  }
  
}

# LOOK AT THE MEANING OF THE GENES PRESENT IN ONE OF THE TEST GROUPS (DEFINED IN THE STATISTICAL PIPELINE)


find_genes <- function(IPA_data,
                       canonical_pathways_file,
                       intermediate_data,
                       miR_mRNA_interaction,
                       altered_miR_database,
                       common_target_threshold,
                       altered_miRNAs_type,
                       adj_p_value_threshold,
                       fold_change_threshold,
                       test_num,
                       top_pathway_num) {
  
  # 1) EXTRACT THE GENE IN A SPECIFIC TEST GROUP
  
  test_group_data <- intermediate_data[[paste0("test_group_", test_num)]] %>% 
    dplyr::filter(targeted_by_miRNA == 1)
  print(unique(test_group_data$miRNA_ID))
  
  # 2) IMPORT THE IPA DATA AND TAKE ONLY THE GENES IN COMMON WITH THE TEST DATA
  
  IPA_data <- IPA_data %>% 
    dplyr::rename(upstream_regulator_gene = Upstream.Regulator,
                  downstream_regulated_gene = Target.Molecules.in.Dataset)
  
  View(IPA_data)
  
  # 3) SELECT ALL THE POSSIBILE MIRNAs-GENEs INTERACTIONS
  
  ### Extract the contrast from the name of the database of the miRNAs
  contrast <- gsub("_miR_sequencing_database", "", deparse(substitute(altered_miR_database)))
  
  ### From the contast, extract the two conditions
  condition_1 <- gsub("_vs_.*","", contrast)
  condition_2 <- gsub(".*_vs_","", contrast)
  
  ### Extract the common target of the miRNA in a specific condition (filter the miRNAs based on the ones present in the test gene data)
  common_target_db <- extract_miRNA_info(altered_miR_database = altered_miR_database,
                                         altered_miRNAs_type = altered_miRNAs_type,
                                         adj_p_value_threshold = adj_p_value_threshold,
                                         fold_change_threshold = fold_change_threshold,
                                         condition_1 = condition_1,
                                         condition_2 = condition_2,
                                         common_target_threshold = common_target_threshold) %>% 
    dplyr::rename(gene_name = common_target_name)
  
  # Check that the common target database is not empty
  if (nrow(common_target_db) == 0) {
    
    msg <- "Error: No miRNA targeted genes found"
    stop(msg)
    
  }
  
  
  # 3) CREATE TWO DIFFERENT TABLES FOR UPSTREAM AND DOWNSTREAM GENES
  
  ## 3.1) Downstream regulated genes
  
  ### 3.1.1) Get the downstream regulated genes
  downstream_regulated_gene_data <- IPA_data %>%
    separate_rows(downstream_regulated_gene, sep = ",") %>% 
    dplyr::distinct(gene_name = downstream_regulated_gene)
  
  ### 3.1.2) Add the information about the miRNAs targeting the downstream regulated genes 
  downstream_regulated_gene_data <- downstream_regulated_gene_data %>% 
    dplyr::filter(gene_name %in% test_group_data$gene_name) %>% # Only the downstream genes present in the test data
    left_join(common_target_db, by = "gene_name") %>% # Add the miRNAs information
    dplyr::filter(!is.na(miRNA_ID)) %>% # Remove possible lines with NA (take only the ones targeted by miRNAs)
    dplyr::select(gene_name, miRNA_ID)
  
  View(downstream_regulated_gene_data)
  
  ### 3.1.4) Check if there are downstream regulated genes
  if (nrow(downstream_regulated_gene_data) == 0) {
    
    msg <- "Error: No downstream regulated genes found that are targeted by miRNAs in the test group"
    stop(msg)
    
  }

  
  ## 3.2) Upstream regulator genes
  
  ## 3.2.1) Take only the upstream genes that are regulating the downstream genes present in the test data (the ones that are regulated)
  upstream_regulator_gene_data <- IPA_data %>% 
    separate_rows(downstream_regulated_gene, sep = ",") %>% 
    dplyr::filter(downstream_regulated_gene %in% downstream_regulated_gene_data$gene_name) 
  
  View(upstream_regulator_gene_data)
  ### 3.2.2) Get the upstream regulator genes
  upstream_regulator_gene_data <- upstream_regulator_gene_data %>% 
    dplyr::distinct(gene_name = upstream_regulator_gene)
    
  ### 3.2.3) Add the information about the miRNAs
  upstream_regulator_gene_data <- upstream_regulator_gene_data %>% 
    left_join(common_target_db, by = "gene_name") %>% # Add the miRNAs information (keep the ones tha are targeted and the ones that are not)
    dplyr::select(gene_name, miRNA_ID)
  
  
  ### 3.2.4) Check if there are upstream regulator genes
  if (nrow(upstream_regulator_gene_data) == 0) {
    
    msg <- "Error: No upstream regulator genes found that regulate the downstream genes"
    stop(msg)
    
  }
  
  # 4) CREATE THE NETWORK WITH NODES AND EDGES
  
  # 4.1) Edges
  # Edges between upstream regulators and downstream regulator
  gene_gene_edges <- IPA_data %>% 
    separate_rows(downstream_regulated_gene, sep = ",") %>%
    dplyr::select(to = downstream_regulated_gene,
                  from = upstream_regulator_gene)

  # Edges between miRNAs and upstream regulator genes (only the ones present in the test group)
  miRNAs_up_st_gene_edges <- upstream_regulator_gene_data %>%
    dplyr::filter(!is.na(miRNA_ID)) %>% 
    dplyr::rename(to = gene_name, 
                  from = miRNA_ID)

  # Edges between miRNAs and downstream regulated genes
  miRNAs_down_st_gene_edges <- downstream_regulated_gene_data %>% 
    dplyr::rename(to = gene_name, 
                  from = miRNA_ID)

  # All the edges
  edges <- rbind(gene_gene_edges, miRNAs_up_st_gene_edges, miRNAs_down_st_gene_edges)
  
  
  # 4.2) Nodes
  # miRNA nodes 
  miRNA_nodes <- data.frame(label = c(miRNAs_up_st_gene_edges$from, miRNAs_down_st_gene_edges$from),
                            shape = "dot",
                            color = "red",
                            level = 1,
                            size = 30,
                            stringsAsFactors = FALSE)

  # miRNA targeted upstream regulator genes nodes
  miRNA_targeted_upstream_regulator_nodes <- data.frame(label = upstream_regulator_gene_data %>% 
                                                          dplyr::filter(!is.na(miRNA_ID)) %>% 
                                                          dplyr::select(label = gene_name),
                                                        shape = "square",
                                                        color = "green",
                                                        level = 3,
                                                        size = 15,
                                                        stringsAsFactors = FALSE)
  print(miRNA_targeted_upstream_regulator_nodes)
  # miRNA non targeted upstream regulator genes nodes
  miRNA_non_targeted_upstream_regulator_nodes <- data.frame(label = upstream_regulator_gene_data%>% 
                                                              dplyr::filter(is.na(miRNA_ID)) %>% 
                                                              dplyr::select(label = gene_name),
                                                            shape = "square",
                                                            color = "green",
                                                            level = 4,
                                                            size = 15,
                                                            stringsAsFactors = FALSE)
    
    
  # Downstream regulated genes
  downstream_target_nodes <- data.frame(label = downstream_regulated_gene_data$gene_name,
                                        shape = "triangle",
                                        color = "lightblue",
                                        level = 2,
                                        size = 20,
                                        stringsAsFactors = FALSE)
       
  # All the nodes 
  nodes <- rbind(miRNA_nodes, miRNA_targeted_upstream_regulator_nodes, miRNA_non_targeted_upstream_regulator_nodes, downstream_target_nodes) %>% 
    dplyr::mutate(id = label,
                  font.size = 15) %>% 
    dplyr::distinct(id, .keep_all = TRUE)
        
  # 4.3) Network
  legend_groups <- data.frame(
      label = c("miRNA", "upstream regulator gene", "downstream target gene"),
      shape = c("dot", "square", "triangle"),
      color = c("red", "green", "lightblue"),
      stringsAsFactors = FALSE
    )
  
  plot <- visNetwork(nodes, edges, height = "1000px", width = "150%") %>%
      visNodes() %>% 
      visEdges(arrows = "to", width = 0.1) %>%
      visHierarchicalLayout(direction = "LR", levelSeparation = 500) %>% 
      visLegend(
            useGroups = FALSE,
            addNodes = legend_groups,
            position = "right",
            ncol = 1,  # Change the number of columns if needed
            stepX = 70, # Horizontal step between items
            stepY = 100  # Vertical step between items
          ) %>% 
    visOptions(highlightNearest = list(enabled = TRUE,
                                       degree = 1,
                                       labelOnly = FALSE)
               )
  
  show(plot)
  
}

# Define the database of altered miRNAs in a specific contrast
# Possible altered_miR_database:
# 1) "OSS_vs_LSS_miR_sequencing_database"
# 2) "OSS_vs_ESS_miR_sequencing_database"
# 3) "ESS_vs_LSS_miR_sequencing_database"
# Write this directly inside the function

# Run the function to get the results (in this case we are interested in the intermediate databases about the groups)
# __ PARAMETERS __

config <- list(
  
  # __ PARAMETERS ABOUT THE FILTERING OF ALTERING MIRNAs __
  altered_miRNAs_type = "UP", # To analyse UP or DOWN regulated miRNAs (used to find the type of altered miRNAs)
  adj_p_value_threshold = 0.05, # Threshold to apply on the adjusted p-values of the miRNAs (used to find the altered miRNAs)
  fold_change_threshold = 1, # Threshold on the fold change (used to find the altered miRNAs)

  # __ PARAMETERS ABOUT THE FINDING OF GENES TARGETED BY ONE OR MORE MIRNA __
  common_target_threshold = 100, # Threshold on the number of common targets (used in the find_common_target function)
    
  # __ PARAMETERS ABOUT THE STATISTICAL TEST __
  miRNAs_threshold = 1, # Theshold to apply on the number of miRNAs targeting a gene (use in the get_graphs function)
  background_group_num = 1, # Background group to use for the statistical test
  about_test = "regulated genes test" # meaning of the test ("regulated genes test" or "KFERQ motif test")
  
)

# __ DATA __

# Data about mRNA and proteomics
mRNA_proteomics_data <- statistics_proteomics_data %>% 
  left_join(expression_database, by = "Target_ID")

# Define the database of altered miRNAs in a specific contrast
# Possible altered_miR_database:
# 1) "OSS_vs_LSS_miR_sequencing_database"
# 2) "OSS_vs_ESS_miR_sequencing_database"
# 3) "ESS_vs_LSS_miR_sequencing_database"
# Write this directly inside the function

# Run the function to get the results
pipeline_statistical_test(altered_miR_database = OSS_vs_LSS_miR_sequencing_database,
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

# Load the intermediate data
load("intermediate_group_data.RData")

find_genes(IPA_data = mRNA_OSS_vs_LSS_contrast_wctf_ur,
           canonical_pathways_file = mRNA_OSS_vs_LSS_contrast_wctf_cp,
           intermediate_data = intermediate_data,
           top_pathway_num = 100,
           test_num = 2,
           miR_mRNA_interaction = miR_mRNA_interaction,
           altered_miR_database = OSS_vs_LSS_miR_sequencing_database,
           common_target_threshold = config$common_target_threshold,
           altered_miRNAs_type = config$altered_miRNAs_type,
           adj_p_value_threshold = config$adj_p_value_threshold,
           fold_change_threshold = config$fold_change_threshold)


