# Script Description:
#
# This script processes and enriches gene expression data from a DESeq2 pipeline.
# It integrates gene annotation and miRNA-target interaction data to produce
# a combined dataset containing both expression and interaction information.
#
# Steps:
# 1. Imports the DESeq2-processed expression data (`mRNA_expresion_deseq2_table`).
# 2. Merges this data with gene annotations from the original mRNA expression table.
#    - Cleans ENSEMBL IDs (removes version numbers like ".1", ".2").
#    - Adds gene symbols and genomic coordinates.
#    - Renames `geneID` to `Target_ID` and reorders columns.
# 3. Integrates the expression data with miRNA-target interaction information
#    from the `all_interactions` dataset.
# 4. Filters for genes that have known miRNA interactions.
# 5. Selects relevant columns for downstream analysis:
#    `miRNA_ID`, `Target_ID`, `gene_name`, `source_database`, and `number_of_database`.

# --- 1. Import the data (from DESeq2 pipeline) ---
expression_database <- mRNA_expresion_deseq2_table

# --- 2. Annotate with gene names and additional info from original expression table ---
expression_database <- expression_database %>% 
  left_join(mRNA_expression_database %>% 
              dplyr::mutate(geneID = gsub("\\.\\d*", "", geneID)) %>% 
              dplyr::select(geneID,
                            gene_name = genesymbol,
                            start,
                            stop),
            by = "geneID") %>% 
  dplyr::rename(Target_ID = geneID) %>% 
  dplyr::relocate(gene_name, .after = Target_ID)

# --- 3. Combine expression data with miRNA-target interaction data ---
miR_mRNA_interaction <- expression_database %>%
  dplyr::left_join(all_interactions, by = "Target_ID") %>% 
  dplyr::filter(!is.na(miRNA_ID)) %>% 
  dplyr::select(miRNA_ID, Target_ID, gene_name, source_database, number_of_database)