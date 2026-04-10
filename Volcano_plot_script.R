library(ggrepel) 
library(ggplot2)
library(annotables)
# genes

# Calculate the values to use as max and min for the volcano plots
min_log2FoldChange = max(min(expression_database$OSS_vs_LSS_paired_log2FoldChange),
                         min(expression_database$ESS_vs_LSS_paired_log2FoldChange))
max_log2FoldChange = min(max(expression_database$OSS_vs_LSS_paired_log2FoldChange),
                         max((expression_database %>% dplyr::filter(!is.na(ESS_vs_LSS_paired_padj)))$ESS_vs_LSS_paired_log2FoldChange))
max_paired_padj = max(min((expression_database %>% dplyr::filter(!is.na(OSS_vs_LSS_paired_padj)))$OSS_vs_LSS_paired_padj),
                      min((expression_database %>% dplyr::filter(!is.na(ESS_vs_LSS_paired_padj)))$ESS_vs_LSS_paired_padj))


# Graham
min_log2FoldChange = max(
  min(expression_database$OSS_vs_LSS_paired_log2FoldChange),
  min(expression_database$ESS_vs_LSS_paired_log2FoldChange),
  min(expression_database$OSS_vs_ESS_paired_log2FoldChange)
)

max_log2FoldChange = min(
  max(expression_database$OSS_vs_LSS_paired_log2FoldChange),
  max((expression_database %>% 
         dplyr::filter(!is.na(ESS_vs_LSS_paired_padj)))$ESS_vs_LSS_paired_log2FoldChange),
  max((expression_database %>% 
         dplyr::filter(!is.na(OSS_vs_ESS_paired_padj)))$OSS_vs_ESS_paired_log2FoldChange)
)

max_paired_padj = max(
  min((expression_database %>% 
         dplyr::filter(!is.na(OSS_vs_LSS_paired_padj)))$OSS_vs_LSS_paired_padj),
  min((expression_database %>% 
         dplyr::filter(!is.na(ESS_vs_LSS_paired_padj)))$ESS_vs_LSS_paired_padj),
  min((expression_database %>% 
         dplyr::filter(!is.na(OSS_vs_ESS_paired_padj)))$OSS_vs_ESS_paired_padj)
)



paired_res_OSS_vs_ESS <- expression_database %>% 
  dplyr::filter(!is.na(OSS_vs_ESS_paired_padj)) %>% 
  dplyr::mutate(color = ifelse(OSS_vs_ESS_paired_log2FoldChange <= -1 & OSS_vs_ESS_paired_padj <= 0.05,
                               "Down-regulated genes",
                               ifelse(OSS_vs_ESS_paired_log2FoldChange >= 1 & OSS_vs_ESS_paired_padj <= 0.05,
                                      "Up-regulated genes",
                                      "Not significant genes")),
                OSS_vs_ESS_paired_padj = ifelse(OSS_vs_ESS_paired_padj < max_paired_padj,
                                                max_paired_padj,
                                                OSS_vs_ESS_paired_padj),
                OSS_vs_ESS_paired_log2FoldChange = ifelse(OSS_vs_ESS_paired_log2FoldChange > max_log2FoldChange,
                                                          max_log2FoldChange,
                                                          ifelse(OSS_vs_ESS_paired_log2FoldChange < min_log2FoldChange,
                                                                 min_log2FoldChange,
                                                                 OSS_vs_ESS_paired_log2FoldChange))
                )
# Graham
OSS_wanted = c("LMEM2", "CYP1B1")
temp <- left_join(paired_res_OSS_vs_ESS, grch38, by = c("Target_ID" = "ensgene"))

OSS_sig_genes <- filter(temp, symbol %in% OSS_wanted) %>% rename("gene_name" = "symbol")
  

# Volcano Plot
ggplot(data = paired_res_OSS_vs_ESS,
       aes(x = OSS_vs_ESS_paired_log2FoldChange,
           y = -log10(OSS_vs_ESS_paired_padj),
           color = color
           )
       ) +
  geom_point(size = 1, shape = 16, alpha = 1) + 
  theme_minimal() +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             color = "black"
             ) + 
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed",
             color = "black") +
  labs(x = "log2(fold change)",
       y = "-log10(adj p-value)") +
  scale_color_manual(values = c("Down-regulated genes" = "darkblue", 
                                "Up-regulated genes" = "darkred",
                                "Not significant genes" = "grey")) +
  geom_text_repel(data = OSS_sig_genes, 
                  aes(label = gene_name ),
                  min.segment.length = 0,
                  max.overlaps = 10,
                  force = 1, 
                  force_pull = 0.2,
                  size = 5, 
                  box.padding = 0.1,
                  seed=7) +
  theme(
    axis.title = element_text(size = 30),   # Increase axis title size
    axis.text = element_text(size = 30),    # Increase axis text size
    legend.position = "right",              # Adjust legend position
    legend.title = element_blank(),         # Remove legend title
    legend.text = element_text(size = 20)   # Increase legend text size
  ) +
  guides(
    size = guide_legend(override.aes = list(shape = 16, color = "black")) #,
#     color = guide_legend(override.aes = list(size = 5))# Adjust size legend
  ) +
  theme(legend.position = "none")

