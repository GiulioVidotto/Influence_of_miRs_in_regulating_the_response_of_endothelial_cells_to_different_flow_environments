## read excel file and check 

dat = read_excel("~/excel_proteomics_first file.xlsx", sheet = 1 )

for (col in grep("Abundance", names(dat), val = T)) {
  show(hist(log2(dat[[col]]), main = col, breaks = 50))
}

# Another sanity check is that (A/B) = (A/C)/(B/C)
# So AR(OSS vs ESS ) = AR(OSS vs LSS) / AR(ESS vs LSS)

dat <- mutate(dat, Ab_OSS_vs_LSS_over_ESS_vs_LSS = protein_OSS_vs_LSS_Abundance_Ratio / protein_ESS_vs_LSS_Abundance_Ratio)

dat2 <- filter(dat, protein_OSS_vs_LSS_Abundance_Ratio != 0.01) %>% filter(protein_OSS_vs_LSS_Abundance_Ratio != 100) %>%
   filter(protein_ESS_vs_LSS_Abundance_Ratio != 0.01) %>%  filter(protein_ESS_vs_LSS_Abundance_Ratio != 100) %>%
   filter(protein_OSS_vs_ESS_Abundance_Ratio != 0.01) %>%  filter(protein_OSS_vs_ESS_Abundance_Ratio != 100)

ggplot(dat2, aes(x = protein_OSS_vs_ESS_Abundance_Ratio, y = Ab_OSS_vs_LSS_over_ESS_vs_LSS)) + 
  geom_point(size = 0.1) + theme_bw(base_size = 14) 

ggplot(dat2, aes(x = protein_OSS_vs_ESS_Abundance_Ratio, y = Ab_OSS_vs_LSS_over_ESS_vs_LSS)) + 
  geom_point(size = 0.1) + theme_bw(base_size = 14) + scale_x_log10() + scale_y_log10()


# output log trans version. 

# Either remove (set to NA) values for FC and pval where FC is 100 or 0.01 (output file "proteomics_3contrasts_log2_table_rm100.txt") or leave them in ("proteomics_3contrasts_log2_table_with100.txt")

dat = read_excel("~/excel_proteomics_first file.xlsx", sheet = 1 )

rm_outliers = TRUE

if (rm_outliers){
  outfile = "proteomics_3contrasts_log2_table_rm100.txt"
} else {
  outfile = "proteomics_3contrasts_log2_table_with100.txt"
}

filterfn <- function(d, v) {
  pv = gsub("_Abundance_Ratio", "_adj_pvalue", v)
  for (to_rm in c(0.01,100)) {
    idx <- which(d[[v]] == to_rm)
    if (length(idx) > 0 ) {
      d[[v]][idx] <- NA 
    }
  }
  # also set corr values to be NA
  
  idx <- which(is.na(d[[v]]))
  d[[pv]][idx] <- NA   
  return(d)
}

abcols = grep("Abundance", names(dat), val = T)

if(rm_outliers) {
  for (col in abcols) {
    dat = filterfn(dat, col)
  }
}

for (col in abcols) {
  dat[[col]] <- log2(dat[[col]]) 
}

names(dat) <- gsub("Ratio", "Log_Ratio", names(dat))
for (col in grep("Abundance", names(dat), val = T)) {
  show(hist(dat[[col]], main = col, breaks = 50))
}
#for (col in grep("pvalue", names(dat), val = T)) {
#  show(hist(dat[[col]], main = col, breaks = 50))
#}
# how mmny NA
colSums(is.na(dat))

# Also ooutput the pval and expression just for the DE proteins (without the 100fold ones again) for DE proteins in OSS_vs_ESS

outdat <- select(dat, Protein_ID, Target_ID, protein_OSS_vs_ESS_Abundance_Log_Ratio,  protein_OSS_vs_ESS_adj_pvalue) %>%
  rename(protein_OSS_vs_ESS_Abundance_Log2_Ratio = protein_OSS_vs_ESS_Abundance_Log_Ratio) %>%
  rename(ensgene = Target_ID) %>%
  filter(! is.na(ensgene)) %>%
  filter(abs(protein_OSS_vs_ESS_Abundance_Log2_Ratio) > 1 & protein_OSS_vs_ESS_adj_pvalue < 0.05) 

library(annotables)
symannot <- select(grch38, ensgene, symbol) %>% unique()
outdat <- left_join(outdat, symannot) %>%
  relocate(symbol, .after = ensgene) %>%
  arrange(desc(protein_OSS_vs_ESS_Abundance_Log2_Ratio))

write_csv(outdat, file = file.path(ipa_dir, "DE_proteins_OSS_vs_ESS.csv"))

