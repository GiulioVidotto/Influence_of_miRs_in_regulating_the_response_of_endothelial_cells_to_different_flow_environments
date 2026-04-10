# IPA workflow

# Load IPA tables for 3 contrast analysis and make excel containing tabs 
#Should need only to change these parts to run it twice, once for mRNA, once for protein

# GRS define stuff

# outfile_name <- "IPA_results_3_contrasts.xlsx"
# ipa_dir <- "~/data/Steve_White/results/Giulio/2026/IPA_output/mRNA_3contrasts_all_IPA_tabs/"
# prefix <- "mRNA_"

outfile_name <- "prot_IPA_results_3_contrasts.xlsx"
ipa_dir <- "" # select the directory
prefix <- "prot_"


contr <- c("OSS_vs_LSS", "ESS_vs_LSS", "OSS_vs_ESS")
shcontr <- gsub("_vs_", "v", contr)

name_lookup <- c("CP", "UR", "CN", "DF", "RE")
filenames <- c("bar_plot", "Upstream_regulators", "Causal_Networks", "Diseases_Functions", "Regulator_Effects")
names(name_lookup) <- filenames

for (fn in filenames)
  cat(name_lookup[fn], " ", fn, "\n")

# Now run it

wb  <- openxlsx::createWorkbook()

for (j in 1:5) {
  IPAtype <- filenames[j]
  for (i in 1:3) {
    ctr <- contr[i]
    sh <- shcontr[i]
    infile <- paste0(prefix, ctr, "_", IPAtype, "_table.txt")
    IPAabbr <- name_lookup[IPAtype]
    outtab <- paste0(sh, "_", IPAabbr)
    cat(infile, " ", outtab, "\n")
    
    # have sorted it out so all files exist
    data <- read_tsv(file.path(ipa_dir, infile), skip = 2)
    
    ##cat("hello"); browser()
    if (IPAabbr == "CP")
      data <- select(data, -6) # get rid of '...6' column
    openxlsx::addWorksheet(wb, outtab)
    openxlsx::writeData(wb, outtab, data)
    
  }
}

outfile <- file.path(ipa_dir, outfile_name)
openxlsx::saveWorkbook(wb, outfile, overwrite = TRUE)
