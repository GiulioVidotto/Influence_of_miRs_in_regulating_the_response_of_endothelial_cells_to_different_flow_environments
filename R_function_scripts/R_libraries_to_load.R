# Libraries to install from Bioconductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Enrichedplot library
BiocManager::install("enrichplot")
library(enrichplot)
# ClusterProfiler library
BiocManager::install("clusterProfiler") # to install clusterProfiler I had to update ggplot2 with the following line of code -> remotes::install_version('ggplot2', version='3.3.6')
library(clusterProfiler)


# Libraries to install from CRAN
# UpSetR library
install.packages("UpSetR")
library(UpSetR)

# Readxl library
install.packages("readxl")
library(readxl)

# Igraph library
install.packages("igraph")
library(igraph)

# Uniprot library
install.packages("UniprotR")
library(UniprotR) 

# Rvest library
install.packages("rvest")
library(rvest)

# Ggh4x library
install.packages("ggh4x")
library(ggh4x)

# Ggnewscale library
install.packages("ggnewscale")
library(ggnewscale)

# Plotly library
install.packages("plotly")
library(plotly)

# VisNetwork library
install.packages("visNetwork")
library(visNetwork)

# Dplyr library
install.packages("dplyr")
library(dplyr)

# Purrr library
install.packages("purrr")
library(purrr)

# Tidyverse library
install.packages("tidyverse")
library(tidyverse)

# Ggpubr library
install.packages("ggpubr")
library(ggpubr)

# Reshape2 library
install.packages("reshape2")
library(reshape2)

# Vcd library
install.packages("vcd")
library(vcd)

# Shapviz library
install.packages("shapviz")
library(shapviz)

