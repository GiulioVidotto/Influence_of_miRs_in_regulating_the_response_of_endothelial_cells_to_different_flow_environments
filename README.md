# Influence of miRs in regulating the response of endothelial cells to different flow environments

This repository contains the code that I wrote for the study that investigated the role of flow-responsive micro-RNAs (miRNAs) in endothelial cells (ECs). The findings of this project provided new insights into the molecular mechanisms governing endothelial function in vascular health and disease. Further details and results will be shared upon the publication of the associated research paper.

## Installation

1. Clone the repository and navigate to the repository:
    ```bash
    git clone https://github.com/GiulioVidotto/Influence_of_miRs_in_regulating_the_response_of_endothelial_cells_to_different_flow_environments.git
    ```
    
2. Navigate to the cloned repository:
    ```bash
    cd Influence_of_miRs_in_regulating_the_response_of_endothelial_cells_to_different_flow_environments.git
    ```
    
## Requisites for the analysis

### System dependencies:

Before running the script, ensure you have the following installed:
1. unzip (to extract `.zip` files):
    - For Ubuntu, run:
    ```bash
    sudo apt update
    sudo apt install unzip
    ```
    - For macOS:
    ```bash
    brew install unzip
    ```
    - For WSL (Windows Subsystem for Linux):
    ```bash
    sudo apt update
    sudo apt install unzip
    ```
2. wget (to download files from the web):
    - This is usually pre-installed on most systems. If not, install it using:
    ```bash
    sudo apt install wget
    ```
2. APAtrap
    - Follow the installation on the APAtrap website (Link to the website: https://sourceforge.net/p/apatrap/wiki/User%20Manual/#jump2, Paper: Congting Ye, Yuqi Long, Guoli Ji, Qingshun Quinn Li, Xiaohui Wu, APAtrap: identification and quantification of alternative polyadenylation sites from RNA-seq data, Bioinformatics, Volume 34, Issue 11, June 2018, Pages 1841–1849, https://doi.org/10.1093/bioinformatics/bty029).
4. ViennaRNA package:
    - follow the installation and tutorial on the ViennaRNA website (Link to the website: https://www.tbi.univie.ac.at/RNA/tutorial/, Paper: Lorenz, R., Bernhart, S.H., Höner zu Siederdissen, C. et al. ViennaRNA Package 2.0. Algorithms Mol Biol 6, 26 (2011). https://doi.org/10.1186/1748-7188-6-26)
    
### R Setup

1. Install R from CRAN (https://cran.r-project.org/mirrors.html) or use your package manager to install it.
2. Optional: download Rstudio (https://posit.co/downloads/)
3. Open R or Rstudio 
4. Set the current working directory to the cloned repository. To set the working directory on R or Rstudio, use:
    ```R
    setwd(list.files("Influence_of_miRs_in_regulating_the_response_of_endothelial_cells_to_different_flow_environments", full.names = TRUE)
    ```
5. This script installs and loads all the required R libraries and function for the project. You only need to run this script once, to make sure your environment is ready. To run the script on R or Rstudio:
    ```R
    source("load_functions_and_libraries_in_R_env.R")
    ```
    
## Download the data

This project depends on public available data and data obtained from lab experiments.
1. Download the public data:
    Make the scripts for downloading public data executable:
    ```bash
    chmod +x download_file.sh            # script for download a specific file
    chmod +x download_public_data.sh     # script containing the links of the files that will be downloaded in defined paths (this script uses download_file.sh to download the files)
    ```
    The script will automatically download files based on the current working directory, which will now be inside the cloned repository (after executing step 3)
    ```bash
    bash ./download_public_data.sh
    ```
    Specifically, this script creates a new directory called "project_data" in the current working directory. Inside this directory, there are other sub-directories to store different files from different websites.
    > **⚠️ Note:** At the moment the link to download the miRNA-mRNA interaction data from miRTarBase is not working as their website is currently unavailable (https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2022/php/index.php).
    ```bash
    project_data/                                        # Main directory for all project-related data
    ├── miR_database/                                    # Contains data from various miRNA databases
    │   ├── miRBase_database/                            # Data from miRBase
    │   │   ├── mature.fa                                # Raw miRNA sequences
    │   │   └── miRBase_human.txt                        # Filtered human miRNA IDs
    │   ├── targetScan_database/                         # TargetScan miRNA binding data
    │   │   └── Conserved_Site_Context_Scores.txt        # miRNA-mRNA interaction data with conserved binding sites
    │   │   └── Nonconserved_Site_Context_Scores.txt     # miRNA-mRNA interaction data with non-conserved binding sites
    │   ├── miRDB_database/                              # data from miRDB
    │   │   └── miRDB_v6.0_prediction_result.txt         # miRNA-mRNA interaction data from miRDB
    │   ├── miRTarBase_database/                         # data from miRTarBase
    │   │   └── miRTarBase_MTI.xlsx                      # miRNA-mRNA interaction data from miRTarBase
    │   └── tarBase_database/                            # data from TarBase
    │       └── Homo_sapiens_TarBase-v9.tsv              # miRNA-mRNA interaction data from TarBase
    └── genome_data/                                     # Contains genomic data files
        ├── bed_and_fasta_file/                          # Genome data in BED and FASTA formats
        └── chr_fasta_file/                              # Chromosome files in FASTA format
            ├── chr1.fa                                  # Chromosome 1 sequence
            ├── chr2.fa                                  # Chromosome 2 sequence
            ├── ...                                      # Other chromosomes (from chr3 to chr22)
            ├── chrX.fa                                  # Chromosome X sequence
            ├── chrY.fa                                  # Chromosome Y sequence
            └── chrM.fa                                  # Mitochondrial DNA sequence
    
2. Download the data obtained in the lab:
    > **⚠️ Note:** The lab data associated with this study will be released alongside the research paper (a folder with this data will be added to the project_data folder)

## Upload the data on R or Rstudio

For this project, four miRNA-mRNA interaction public repositories were considered (miRTarBase v9.0, targetScan v8.0, miRDB v6.0, tarBase v9.0). In the previous step, the miRNA-mRNA interaction files from these repositories were downloaded. The scripts contained in the folder called "upload_public_databases_R_scripts" load the public data on R or Rstudio. Alongside loading the data, the scripts will run a check on the miRNA notation based on the one used in miRBase (Release 22.1).
To upload the data, open R or Rstudio and follow the instructions on the scripts. To upload the data from miRDB and miRTarBase, it is required to use Biomart (https://useast.ensembl.org/info/data/biomart/index.html) to obtain columns not present in the original datasets. At the end of each script, a new file with the repository data and the controlled miRNA notation is stored in a folder called "final_outputs". This new folder is created as a subfolder of the miR_databases folder.

## miRNA-mRNA interaction analysis

### Calculate the minimum free energy (MFE) of each miRNA-mRNA interaction

To calculate the minimum free energy (MFE) of each interaction, we considered the sequences of miRNAs and the 3'UTR regions of mRNAs.
- miRNA Sequences: Obtained from miRBase (Release 22.1). The sequences are stored in the `mature.fa` file in `./project_data/miR_database/miRBase_database`.
- 3'UTR Regions of mRNAs: Extracted using APAtrap (https://sourceforge.net/p/apatrap/wiki/User%20Manual/#jump2).

Step-by-Step Process:
1. Extracting 3'UTR Regions:
    Use bedgraph files as input for the APAtrap tool to identify poly-A sites. Create the following folder `./project_data/genomic_data/APA_output` and save within it the APAtrap output file as `all_APA_output_chr.txt`.
2. Updating Gene Coordinates:
    Run the R script updating_gene_coordinates_based_on_poly_A_location.R in R or RStudio.
    - Input: APAtrap output file stored in `./project_data/genome_data/APA_output/all_APA_output_chr.txt`.
    - Output: Updated BED files with 3'UTR coordinates, saved in `./project_data/genome_data/bed_and_fasta_file/utr_with_poly_A_bed_file/`. The script generates multiple BED files, each corresponding to a specific chromosome with the following file naming format `3UTR_poly_A_regions_chr[chromosome number].bed`. Each BED file contains the updated 3'UTR coordinates for mRNAs based on their chromosomal location.
3. Obtaining DNA Sequences:
    Execute the shell script get_sequence.sh.
    - Input: BED files from Step 2.
    - Output: DNA sequences for the 3'UTR regions that will be stored in `./project_data/genome_data/bed_and_fasta_file/utr_with_poly_A_output/3UTR_with_poly_A_sequences_mRNA/mRNA_${file}`
    ```bash
    bash ./MFE_scripts/get_sequence.sh
    ```
4. Converting DNA to mRNA:
    Use the script DNA_to_mRNA.sh to convert DNA sequences into mRNA sequences.
    - Input: DNA sequence files from Step 3.
    - Output: mRNA sequences.
    Example command:3'UTR regions. To run the script:
    ```bash
    bash ./MFE_scripts/DNA_to_mRNA.sh <DNA_SEQUENCE_FILE> # Substitute <DNA_SEQUENCE_FILE> with the name of the file 3UTR_with_poly_A_sequences_on_chr[chromosome number]
    ```
5. Run the script get_MFE_function.sh to calculate MFE values using tools from the ViennaRNA package.
    Input:
    - <TARGET_QUERY_FILE>: Target sequences.
    - <GENE_SEQ_FILE>: mRNA sequences.
    - <MIRNA_SEQ_FILE>: miRNA sequences.
Logic:
For 3'UTR sequences shorter than 2000 bp, the RNAup function was used.
For longer sequences, a combination of RNAplfold and RNAplex functions was employed.
Example command:
    ```bash
    bash ./MFE_scripts/get_MFE_function.sh <TARGET_QUERY_FILE> <GENE_SEQ_FILE> <MIRNA_SEQ_FILE>
    ```


