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
1. **`unzip`** (to extract `.zip` files):
    - For **Ubuntu**, run:
    ```bash
    sudo apt update
    sudo apt install unzip
    ```
    - For **macOS**:
    ```bash
    brew install unzip
    ```
    - For **WSL (Windows Subsystem for Linux)**:
    ```bash
    sudo apt update
    sudo apt install unzip
    ```

2. **`wget`** (to download files from the web):
    - This is usually pre-installed on most systems. If not, install it using:
    ```bash
    sudo apt install wget
    ```
    
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

To calculate the minimum free energy (MFE) of each interaction, the sequences of the miRNAs and of the 3'UTR regions on the mRNAs were considered. The sequences of the miRNAs were obtained from miRBAse (Release 22.1) while the sequences of the 3'UTR regions of the mRNAs were obtained using a software called APAtrap (https://sourceforge.net/p/apatrap/wiki/User%20Manual/). To obtain the wanted regions, we used the bedgraph files as input and followed the instruction on the APAtrap website. After obtaining the sequences of both biological molecules, we used the scripts contained in the folder called MFE_scripts to calculate the MFE of each iteraction. Specifically:
1. Use the bedgraph files as inputs for the APAtrap function and save the output file as "all_APA_output_chr.txt". This file must be saved in a new subfolder of "genome_data" folder called "APA_output".
2. Open R or Rstudio and run the script "updating_gene_coordinates_based_on_poly_A_location.R". This script loads the APA output file (generated using APAtrap) and processes it to update the stop coordinates of each mRNA based on the furthest poly-A site. The script outputs a separate BED graph file for each chromosome, storing the updated coordinates of the 3'UTR region for all mRNAs on that chromosome. These bedgraph files are stored in a new folder ("./project_data/genome_data/bed_and_fasta_file/utr_with_poly_A_bed_file") where each file has the following name: "3UTR_poly_A_regions_chr[chromosome number]".
3. Run the script called "get_sequence.sh" from the shell. This script takes as input the bed files obtained at the previous step to obtain the DNA sequences corresponding to the 3'UTR region coordinates of each mRNA. To run this script:
    ```bash
    bash ./MFE_scripts/get_sequence.sh
    ```
4. Run the script called "DNA_to_mRNA.sh" to tranform the DNA sequences obtained at the previous step into mRNA sequences. This script accepts as input the name of the files that stores the DNA sequences of the 3'UTR regions. To run the script:
    ```bash
    bash ./MFE_scripts/DNA_to_mRNA.sh <DNA_SEQUENCE_FILE> # Substitute <DNA_SEQUENCE_FILE> with the name of the file 3UTR_with_poly_A_sequences_on_chr[chromosome number]
    ```
5. Lastly, run the script "get_MFE_function.sh" to obtain the MFE values of each interaction. To run the script:
    ```bash
    bash ./MFE_scripts/get_MFE_function.sh <TARGET_QUERY_FILE> <GENE_SEQ_FILE> <MIRNA_SEQ_FILE>
    ```


