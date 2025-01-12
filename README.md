# Influence of miRs in regulating the response of endothelial cells to different flow environments

This repository contains the code that I wrote study investigated the role of flow-responsive micro-RNAs (miRNAs) in endothelial cells (ECs). The findings of this project provided new insights into the molecular mechanisms governing endothelial function in vascular health and disease. Further details and results will be shared upon the publication of the associated research paper.

## Installation

1. Clone the repository and navigate to the repository:
    ```bash
    git clone https://github.com/GiulioVidotto/Influence_of_miRs_in_regulating_the_response_of_endothelial_cells_to_different_flow_environments.git
    ```
2. Download the data:
    This project depends on public available data and data obtained from lab experiments.
    2.1 Download the data obtained in the lab:
    > **⚠️ Note:** The lab data associated with this study will be released alongside the research paper.
    2.2 Download the public data:
    ```bash
    bash download_public_data.sh
    ```
    This script creates a new directory in your home directory (`$HOME/project_data`). Inside this directory, there are other sub-directories to store different files from different websites.
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
    ├── genome_data/                                     # Contains genomic data files
    │   ├── bed_and_fasta_file/                          # Genome data in BED and FASTA formats
    │   └── chr_fasta_file/                              # Chromosome files in FASTA format
    │       ├── chr1.fa                                  # Chromosome 1 sequence
    │       ├── chr2.fa                                  # Chromosome 2 sequence
    │       ├── ...                                      # Other chromosomes (from chr3 to chr22)
    │       ├── chrX.fa                                  # Chromosome X sequence
    │       ├── chrY.fa                                  # Chromosome Y sequence
    │       └── chrM.fa                                  # Mitochondrial DNA sequence
    ```
4. 
    
