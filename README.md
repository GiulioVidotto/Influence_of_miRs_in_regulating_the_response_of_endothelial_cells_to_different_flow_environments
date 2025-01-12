# Influence of miRs in regulating the response of endothelial cells to different flow environments

This study investigated the role of flow-responsive micro-RNAs (miRNAs) in endothelial cells (ECs). The findings of this project provided new insights into the molecular mechanisms governing endothelial function in vascular health and disease. Further details and results will be shared upon the publication of the associated research paper.

## Installation

1. Clone the repository and navigate to the repository:
    ```bash
    git clone https://github.com/GiulioVidotto/Influence_of_miRs_in_regulating_the_response_of_endothelial_cells_to_different_flow_environments.git
    ```
2. Download the public data:
    ```bash
    bash download_public_data.sh
    ```
   This script creates a new directory in the home called "project_data". Inside this directory, there are other sub-directories to store different files from different websites.
   - **`miR_database/`**
    - **Purpose**: Contains data from various miRNA databases used in the project.
    - **Subfolders**:
      - `miRBase_database/`: Stores the miRNA data from **miRBase**. The `mature.fa` file is downloaded and processed to include only human miRNA IDs, which are saved in the `miRBase_human.txt` file.
      - `targetScan_database/`: Contains the files downloaded from **TargetScan**, including context scores for conserved and non-conserved sites.
      - `miRDB_database/`: Stores the miRNA prediction results from **miRDB**.
      - `miRTarBase_database/`: Stores the miRNA-mRNA interaction data from **miRTarBase**.
      - `tarBase_database/`: Stores the miRNA-mRNA interaction data from **TarBase**.

- **`genome_data/`**
    - **Purpose**: Contains genomic data files.
    - **Subfolders**:
      - `bed_and_fasta_file/`: Stores genome data in **FASTA** format and other genomic files.
        - `chr_fasta_file/`: Contains the chromosome files (1-22, X, Y, and mitochondrial DNA) in **FASTA** format. These are downloaded from UCSC and stored here.

> **⚠️ Note:** The databases associated with this study will be released alongside the research paper.
