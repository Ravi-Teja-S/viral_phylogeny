# 🧬 Viral Phylogeny: SARS-CoV-2 Genomic Analysis

[![Python](https://img.shields.io/badge/Python-3.x-blue.svg)](https://www.python.org/)
[![Biopython](https://img.shields.io/badge/Biopython-1.78%2B-green.svg)](https://biopython.org/)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/Ravi-Teja-S/viral_phylogeny/graphs/commit-activity)

##  Project Overview

This project performs a comprehensive **phylogenomic analysis of SARS-CoV-2 (COVID-19) genome sequences** using Python. By leveraging the **Biopython** library, this workflow automates the retrieval of viral sequences, aligns them to identify mutations, and constructs phylogenetic trees to visualize evolutionary relationships between different viral isolates.

This tool is useful for understanding the genetic divergence of SARS-CoV-2 variants and tracking evolutionary lineages.

##  Repository Link
[https://github.com/Ravi-Teja-S/viral_phylogeny](https://github.com/Ravi-Teja-S/viral_phylogeny)

##  Objectives

- **Sequence Retrieval:** Automate fetching of SARS-CoV-2 genomes from the NCBI Nucleotide database.
- **Multiple Sequence Alignment (MSA):** Align sequences to identify conserved regions and variations.
- **Phylogenetics:** Construct evolutionary trees using distance-based methods (e.g., UPGMA, Neighbor-Joining).
- **Visualization:** Render and analyze the phylogenetic tree to interpret viral evolution.

##  Tech Stack

- **Language:** Python 3
- **Bioinformatics:** Biopython (`Bio.Entrez`, `Bio.Phylo`, `Bio.Align`)
- **Data Source:** NCBI GenBank (FASTA format)
- **Visualization:** Matplotlib / Biopython Phylo module

##  Getting Started

Follow these instructions to set up the project on your local machine.

### Prerequisites

Ensure you have Python 12 installed. You can verify this by running:
```bash
python --version
```

### Clone the Repository
```bash
git clone https://github.com/Ravi-Teja-S/viral_phylogeny.git
cd viral_phylogeny
```

### Create a virtual environment 
```bash
python -m venv venv
venv\Scripts\activate
```

### Install Requirements
```bash
pip install -r requirements.txt
```

##  Steps for editing the repository [Required before making changes]

1. Update the main branch and create a new feature branch
   ```bash
   git checkout main
   git pull origin main
   git checkout -b feature-update-your_branch_name
   ```
   
2. Now you can make changes in the code as needed.

3. Stage, commit, and push your changes in your branch
    ```bash
    git add .
    git commit -m "Your commit message"
    git push origin feature-update-your_branch_name
    ```
    
4. Go to the repository on GitHub and click “Compare & pull request” to submit your changes for review.

### NOTE: do not push to main branch.
