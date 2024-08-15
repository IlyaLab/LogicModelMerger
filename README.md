# A Workflow for Merging Logical Gene Regulatory Network Models

## Overview

This repository contains the code and resources for merging gene regulatory network (GRN) models as described in the manuscript **"Merging Logical Models: An Application in Acute Myeloid Leukemia Modeling"**. The repository provides a complete workflow for merging logic models, enhancing our understanding of complex biological systems, and demonstrate its application in Acute Myeloid Leukemia (AML).

## Repository Structure

This repository is organized according to the workflow described in the manuscript:

1. **`Standardizing & Annotating Models`**:
    - To fetch HGNC gene symbols for input SBML-qual models and update the gene names after manual verification.
2. **`Reproducing Selected Models`**: 
    - Folders contain reproducibility check notebooks for each collected model, including:
        - `Bonzanni2013`
        - `Ikonomi2020`
        - `Krumsiek2011`
        - `Latini2023`
        - `Palma2021`
3. **`Composing Models`**: 
    - Contains code for merging logical models, including the notebook that merges models using the OR, AND, and Inhibitor Wins methods.
    - Support model input in:
        - Text files using a  EBNF description
        - SBML-qual files
4. **`Evaluating the Merged Model`**: 
    - Contains notebooks for various evaluation tasks:
        - **Coverage**: Assessing the coverage of AML patients with each mutation profiles using BeatAML, TCGA, AMLSG and cBioPortal data.
        - **Stable States Heatmap**: Visualizing stable states of the merged models and clustering them with individual models.
        - **Correlation with HSC Expression**: Analyzing the correlation of model predictions with hematopoietic stem cell expression data.
        - **Correlation with Clinical Outcomes**: Separate notebooks for evaluating correlations with different clinical indicators/datasets:
            - Blast percentages from the BeatAML data (using different approaches)
            - Blast percentages from the TCGA data
            - Hazard ratio for death from the AMLSG data

- **`Data`**: Contains datasets used for model evaluation.
   
- **`Models`**: Stores models in text files and SBML-qual format.

## Workflow and Tools

The workflow presented in this repository involves sequential steps:

1. **Identifying Candidate Models**: This step involves reviewing existing literature, repositories, and databases to identify models with shared components, such as shared genes in the GRNs.
2. **Standardizing and Annotating Models**: Standardizing the gene names using international standards like HGNC approved symbols and ensuring models are reproducible.
3. **Reproducing Selected Models**: Verifying that selected models replicate the behaviors described in their original publications.
4. **Merging Models**: Using the provided code to merge models with different logical combination methods (`OR`, `AND`, `Inhibitor Wins`).
5. **Evaluating the Merged Models**: Comparing the predictive accuracy and robustness of the merged models against the original models and applying the merged models to new, untested scenarios.

### CoLoMoTo Interactive Notebook

All Jupyter notebooks in this repository were conducted using the CoLoMoTo Interactive Notebook with the Docker image `colomoto/colomoto-docker:2024-03-01`. The CoLoMoTo notebook provides a unified environment to edit, execute, share, and reproduce analyses of qualitative models of biological networks.

For more information about CoLoMoTo Interactive Notebook, visit [CoLoMoTo](http://www.colomoto.org/notebook/).

## How to Use This Repository

1. **Clone the repository**:
   ```bash
   git clone https://github.com/IlyaLab/LogicModelMerger.git

2. **Install and open the CoLoMoTo Notebook** (Optional):
   Please refer to the usage guide on their website.

3. **Run notebooks:**
   Navigate to the relevant directory and open the Jupyter notebooks using the CoLoMoTo notebook or your preferred Jupyter environment.

### Reproducing the Figures in the Manuscript

Figures from the manuscript can be reproduced using the following notebooks:
- Figure 3: `Evaluating the merged model/Stable states heatmap.ipynb`
- Figure 4: `Evaluating the merged model/Correlation with HSC expression.ipynb`
- Figure 5a,b: `Evaluating the merged model/Correlation with clinical outcome_BeatAML_Palma approach.ipynb`
- Figure 5c,d: `Evaluating the merged model/Correlation with clinical outcome_all mutation.ipynb`
- Figure S1: `Evaluating the merged model/Stable states heatmap.ipynb`
- Figure S2: `Evaluating the merged model/Correlation with HSC expression.ipynb`
- Figure S3, S4: `Reproducing selected models/Palma2021/Palma2021.ipynb`
- Figure S5: `Evaluating the merged model/Correlation with clinical outcome_BeatAML_Palma approach.ipynb`
- Figure S6: `Evaluating the merged model/Correlation with clinical outcome_all mutation.ipynb`
- Figure S7: `Evaluating the merged model/Coverage.ipynb`
