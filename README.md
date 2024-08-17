# A Workflow for Merging Logical Gene Regulatory Network Models

## Overview

This repository contains the code and resources for merging gene regulatory network (GRN) models as described in the manuscript **"Merging Logical Models: An Application in Acute Myeloid Leukemia Modeling"**. The repository provides a complete workflow for merging logic models, enhancing our understanding of complex biological systems, and demonstrate its application in Acute Myeloid Leukemia (AML).

The workflow presented in this repository involves sequential steps:

1. **Identifying Candidate Models**: This step involves reviewing existing literature, repositories, and databases to identify models with shared components, such as shared genes in the GRNs.
2. **Standardizing and Annotating Models**: Standardizing the gene names using international standards like HGNC approved symbols and ensuring models are reproducible.
3. **Reproducing Selected Models**: Verifying that selected models replicate the behaviors described in their original publications.
4. **Merging Models**: Using the provided code to merge models with different logical combination methods (`OR`, `AND`, `Inhibitor Wins`).
5. **Evaluating the Merged Models**: Comparing the predictive accuracy and robustness of the merged models against the original models and applying the merged models to new, untested scenarios. 

## Repository Structure

This repository is organized according to the workflow described in the manuscript:

1. **`Standardizing & Annotating Models`**:
    - [Standardization](Standardizing%20and%20annotating%20models/Standardization.ipynb): Converting models in text file to SBML-qual format.
    - [Annotation](Standardizing%20and%20annotating%20models/Annotation.ipynb): Fetching HGNC gene symbols for input SBML-qual models and updating the gene names after manual verification.
2. **`Reproducing Selected Models`**: 
    - Reproducibility check for each collected model, including:
        - [`Bonzanni2013`](Reproducing%20selected%20models/Bonzanni2013)
        - [`Ikonomi2020`](Reproducing%20selected%20models/Ikonomi2020)
        - [`Krumsiek2011`](Reproducing%20selected%20models/Krumsiek2011)
        - [`Palma2021`](Reproducing%20selected%20models/Palma2021) *(As Fig S3, S4 in the manuscript)*
3. **`Composing Models`**: 
    - [Merge logical models](Composing%20models/Merge%20logical%20models.ipynb): Merging logical models, including the OR, AND, and Inhibitor Wins methods.
    - Support model input in:
        - Text files using a  EBNF description
        - SBML-qual files
4. **`Evaluating the Merged Model`**: 
    - Contains notebooks for various evaluation tasks:
        - [**Coverage**](Evaluating%20the%20merged%20model/Coverage.ipynb): Assessing the coverage of AML patients with each mutation profiles using BeatAML, TCGA, AMLSG and cBioPortal data.  *(Fig S7)*
        - [**Stable States Heatmap**](Evaluating%20the%20merged%20model/Stable%20states%20heatmap.ipynb): Visualizing stable states of the merged models and clustering them with individual models. *(Fig 3)*
        - [**Correlation with HSC Expression**](Evaluating%20the%20merged%20model/Correlation%20with%20HSC%20expression.ipynb): Analyzing the correlation of model predictions with hematopoietic stem cell expression data. *(Fig 4, Fig S2)*
        - **Correlation with Clinical Outcomes**: Separate notebooks for evaluating correlations with different clinical indicators/datasets:
            - Blast percentages from the BeatAML data
                - [Using approach similar to Palma et al.](Evaluating%20the%20merged%20model/Correlation%20with%20clinical%20outcome_BeatAML_Palma%20approach.ipynb) *(Fig 5a,b, Fig S5)*
                - [Using all mutations](Evaluating%20the%20merged%20model/Correlation%20with%20clinical%20outcome_BeatAML_all%20mutation.ipynb) *(Fig 5c,d, Fig S6)*
            - [Blast percentages from the TCGA data](Evaluating%20the%20merged%20model/Correlation%20with%20clinical%20outcome_TCGA.ipynb)
            - [Hazard ratio for death from the AMLSG data](Evaluating%20the%20merged%20model/Correlation%20with%20clinical%20outcome_AMLSG.ipynb)

- [**`Data`**](Data): Contains datasets used for model evaluation.
   
- [**`Models`**](Models): Stores models in text files and SBML-qual format.


## How to Use This Repository

1. **Clone the repository**:
   ```bash
   git clone https://github.com/IlyaLab/LogicModelMerger.git

2. **Install and open the CoLoMoTo Notebook** (Optional):
   Please refer to the usage guide on their website.

3. **Run notebooks:**
   Navigate to the relevant directory and open the Jupyter notebooks using the CoLoMoTo notebook or your preferred Jupyter environment.


### CoLoMoTo Interactive Notebook

All Jupyter notebooks in this repository were conducted using the CoLoMoTo Interactive Notebook with the Docker image `colomoto/colomoto-docker:2024-03-01`. The CoLoMoTo notebook provides a unified environment to edit, execute, share, and reproduce analyses of qualitative models of biological networks.

For more information about CoLoMoTo Interactive Notebook, visit [CoLoMoTo](http://www.colomoto.org/notebook/).
