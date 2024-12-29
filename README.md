# A Workflow for Merging Logical Gene Regulatory Network Models

## Overview

This repository provides a semi-automated workflow for merging logic models, and demonstrate its application using a set of Acute Myeloid Leukemia (AML) models.  

The workflow is described in the manuscript **"LM-Merger: A workflow for merging logical models with an application to gene regulation"**. ([link](https://www.biorxiv.org/content/10.1101/2024.09.13.612961v3))
>Li LX, Aguilar B, Gennari J, Qin G. LM-Merger: A workflow for merging logical models with an application to gene regulation. bioRxiv. 2024 Jan 1;2024.09.13.612961. 


The workflow involves sequential steps:

1. **Finding Models**: This step involves reviewing existing literature, repositories, and databases to identify models with shared components.
2. **Standardizing and Annotating Models**: Convert models to SBML-qual format, and annotate gene names using HGNC approved symbols.
3. **Reproducing Selected Models**: Verifying that selected models replicate the behaviors described in their original resources.
4. **Merging Models**: Using the provided tool to automatically merge models with different logical combination methods (`OR`, `AND`, `Inhibitor Wins`).
5. **Evaluating the Merged Models**: Comparing the predictive accuracy and robustness of the merged models against the original models and applying the merged models to new, untested scenarios. 

## Repository Structure

This repository is organized according to the workflow described in the manuscript:

1. **`Standardizing & Annotating Models`**:
    - [Standardization](Standardizing%20and%20annotating%20models/Convert%20model%20in%20text%20file%20to%20SBML-qual.ipynb): Converting models in text file to SBML-qual format.
    - [Annotation](Standardizing%20and%20annotating%20models/Standardize%20gene%20names%20to%20HGNC%20symbol.ipynb): Fetching HGNC gene symbols for input SBML-qual models and updating the gene names after manual verification.
2. **`Reproducing Selected Models`**: 
    - Reproducibility check for each collected model, including:
        - [Bonzanni2013](Reproducing%20selected%20models/Bonzanni2013/Bonzanni2013.ipynb)
        - [Ikonomi2020](Reproducing%20selected%20models/Ikonomi2020/Ikonomi2020.ipynb)
        - [Krumsiek2011](Reproducing%20selected%20models/Krumsiek2011/Krumsiek2011.ipynb)
        - [Palma2021](Reproducing%20selected%20models/Palma2021/Palma2021.ipynb) *(As Table S3 in the manuscript)*
3. **`Composing Models`**: 
    - [Merge logical models](Composing%20models/Merge%20logical%20models.ipynb): Automatically merging logical models, including the OR, AND, and Inhibitor Wins methods.
    - Support models in:
        - Text files using a EBNF description as in `Boolnet`
        - SBML-qual files
4. **`Evaluating the Merged Model`**: 
    - [Functions](Evaluating%20the%20merged%20model/Helper%20functions.ipynb): Provides some helper functions for evaluating logical models.
    - Contains notebooks for various evaluation tasks:
        - [Coverage](Evaluating%20the%20merged%20model/Coverage.ipynb): Assessing the coverage of AML patients with each mutation profiles using BeatAML, TCGA, AMLSG and cBioPortal data.
        - Stable States Heatmap: Visualizing stable states of the merged models and clustering them with individual models.
            - [Stable States Heatmap - Asynchronous update](Evaluating%20the%20merged%20model/Stable%20states%20heatmap_asynchronous.ipynb) *(Fig 2B, S1)*
            - [Stable States Heatmap - Synchronous update](Evaluating%20the%20merged%20model/Stable%20states%20heatmap_synchronous.ipynb) *(Fig 3B, S1)*
        - [Correlation with HSC Expression](Evaluating%20the%20merged%20model/Correlation%20with%20HSC%20expression.ipynb): Analyzing the correlation of model predictions with hematopoietic stem cell expression data. *(Fig 2C-E, Fig S2)*
        - Correlation with Clinical Outcomes: Separate notebooks for evaluating correlations with different clinical indicators/datasets:
            - Blast percentages from the BeatAML data
                - [Using approach similar to Palma et al.](Evaluating%20the%20merged%20model/Correlation%20with%20clinical%20outcome_BeatAML_Palma%20approach.ipynb) *(Fig 3C-D, Fig S4)*
                - [Using all mutations](Evaluating%20the%20merged%20model/Correlation%20with%20clinical%20outcome_BeatAML_all%20mutation.ipynb) *(Fig 3E-F, Fig S4)*
            - [Blast percentages from the TCGA data](Evaluating%20the%20merged%20model/Correlation%20with%20clinical%20outcome_TCGA.ipynb) *(TableS3)*
            - [Hazard ratio for death from the AMLSG data](Evaluating%20the%20merged%20model/Correlation%20with%20clinical%20outcome_AMLSG.ipynb)

- [**`Data`**](Data): Contains datasets used for model evaluation.
   
- [**`Models`**](Models): Stores models in text files and SBML-qual format.


## How to Use This Repository

1. **Clone the repository**:
   ```bash
   git clone https://github.com/IlyaLab/LogicModelMerger.git
    ```
2. **Install and open the CoLoMoTo Notebook** :  

    While running the notebooks in a standard Jupyter environment is possible, we recommend using the CoLoMoTo Interactive Notebook for better compatibility and streamlined integration with qualitative modeling tools.  

   ```bash
   pip install -U colomoto-docker 
   colomoto-docker # start the notebook in a terminal
    ```
   Please refer to the usage guide on their [website](https://colomoto.github.io/colomoto-docker/).

3. **Run notebooks:**  
   Navigate to the relevant directory and open the Jupyter notebooks using the CoLoMoTo notebook or your preferred Jupyter environment.


### CoLoMoTo Interactive Notebook

Jupyter notebooks in this repository were conducted using the CoLoMoTo Interactive Notebook with Docker image `colomoto/colomoto-docker:2024-03-01`.
  
 The CoLoMoTo notebook provides a unified environment to edit, execute, share, and reproduce analyses of qualitative models of biological networks.

For more information, visit [CoLoMoTo](http://www.colomoto.org/notebook/).
