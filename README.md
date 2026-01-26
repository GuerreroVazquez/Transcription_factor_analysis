# Transcription Factor Enrichment Analysis
Transcription factor analysis for RNAseq data.
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- Add DOI badge here if archiving via Zenodo -->

## Overview

This repository contains the computational pipeline for identifying **Enriched Transcription Factors (ETFs)** from RNA-Seq data. It utilises the **decoupleR** package to infer transcription factor activity based on gene expression changes, leveraging the **CollecTRI** gene regulatory network resource.

This workflow corresponds to the **Transcription Factor Enrichment** module described in **Chapter 4** of the PhD thesis *"Computational approaches for therapeutic target discovery to ameliorate muscle wasting during ageing"* (University of Galway, 2025).

## Methodology

As detailed in **Section 4.2.1** of the thesis:

1.  **Input Data:** Accepts Differential Expression statistics (typically Wald test statistics or Log2 Fold Changes) from DESeq2.
2.  **Network Source:** Utilises the **CollecTRI** database (filtered for Human interactions), a comprehensive resource of signed TF-gene interactions.
3.  **Enrichment Algorithm:** Applies the **Univariate Linear Model (ULM)** method via `decoupleR`. This fits a linear model for each TF to estimate its activity score based on the expression of its target genes.
4.  **Filtering:** TFs are designated as "Enriched" (ETF) if they achieve an enrichment p-value < 0.01.

## Dependencies

*   Python 3.11+
*   [decoupleR](https://github.com/saezlab/decoupler-py)
*   [OmniPath](https://omnipathdb.org/) (for retrieving CollecTRI)
*   Pandas
*   NumPy

## Usage

The environment is build from the [Dockerfile](./container/Dockerfile) and the container are managed using `docker compose`

To build and run the environment, run 'docker compose up' from the repository folder. The first time, it may take some minutes to build the image.

To access jupyter go to: 'localhost:8860'

### 1. Installation
```bash
pip install decoupler pandas numpy omnipath
```


### 2. Running the Analysis

The main script takes a CSV file containing gene symbols and their corresponding statistics (e.g., `stat` from DESeq2).

```python
import decoupler as dc
import pandas as pd

# 1. Load Data
# Data should have columns: 'Gene' and 'stat'
df = pd.read_csv('deseq2_results.csv', index_col=0) 
mat = df[['stat']].T # Transpose to match decoupler format (samples x genes)

# 2. Retrieve CollecTRI Network
net = dc.get_collectri(organism='human', split_complexes=False)

# 3. Run ULM (Univariate Linear Model)
results = dc.run_ulm(
    mat=mat,
    net=net,
    source='source',
    target='target',
    weight='weight',
    verbose=True
)

# 4. Extract and Filter Results
# p-values are stored in results[1]
tf_acts = results[1]
significant_tfs = tf_acts[tf_acts < 0.01].dropna(axis=1)

print(f"Identified {significant_tfs.shape[1]} significant TFs.")
significant_tfs.to_csv("enriched_tfs.csv")
```

## Outputs

The analysis generates a matrix of enrichment scores and p-values. In the context of the thesis (miRKat pipeline), the output is used to weight nodes in the regulatory network:
*   **Positive Score:** Indicates TF activation.
*   **Negative Score:** Indicates TF inhibition.



