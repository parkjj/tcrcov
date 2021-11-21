# T cell receptor repertoire signatures associated with COVID-19 severity

This repository contains code accompanying the manuscript "T cell receptor repertoire signatures associated with COVID-19 severity."

## Project Structure
```
.
├── ml
│   ├── ML_models.py
│   └── significant_features.py
├── motif
│   ├── gliph2
│   ├── kmer
│   └── olga
├── stats
└── transcriptome

```

## Dependencies

```
Python 3.8.6
numpy 1.12.1
pandas 1.2.4
scikit-learn 0.23.1
Plotly 5.1.0
R 4.1.0
Immunarch 0.6.6
Seurat 4.0.4
GLIPH2 0.01 
OLGA 1.2.4
```
## Pipeline

### Stats
This folder contains code for calculating diversity metrics, CDR3 gene usage analyses, VDJ gene usage comparisons, and statistics. Most of the R scripts are run by dataset using the Immunarch package for uniform loading of TCR repertoires. The "stats_plots.R" script contain code for processing outputs into visualizations.

### Motif
This folder contains code for k-mer analyses (including frequency matrix generation, principal components analysis, and heatmaps), GLIPH2 analyses, and OLGA generation probability calculations. External tools for running GLIPH2 and OLGA can be found at http://50.255.35.37:8080/ and https://github.com/statbiophys/OLGA respectively; the R scripts are used for processing or downstream analysis. 

### Transcriptome
This folder contains code for analysis of single cell RNA sequencing datasets using Seurat. R scripts with the suffix "process.R" are used to process transcriptomes and merge them into a single Seurat object, which is then exported to analysis with the remaining R scripts on a high performance computing cluster. 

### ML
This folder contains code for training and evaluation of machine learning models. "ML_models.py" contains code that trains classical machine learning models and visualizes their ROC curves using the Plotly library. ML models were trained using k-mer representations of TCR repertoire data. Models were implemented with scikit-learn.


## Contact

Please see manuscript for further details. For any questions, contact jonathan [dot] park [at] yale [dot] edu.
This study was conducted by the Sidi Chen lab at Yale University (https://sidichenlab.org/).
