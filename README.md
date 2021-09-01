# RasperGade16S: Prediction and correction of 16S rRNA GCN
`RasperGade16S` is a specialized R package that predicts the 16S rRNA GCN of bacteria given their 16S rRNA gene sequence.

In addition to what `RasperGade` can do, `RasperGade16S` has some functions dedicated to 16S rRNA GCN:

1. addressing rate heterogeneity with a binary rate model

2. providing an established pipeline to predict 16S rRNA GCN 

3. correcting 16S rRNA GCN variation in analyses of bacterial composition with confidence provided

Detailed analyses are described in the preprint: `Accounting for 16S rRNA copy number prediction uncertainty and its implications in microbial diversity analyses`

doi: --

## System requirements
For full functionality of RasperGade16S, POSIX-compliant system (e.g., Linux or macOS) is required.

The package has been tested on Ubuntu and macOS.

The following software is required:

HMMER 3 (available from http://hmmer.org/)

EPA-ng (available from  https://github.com/Pbdas/epa-ng)

The following packages are required:
`ape`,`castor`,`RasperGade`,`phyloseq`,`vegan`,`seqinr`,`treeio`,`microbiome`

## Installation
RasperGade16S can be installed using the following command in R
```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github(repo = "wu-lab-uva/RasperGade16S")
```
`treeio`,`phyloseq`,`microbiome` can be installed from Bioconductor
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("phyloseq")
BiocManager::install("treeio")
BiocManager::install("microbiome")
```
Other packages can be installed from CRAN.
## Data format
RasperGade16S can start from 16S rRNA sequences (e.g., representative sequences for OTUs)
## Reproduction of results in the preprint
To reproduce the results in the preprint, follow the instrucitons in
