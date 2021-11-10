# RasperGade16S: Prediction and correction of 16S rRNA GCN
`RasperGade16S` is a specialized R package that predicts the 16S rRNA GCN of bacteria given their 16S rRNA gene sequence.

In addition to what `RasperGade` can do, `RasperGade16S` has some functions dedicated to 16S rRNA GCN:

1. addressing rate heterogeneity with a binary rate model

2. providing an established pipeline to predict 16S rRNA GCN 

3. correcting 16S rRNA GCN variation in analyses of bacterial composition with confidence provided

Detailed analyses are described in the preprint: `Accounting for 16S rRNA copy number prediction uncertainty and its implications in microbial diversity analyses`

doi: --

## System requirements
For full functionality of RasperGade16S, POSIX-compliant system (e.g., Linux or macOS but not Windows) is required.

The package has been tested on Ubuntu and macOS.

The following software is required:

HMMER 3 (available from http://hmmer.org/)

EPA-ng (available from  https://github.com/Pbdas/epa-ng)

For macOS, you can install HMMER3 and EPA-ng using brew (https://brew.sh)
    brew install hmmer
    brew install easel
    brew install epa-ng

The following packages are required:
`ape`,`castor`,`RasperGade`,`phyloseq`,`vegan`,`seqinr`,`treeio`,`microbiome`

## Installation
`RasperGade` and `RasperGade16S` can be installed using the following command in R
```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github(repo = "wu-lab-uva/RasperGade")
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
Other packages can be installed from CRAN
```
install.packages("ape")
install.packages("castor")
install.packages("vegan")
install.packages("seqinr")
```
## Data format
RasperGade16S can start from 16S rRNA sequences (e.g., representative sequences for OTUs)

## Demo 
Once installed, a small demo can be run in R to check if RasperGade16S operates properly:
```
library(RasperGade16S)
pred.GCN = predict_16SGCN_from_sequences()
print(pred.GCN$tab)
```
If everything works well, the following output should be expected
```
              label x      probs
0.1 GY203941.1.1493 4 0.99999796
0.2 GY324971.1.1500 9 0.82278712
0.3 JQ765433.1.1505 8 0.04355867
0.4 JQ765578.1.1444 6 0.83112359
0.5 JQ766308.1.1248 2 0.99966469
```
To run on your own sequence, supply the path to the FASTA file via the seqs= argument
```
pred.GCN = predict_16SGCN_from_sequences(seqs=PATH_TO_YOUR_SEQUENCE)
```
