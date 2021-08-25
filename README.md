# RasperGade16S: Prediction and correction of 16S rRNA GCN with RasperGade in 16S-profiled bacterial communities
`RasperGade16S` is a specialized R package that predicts ancestral and hidden states while accounting for pulsed evolution and time-independent variation.
In addition to what `RasperGade` can do, `RasperGade16S` has some functions dedicated to 16S rRNA GCN:

1. addressing rate heterogeneity with a binary rate model

2. providing an established reference to predict 16S rRNA GCN 

3. correcting 16S rRNA GCN variation in analyses of bacterial composition with confidence provided

Detailed analyses are described in the preprint: `16S rRNA copy number prediction with confidence estimate and its implications in microbial diversity analyses`

doi: 

## System requirements

## Installation (Not available at the moment)
```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github(repo = "wu-lab-uva/RasperGade16S")
```
## Data format

## Reproduction of results in the preprint
To reproduce the results in the preprint, follow the instrucitons in
