#  Selection and exploration of Highly Associated Variants (HAVs) 

## Required R / Bioconductor packages

* [topGO](https://bioconductor.org/packages/release/bioc/html/topGO.html)

* [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) 

* [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html) 


## Goals 

* Load and plot all the p-values from GLM and pGLS association tests and compares these agains p-values from simulations (Fig. 2A). 

* Select Highly Associated Variants (HAVs) and find near which genes they are

* Run Gene Ontology (GO) analysis

* Make a violin plot for association of the top SNP (Fig. 2D).
 

## Usage

All the analyses are documented in R script `GWAS_furtherSNPinvestigation.R` and the input files in this folder. 

