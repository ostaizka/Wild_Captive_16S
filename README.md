# Diversity and compositional changes in the gut microbiota of wild and captive vertebrates: a meta-analysis
This repository contains the bioinformatic resources related to the meta-analysis of gut microbiota features of captive and wild specimens

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5594740.svg)](https://doi.org/10.5281/zenodo.5594740)


## Raw data

Accession numbers to the raw sequencing data can be found in the column "Accession" of the file Data/metadata.filtered.csv.

## Supplementary code A
### Bioinformatic processing pipeline

1) Data input
2) Transfer files from repository to server
3) Get data properties
4) Check whether the file contains primers
5) Quality-filter and collapse
6) Remove primers
7) Convert to fasta
8) Add sample name
9) Calculate minimum cutoff value
10) Dereplicate
11) Denoise
12) Create zotu table
13) Assign taxonomy
14) Aggregate taxonomy

## Supplementary code B
### Diversity analysis pipeline

 1) Prepare working abundance-tables
 2) Across-species differences
 3) Diversity summaries
 4) Overall wild vs captive diversity meta-analyses
 5) Taxon-specific diversity meta-analyses
 6) Compositional differences between wild and captive animals (diversity partitioning)
 7) Compositional differences (analysis of variance)
 8) Distribution of the origin of the detected genera
 9) Pairwise dissimilarity correlation between wild and captive
 10) Nearest dataset
 11) Hierarchical clustering and topological differences
 12) NMDS plot

## File archive
### Documents used for and generated by the diversity analysis pipeline
- **Data**: input data for the diversity analyses
- **Results**: intermediate and final results of the diversity analyses
  - **Plots**: various plots generated in the diversity analyses
  - **RDS**: matrices and other large objects generated in the diversity analyses
- **Species**: read-abundance tables of individual samples.
- **Tables**: abundance tables of species.
