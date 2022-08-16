[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6997675.svg)](https://doi.org/10.5281/zenodo.6997675)  [![DOI](https://img.shields.io/badge/bioR%CF%87iv%20DOI-10.1101%2F2022.04.26.489598%20-BC2635)](https://doi.org/10.1101/2022.04.26.489598)

# A critique of the hypothesis that CA repeats are primary targets of neuronal MeCP2 

This repository contains source code and files required to replicate the findings from our manuscript.

> [Chhatbar, Kashyap, John Connelly, Shaun Webb, Skirmantas Kriaucionis, and Adrian Bird. "A critique of the hypothesis that CA repeats are primary targets of neuronal MeCP2." bioRxiv (2022).](https://doi.org/10.1101/2022.04.26.489598)

Contrary to the claims of a recent publication,

> [Ibrahim, Abdulkhaleg, Christophe Papin, Kareem Mohideen-Abdul, Stéphanie Le Gras, Isabelle Stoll, Christian Bronner, Stefan Dimitrov, Bruno P. Klaholz, and Ali Hamiche. "MeCP2 is a microsatellite binding protein that protects CA repeats from nucleosome invasion." Science 372, no. 6549 (2021): eabd5581.](https://doi.org/10.1126/science.abd5581)

which proposes that arrays of tandemly repeated CA containing either methylated or hydroxymethylated cytosine are the primary targets for MeCP2 binding and function. We investigated these predictions using a range of published data sets and found **no support for the hypothesis that CA repeats are key mediators of MeCP2 function**.

Before starting the analysis, fulfill the [installation requirements](prerequisities.md). They should install the required software necessary for reproducing the analysis.

The repository is divided into multiple directories containing step-by-step information formatted as jupyter notebooks. The execution order of the notebooks is as follows:

1. CA repeats in the genome
    1. [Extract CA repeat loci in mammalian genomes](methods/01_CA_repeat_loci.ipynb)    
    2. [Count number of CAC occurences in mammalian genomes](methods/02_CAC_counts_genome.ipynb)
    
2. Quantification of DNA methylation in mammalian brain
    1. [Download processed Whole genome bisulfite sequencing data sets](wgbs/01_download.ipynb)
    2. [Quantify DNA methylation for different cytosine contexts in mammalian brain tissues](wgbs/02_postprocess.ipynb)

3. Quantification of MeCP2 enrichment in mammalian brain
    1. [Download and process raw ChIPseq data sets](chipseq/01_download_process.ipynb)
    2. [Divide the genome into 1kb windows](chipseq/02_windows.ipynb)

4. Comparative analysis of DNA methylation and MeCP2 enrichment in mammalian brain
    1. [MeCP2 enrichment from ChIP-seq vs DNA methylation of genomic windows in mammalian brain](chipseq/03_analysis.ipynb)

5. mRNA abundance quantification in MeCP2 wildtype (WT) and MeCP2 Knockout (KO) brain
    1. [Download and process raw RNAseq data sets](rnaseq/01_download_process.ipynb)

6. Comparative analysis of DNA methylation and transcriptional regulation in the absence of MeCP2
    1. [Differential gene expression vs DNA methylation of gene bodies in mammalian brain](rnaseq/02_analysis.ipynb)


The notebooks to generate all the figures can be found in this [directory](figures).
