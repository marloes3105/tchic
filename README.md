# tchic
to share all things tchic

## Overview:
1. Snakemake-workflow
Here you can find the Snakefile and config.json for running the first steps on the cluster, from fastq files to bam and csv files (and a loom file for transcriptome analysis).


2. transcriptome-notebook
Here you can find the jupyter notebook for analysis of the transcriptome. For this we use Scanpy and scVelo. Here you can also make a csv file with cell names x cell types for splitting of the chic bam file into pseudobulk files.


3. chic-analysis
Here you can find everything for analysing the chic data. Currently added:
- peak-calling: to run SEACR and call peaks. We need to check if this works at all and gives reliable results.
- single-cell: currently containing a notebook to plot chic coverage of specific genome areas on top of the transcriptome umap. In the example, we plot TSS coverage of k4me3 and k27me3 for specific genes.

to be added:
- pseudobulk analysis
- how to get TSS/gene body coverage csv count tables.
