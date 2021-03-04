# tchic
to share all things tchic

## Overview:
1. Snakemake-workflow:<br/>
Here you can find the Snakefile and config.json for running the first steps on the cluster, from fastq files to bam and csv files (and a loom file for transcriptome analysis).


2. transcriptome-notebooks:<br/>
Here you can find the jupyter notebook for analysis of the transcriptome. For this we use Scanpy and scVelo. Here you can also make a csv file with cell names x cell types for splitting of the chic bam file into pseudobulk files.


3. transcriptome-data: <br/>
Here you can find relevant data files, e.g. cell type cluster csv files to split bam files into pseudobulks.


4. chic-analysis: <br/>
Here you can find everything for analysing the chic data. Currently added:
- peak-calling: to run SEACR and call peaks. We need to check if this works at all and gives reliable results.
- single-cell: currently containing: (1) a notebook to plot chic coverage of specific genome areas on top of the transcriptome umap. In the example, we plot TSS coverage of k4me3 and k27me3 for specific genes. (2) a notebook to cluster chic data using scanpy. this is a test and still needs some work, for example clustering is dominated by read count per cell. Could be a good check to see if there is any type of structure in the data (if you get separate clusters of cells).


5. merged-chic-trans: <br/>
Here you can find analyses that merge transcriptome and chic data. Currently contains a notebook on how to merge chic dataframes into the transcriptome adata file. This adata file is the result of the "QC_cellTyping" notebook found in transcriptome-notebooks. Chic input files are TSS tables. Merging these chic TSS tables with the scanpy adata file allows us to plot spliced, unspliced, and chic coverage with the same plots supported by scanpy (e.g. coverage on the umap, dot plots, violin plots and so on). To be added: (1) a notebook with plotting examples, (2) a way to normalise chic TSS tables and merge these too.


## to be added:<br/>
- pseudobulk analysis
- how to get TSS/gene body coverage csv count tables.


## Necessary dependencies/packages:
For running the snakemake workflow on the cluster: <br/>
Please install a conda or python virtual environment and install all necessary dependencies in this environment. Note that this needs to be activated before running the snakemake. The following dependencies are necessary: [SingleCellMultiOmics](https://github.com/BuysDB/SingleCellMultiOmics) package, cutadapt, BWA or bowtie, samtools, STAR, velocyto.

For local analysis/running code in jupyter notebooks: <br/>
All notebooks are python-based and therefore require basic packages such as matplotlib, seaborn, pandas, and numpy. <br/>
All RNAseq analyses are based on Scanpy. For more information, see [here](https://scanpy.readthedocs.io/en/stable/). Tutorials, various ways of visualising the data and more can be found [here](https://scanpy.readthedocs.io/en/stable/tutorials.html). <br/>
RNA velocity is computed using scVelo, for more information see their [documentation](https://scvelo.readthedocs.io/). Especially the tutorials can be useful for [basic velocity analyis](https://scvelo.readthedocs.io/VelocityBasics.html) and [various ways of visualising and modelling](https://scvelo.readthedocs.io/DynamicalModeling.html).<br/>
Trajectory inference using PAGA is based on [this tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/paga-paul15.html). We also use wishbone, for more information see [here](https://scanpy.readthedocs.io/en/stable/external/scanpy.external.tl.wishbone.html) and [here](https://github.com/ManuSetty/wishbone).




