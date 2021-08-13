# TChIC
to share all things tchic

## Overview:
1. **workflows:** <br>
This folder contains 3 sub-folders: 1. snakemake-workflow, 2. cell-selection, 3. snakemake-workflow. <br>
The first snakemake processes both chic and transcriptome data from fastq until count tables (or for transcriptome: .loom files with RNA velocity information). After this, the transcriptome data can be analysed locally using the notebooks in the analysis/ folder. The next two steps are specific for chic data. Step 2 uses information from the transcriptome to define similar (neighbour) cells and define bad cells in the chic data. Step 3 filters the chic data and removes bad cells, produces TSS and gene body count tables and creates pseudobulk bigwig files using your favourite clusters (e.g. cell types, time points, replicates).


2. **analysis:** <br>
a. *transcriptome-notebooks:* <br/>
Here you can find the jupyter notebook for analysis of the transcriptome. For this we use Scanpy and scVelo. Here you can also make a csv file with cell names x cell types for splitting of the chic bam file into pseudobulk files, and a csv file with umap coordinates for filtering out bad cells in the chic data. <br>
b. *chic-analysis:* <br/>
Here you can find everything for analysing the chic data. Currently added:
Test notebooks to plot chic data on top of the umap and cluster chic data, notebook to normalise chic data using percentile normalisation. Normalised data is used as input for the merged chic+transcriptome analysis. <br>
c. *merged-chic-trans:* <br/>
Here you can find analyses that merge transcriptome and chic data. Currently contains a notebook on how to merge chic dataframes into the transcriptome adata file. This adata file is the result of the "QC_cellTyping" notebook found in transcriptome-notebooks. Chic input files are raw and normalised TSS tables (raw data:output from workflow 3, normalised data: output from chic-analysis normalisation notebook). Merging these chic TSS tables with the scanpy adata file allows us to plot spliced, unspliced, and chic coverage with the same plots supported by scanpy (e.g. coverage on the umap, dot plots, violin plots and so on). After merging together chic and transcriptome data, you can use the notebooks for plotting, finding marker genes and making smooth line plots.

3. **data:**<br/>
Here you can find relevant data files, e.g. cell type cluster csv files to split bam files into pseudobulks.

<br>


## Necessary dependencies/packages:
For running the snakemake workflow on the cluster: <br/>
Please install a conda or python virtual environment and install all necessary dependencies in this environment. Note that this needs to be activated before running the snakemake. The following dependencies are necessary: [SingleCellMultiOmics](https://github.com/BuysDB/SingleCellMultiOmics) package, cutadapt, BWA or bowtie, samtools, STAR, velocyto.

For local analysis/running code in jupyter notebooks: <br/>
All notebooks are python-based and therefore require basic packages such as matplotlib, seaborn, pandas, and numpy. <br/>
All RNAseq analyses are based on Scanpy. For more information, see [here](https://scanpy.readthedocs.io/en/stable/). Tutorials, various ways of visualising the data and more can be found [here](https://scanpy.readthedocs.io/en/stable/tutorials.html). <br/>
RNA velocity is computed using scVelo, for more information see their [documentation](https://scvelo.readthedocs.io/). Especially the tutorials can be useful for [basic velocity analyis](https://scvelo.readthedocs.io/VelocityBasics.html) and [various ways of visualising and modelling](https://scvelo.readthedocs.io/DynamicalModeling.html).<br/>
Trajectory inference using PAGA is based on [this tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/paga-paul15.html). We also use wishbone, for more information see [here](https://scanpy.readthedocs.io/en/stable/external/scanpy.external.tl.wishbone.html) and [here](https://github.com/ManuSetty/wishbone).




