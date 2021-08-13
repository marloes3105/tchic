# Transcriptome notebooks

## Overview of added notebooks:
1. Scanpy-scVelo-exampleNotebook:<br/>
This is one notebook combining Scanpy preprocessing, plotting and processing, followed by scVelo for velocity and pseudotime analysis included in the scVelo package. Files required: .loom output files produced by velocyto in the snakemake workflow.

2. QC_cellTyping:<br/>
This notebook does Scanpy preprocessing, clustering, plotting and defining cell types. The adata file is then saved and can be imported for specific analyses, e.g. trajectory inference or for merging chic data and plotting all measurements together. Files required: .loom output files produced by velocyto in the snakemake workflow.


3. Trajectories:<br/>
This notebook takes the adata file produced by the QC_cellTyping notebook as input. Here we can calculate trajectories and pseudotime using PAGA and Wishbone - these are external Scanpy packages, meaning they need to be installed separately but are integrated with scanpy datafiles (such as adata). We can also compute RNA velocity using scVelo and hopefully find some interesting genes in terms of velocity, or look at differences across the calculated pseudotime in terms of spliced/unspliced transcripts.


## Required packages
All notebooks are python-based and therefore require basic packages such as matplotlib, seaborn, pandas, and numpy.
All RNAseq analyses are based on Scanpy. For more information, see [here](https://scanpy.readthedocs.io/en/stable/). Tutorials, various ways of visualising the data and more can be found [here](https://scanpy.readthedocs.io/en/stable/tutorials.html). 
RNA velocity is computed using scVelo, for more information see their documentation. https://scvelo.readthedocs.io/ Especially the tutorials can be useful for [basic velocity analyis](https://scvelo.readthedocs.io/VelocityBasics.html) and [various ways of visualising and modelling](https://scvelo.readthedocs.io/DynamicalModeling.html).
Trajectory inference using PAGA is based on [this tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/paga-paul15.html). We also use wishbone, for more information see [here](https://scanpy.readthedocs.io/en/stable/external/scanpy.external.tl.wishbone.html) and [here](https://github.com/ManuSetty/wishbone).


