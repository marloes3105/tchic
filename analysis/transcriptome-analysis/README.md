# Transcriptome notebooks

## Overview of added notebooks:
1. Scanpy-scVelo-exampleNotebook:<br/>
This is one notebook combining Scanpy preprocessing, plotting and processing, followed by scVelo for velocity and pseudotime analysis included in the scVelo package. Files required: .loom output files produced by velocyto in the snakemake workflow.

2. QC_cellTyping:<br/>
This notebook does Scanpy preprocessing, clustering, plotting and defining cell types. The adata file is then saved and can be imported for specific analyses, e.g. trajectory inference or for merging chic data and plotting all measurements together. Files required: .loom output files produced by velocyto in the snakemake workflow.


3. Trajectories:<br/>
This R-script loads the necessary files pulled from the Scanpy adata file and runs monocle (standard settings) for trajectory inference/pseudotime analysis. To generate the files necessary for this R script you can use the following:
```
#adata_var
adata.var.to_csv("adata_var.csv")

#adata_obs
adata.obs.to_csv("adata_obs.csv")

#umap coordinates
pd.DataFrame(adata.obsm['X_umap']).to_csv("umap_coord.csv")

#count matrix - sparse!!
## import packages
import scipy
import scipy.sparse as sparse
import scipy.io as sio
import scipy.stats as stats
import numpy as np
## extract count table and make sure it is integers
counts_int = pd.DataFrame(adata.layers['matrix'].toarray())
counts_int = counts_int.astype(int)
## format colnames, rownames etc (not sure if this matters?) and transpose
counts_int.index = adata.obs.index
counts_int.columns = adata.var.index
counts_int = pd.DataFrame(counts_int).T
## make a scipy sparse matrix (!!) and write
counts_int_sparse = scipy.sparse.csr_matrix(counts_int)
sio.mmwrite("sparse_matrix_int.mtx",counts_int_sparse)
```


## Required packages
All notebooks are python-based and therefore require basic packages such as matplotlib, seaborn, pandas, and numpy.
All RNAseq analyses are based on Scanpy. For more information, see [here](https://scanpy.readthedocs.io/en/stable/). Tutorials, various ways of visualising the data and more can be found [here](https://scanpy.readthedocs.io/en/stable/tutorials.html). 
RNA velocity is computed using scVelo, for more information see their documentation. https://scvelo.readthedocs.io/ Especially the tutorials can be useful for [basic velocity analyis](https://scvelo.readthedocs.io/VelocityBasics.html) and [various ways of visualising and modelling](https://scvelo.readthedocs.io/DynamicalModeling.html).
Trajectory inference using PAGA is based on [this tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/paga-paul15.html). We also use wishbone, for more information see [here](https://scanpy.readthedocs.io/en/stable/external/scanpy.external.tl.wishbone.html) and [here](https://github.com/ManuSetty/wishbone).


