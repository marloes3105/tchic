# Analysis of both transcriptome and chic

To merge the two data types, please first use the "AddChICToAdataLayers" notebook. This notebook requires (at least) 2 input files: 
(1) the adata object produced by the QC_cellTyping notebook (see tchic/transcriptome-notebooks for more information) and 
(2) TSS coverage count table(s) of your favourite chic data. <br/>
To add normalised data, please see chic-analysis/single-cell/20210305_normaliseChICData.ipynb

This notebook formats the chic count tables to make sure it has the same dimensions and ordering as compared to the transcriptome data, and merges these count tables with the adata file in a different layer. 
Scanpy supports plotting of different layers by adding the following argument: layer = 'x', where 'x' is 'spliced', 'unspliced' (for spliced/unspliced counts), 
or one of the names of the added chic layers (default: 'k4_raw' and 'k27_raw', for raw counts of TSS tables from k4me3 and k27me3, respectively). 


Various plotting examples for all layers (spliced, unspliced, chic histone presence) can be found in the notebook "plotsAdataWithChICLayers". This notebook requires the adata file produced by the "AddChICToAdataLayers" notebook as input file.


If you systematically want to define marker genes (top differentially expressed genes) between e.g. cell types or leiden clusters, you can use the "trans+ChIC_findMarkerGenes" notebook. You can define marker genes based on the transcriptome (standard), but by defining the layer you can also extend this to specifically spliced or unspliced counts, or define marker genes based on chic data. Simply specify which data type you want to use by adding the argument "layer = x" (similar to how you can plot these different data types). 
<br/>
