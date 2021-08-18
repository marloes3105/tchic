# Analysis of both transcriptome and chic

To merge the two data types, please first use the "AddChICToAdataLayers" notebook. This notebook requires (at least) 2 input files: 
(1) the adata object produced by the QC_cellTyping notebook (see tchic/transcriptome-notebooks for more information) and 
(2) TSS coverage count table(s) of your favourite chic data. <br/>
To add normalised data, please see chic-analysis/single-cell/20210305_normaliseChICData.ipynb

This notebook formats the chic count tables to make sure it has the same dimensions and ordering as compared to the transcriptome data, and merges these count tables with the adata file in a different layer. 
Scanpy supports plotting of different layers by adding the following argument: layer = 'x', where 'x' is 'spliced', 'unspliced' (for spliced/unspliced counts), 
or one of the names of the added chic layers (default: 'k4_raw' and 'k27_raw', for raw counts of TSS tables from k4me3 and k27me3, respectively). 

