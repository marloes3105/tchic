# Processing of TChIC data

This folder contains 3 sub-folders:
1. snakemake-workflow
2. cell-selection
3. snakemake-workflow

<br>

## summary
(1) The first snakemake processes both chic and transcriptome data from fastq until count tables (or for transcriptome: .loom files with RNA velocity information). After this, the transcriptome data can be analysed locally using the notebooks in the next step. 

(2) Step 2 splits the ChIC and transcriptome read-outs: the transcriptome .loom files (of all cells together, not split per histone mark) can be read in using the script '1_TCHC_ReadInPlates'. Then, basic QC is performed on the transcriptome fraction using scanpy, and principle components are calculated and exported in notebook 2a. These will be used to perform QC on the chic fraction.
In notebook 2b, QC is performed on the ChIC modality. This notebook uses information from the transcriptome to define similar (neighbour) cells and define bad cells in the chic data. This should be performed on each ChIC modification separately.
In notebook 2c, we go back to the transcriptome modality, filter out the QCfail cells and generate the UMAP.

(3) Step 3 filters the ChIC data and removes bad cells using the output of step 2, produces TSS and gene body count tables and creates pseudobulk bigwig files using your favourite clusters (e.g. cell types, time points, replicates).
Here we also process the transcriptome further, perform cell typing, calculate pseudotime etc.
