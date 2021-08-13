# Processing of TChIC data

This folder contains 3 sub-folders:
1. snakemake-workflow
2. cell-selection
3. snakemake-workflow

<br>

## summary
The first snakemake processes both chic and transcriptome data from fastq until count tables (or for transcriptome: .loom files with RNA velocity information). After this, the transcriptome data can be analysed locally using the notebooks in the analysis/ folder. The next two steps are specific for chic data. Step 2 uses information from the transcriptome to define similar (neighbour) cells and define bad cells in the chic data. Step 3 filters the chic data and removes bad cells, produces TSS and gene body count tables and creates pseudobulk bigwig files using your favourite clusters (e.g. cell types, time points, replicates).
