## Peak calling using SEACR

1. Make sure you have SEACR installed (https://github.com/FredHutch/SEACR#readme)
Alternatively, use the following conda environment:
```
/hpc/hub_oudenaarden/Marloes/miniconda3/envs/cutruntools/
```

2. Make pseudobulks of each cell type 

3. Go to folder with pseudobulk bam files, copy files from this folder (sizes.genome and peak_calling.sh) and submit the bash script:
```
submit peak_calling.sh
```
Output files are called $sample.stringent.bed
