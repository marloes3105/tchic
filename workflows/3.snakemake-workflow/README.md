# Processing of chic data

## necessary inputs:
- merged bam file to be processed
- output of cell selection script (step 2)
- csv file with cell names x information of clusters to split on - this can be cell types, days, replicates etc
- TSS and gene body bed files

<br> 

## summary
This workflow takes one merged bam, filters out bad cells and from the 'good cells' fraction it makes chic coverage tables of transcription start sites (TSS) and gene bodies. It also splits the bam file into your preferred pseudobulk clusters (e.g. cell types) and makes bigwig files.

<br>

## step-by-step overview of workflow
1. Remove bad cells. <br>
Using the output from step 2, we can now remove cells that are classified as 'bad'. We simply split the merged bam file into 2: one with bad cells, one with good cells. 

Example code: <br>
```
bamExtractSamples.py merged_tagged.bam k4me3_QCselector_round2.csv -o processed/merged_tagged_.bam

```

<br>

2. Split bam file into clusters <br>
The filtered bam file is now split into your preferred clusters, based on your transcriptome data. These can be your defined cell types, days, replicates etc. For each cluster, it outputs an indexed bam file.

Example code: <br>
```
bamExtractSamples.py merged_tagged_good.bam celltypes_updated.csv -o processed/bigwigs/subset_.bam
```

<br>

3. Make bigwigs out of split bam files <br>
Each 'pseudobulk' bam file is now transformed into a bigwig file. These can easily be loaded into IGV to look at cluster-specific coverage of your favourite genomic locations. Parameters can be changed in the config file. At the end of this step, the split bam files are removed again.

Example code: <br>
```
bamCoverage --ignoreDuplicates --binSize 200 --smoothLength 1000 --minMappingQuality 50 --centerReads --normalizeUsing CPM -b subset.bam -o subset.bw
```

<br>

4. Make feature tables <br>
The filtered merged bam file is also used to make coverage tables of your favourite feature. Currently, we make these tables for gene bodies and transcription start sites (TSS). These csv tables can be used in the downstream analysis by loading them into the transcriptome adata file (see jupyter notebooks for merged chic and transcriptome analysis). 

Example code for TSS: <br>
```
bamToCountTable.py --filterXA -minMQ 50 merged_tagged_good.bam -o count_table_TSS.csv -sampleTags SM -joinedFeatureTags reference_name --dedup --r1only -bedfile /hpc/hub_oudenaarden/Peter/annotations/TSS_10kbwindow_mouse_ensmble97.bed
```

