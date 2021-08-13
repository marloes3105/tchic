# Overview

This snakemake takes TChIC fastq files, splits them into chic and transcriptome libraries, does trimming, mapping tagging and makes count tables (or for transcriptome: loom files). ChIC output files are stored in the folder `processed_chic/`, transcriptome output fules are stored in the folder `processed_transcriptome/`. The most convenient way of running this snakemake is pooling all fastq files per chic mark into one folder and run the worflow for each mark separately. This way, the script generates merged files per mark. <br>

Chic data processing is based on the existing [Chic pipeline](https://github.com/BuysDB/SingleCellMultiOmics/wiki/scCHIC-data-processing) <br>
Transcriptome data processing is based on the existing [vasa pipeline](https://github.com/BuysDB/SingleCellMultiOmics/tree/master/singlecellmultiomics/snakemake_workflows/vasa) <br>

All steps are described in more detail below. <br>

### to run this snakemake:
- Edit the config file (config.json), make sure to set the right paths to the reference files and check if the parameters are correct.
- Submit the workflow to the cluster. Important: when running from a conda environment, do not forget to use --use-conda! <p>
  Example submission using submission.py and a conda environment: 
  ```python
  submission.py "source activate conda;snakemake --use-conda --cluster slurm_wrapper.py --jobs 20 --restart-times 3" -y -time 50 -m 2 -sched slurm -N tchic
  ```

<br>


## ChIC pipeline step-by-step
  
This pipeline contains the following steps:
1. Demultiplexing.
2. Adapter trimming using cutadapt.
3. Mapping using BWA or bowtie2.
4. Tagging the bam file with relevant information.
5. Filtering the bam file to throw out any low-quality reads and contamination.
6. Removing VASA bleedthrough.
7. Generating library quality statistics plots.
8. Generating binned count tables of filtered data.
9. Merging all chic libraries into one merged_tagged.bam, calculating library statistics and making a binned count table of the merged data. 

<br> 
  
**1. Demultiplexing** <br>
Demultiplexing adds UMI, cell, sample and sequencing index information to the header of the FastQ file. This information is later used to figure out which reads belong to the same molecule and to what cell every molecule belongs. The demultiplexer uses barcodes available in the [barcodes folder](https://github.com/BuysDB/SingleCellMultiOmics/tree/master/singlecellmultiomics/modularDemultiplexer/barcodes). The right set of barcodes is automatically detected if `-use` is not supplied.
  
Example code:
```
demux.py -merge _ *fastq.gz -use scCHIC384C8U3l -hd 0 -o processed_chic --y 
```
  
  <br>
  
**2. Adapter trimming.** <br>
Removes residual primers, this improves mapping performance.
  
Example code:
```
cutadapt -o trimmed.R1.fastq.gz -p trimmed.R2.fastq.gz demultiplexedR1.fastq.gz demultiplexedR2.fastq.gz -m 3 -a "IlluminaSmallAdapterConcatBait=GGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTT" -a "IlluminaIndexAdapter=GGAATTCTCGGGTGCCAAGGAACTCCAGTCACN{6}ATCTCGTATGCCGTCTTCTGCTTG"  -A "IlluminaPairedEndPCRPrimer2.0=AGATCGGAAGAGCGN{6}CAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG;min_overlap=5" -A "universalPrimer=GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT;min_overlap=5" -a  "IlluminaGEX=TTTTTAATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGACGATC;min_overlap=5" -a "IlluminaMultiplexingPCRPrimer=GGAACTCCAGTCACN{6}TCTCGTATGCCGTCTTCTGCTTG;min_overlap=5" -A "Aseq=TGGCACCCGAGAATTCCA" -a "Aseq=TGGCACCCGAGAATTCCA"  -a "illuminaSmallRNAAdapter=TCGTATGCCGTCTTCTGCTTGT"
```
  <br>
  
**3. Mapping** <br>
The trimmed fastq files can now be mapped using a mapper of your choice. In this example we will use BWA. We have tested mapping with bwa and bowtie2 and although both are possible, bowtie2 does not perform as well as bwa mapping this data. Importantly, when we map with bwa, we add the `-I` flag, defining the expected insert size of our amplicon. Since we do not expect fragments larger than 1000bp, we set the `-I` flag at 1000. After mapping, we sort and index the bam file.

Example code:
```
bwa mem -M -I 1000 -t 8 /hpc/hub_oudenaarden/group_references/ensembl/97/mus_musculus/primary_assembly_97_129B6Masked_ERCC92.fa trimmed.R1.fastq.gz trimmed.R2.fastq.gz
```
  <br>
  
**4. Tagging** <br>
The mapped bam file is tagged using bamtagmultiome.py `-method chict`. This places all information that was placed in the read name for alignment back into the bamfile as tags. You can also enable or disable the use of a blacklist. It is recommended to use a blacklist, this speeds up the script and removes reads mapped to uninformative genomic regions.
  
Example code:
```
bamtagmultiome.py --multiprocess -tagthreads 8 -blacklist /hpc/hub_oudenaarden/bdebarbanson/ref/500_plus_coverage_spikes_blacklist_CHiC_mm10.bed -introns /hpc/hub_oudenaarden/group_references/ensembl/97/mus_musculus/introns.gtf.gz -exons /hpc/hub_oudenaarden/group_references/ensembl/97/mus_musculus/exons.gtf.gz -method chict sorted.bam -o tagged.bam
```
  <br>
  
**5. Filtering** <br>
In this step we filter out low-quality read pairs. Tag `fS` indicates the fragment size, and we use this tag to remove read pairs that span a region of >1000 bp and reads that are not properly paired. Tag `MQ` indicates the mapping quality and we use this tag to remove read pairs with a mapping quality below 50. These values can be changed in the config file by changing the `insertSize` and `counting_min_mq` parameters, resp.
  
Example code:
```
bamFilter.py tagged.bam -o tagged_filtered.bam 'r.has_tag("fS") and (r.get_tag("fS") < 1000) and r.has_tag("MQ") and (r.get_tag("MQ") >= 50 )' 

```
  <br>
  
**6. Removing VASA bleedthrough.** <br>
We still might have a small fraction of transcriptome reads as bleed-through in the chic fraction. There's several ways of removing these reads, and this pipeline includes 2 options which you can choose by setting the `cleanup` parameter in the config file. If you choose `simple`, the cleanup is removed based on the ligation motif only. Chic reads have a ligation motif of TA or TT. However, ~40% of bleed-through transcriptome reads have this ligation motif too. The remainder of the transcriptome bleed-through mainly contains ligation motifs CA, CT, CC, TC, or TG. We use this information to select reads with TA or TT ligation motif, and remove reads with transcript-specific ligation motifs.
  
Example code:
```
bamFilter.py tagged_filtered.bam -o chic_classified.bam 'r.has_tag("lh") and (r.get_tag("lh") == "TA" or r.get_tag("lh") == "TT") or (r.get_tag("lh") != "CT" and r.get_tag("lh") != "CC" and r.get_tag("lh") != "CA" and r.get_tag("lh") != "TC" and r.get_tag("lh") != "TG")' '''
      
```

If you set the `cleanup` option to `training`, a random forest classifier is used instead. This classifier needs a chic-only (`chic_training_file`) and transcriptome-only (`vasa_training_file`) tagged bam file. It uses this to classify which features define a transcript read and which features define a chic read, and splits the bam file based on this classification. By default, it uses all features for this. You can however also pre-select which features to use by passing a list to the `-selection` option. 
  
Example code:
```
python /hpc/hub_oudenaarden/Marloes/bin/tchic/chic_transcriptome_split_updated.py -chic_bam /hpc/hub_oudenaarden/Marloes/bin/tchic/train_data/PZ-TChIC-K562-SC-H3K27me3-C/tagged.bam -vasa_bam /hpc/hub_oudenaarden/Marloes/bin/tchic/train_data/PZ-TChIC-K562-SC-H3K27me3-T/tagged.bam -tchic_bam tagged_filtered.bam -o processed_chic/{library}/ -selection TA,TC,TG,TT,CT,CC,CA
```
  
  <br>
  
**7. Generating library quality statistics plots.** <br>
Read the [Library-statistics-plots](https://github.com/BuysDB/SingleCellMultiOmics/wiki/Library-statistics-plots) page for information about library quality statistics.  
  
Example code:
```
libraryStatistics.py processed_chic/{library} -tagged_bam chic_classified.bam
```
  <br>
  
**8. Generating binned count tables.** <br>
The mapped and tagged bam file is converted to a count table using `bamToCountTable.py`. The minMQ parameter defines what the lowest mapping quality for a molecule should be to be counted. The `--filterXA` filters reads which map to multiple locations (excludes alternative loci). The `--dedup` flag makes sure every molecule is counted only once, if this flag is not specified every fragment in the BAM file with a `DS` tag (Which is a valid scCHIC site) will be counted. It takes `-joinedFeatureTags` as rows and `-sampleTags` as columns for this count table. For columns, we use the `SM` tag (samplename). For rows, we use the `reference_name` tag.  
  
Example code:
```
bamToCountTable.py -bin 50_000 -minMQ 50 --noNames chic_classified.bam -sampleTags SM -joinedFeatureTags reference_name -binTag DS -o count_table.csv --dedup --r1only

```
  <br>
  
**9. Merging all chic libraries.** <br>
Here we merge all chic libraries into one bam file. We also make a binned count table and generate library statistics plots of the merged file (as described above). This is convenient for downstream processing. 
  
Example code:
```
samtools merge -c merged_tagged.bam {library}/chic_classified.bam ; samtools index merged_tagged.bam"
```
  <br>
  
  



## Transcriptome pipeline step-by-step
  
This pipeline contains the following steps:
1. Demultiplexing.
2. Trimming away polyA tails.
3. Adapter trimming using cutadapt.
4. Mapping using STAR.
5. Tagging the bam file with relevant information.
6. Generating library quality statistics plots.  
7. Reformatting bam file headers for Velocyto.
8. Deduplicating.
9. Running Velocyto for RNA velocity.
10. Making exonic and intronic count tables.
  
<br> 
  
**1. Demultiplexing** <br>
Demultiplexing adds UMI, cell, sample and sequencing index information to the header of the FastQ file. This information is later used to figure out which reads belong to the same molecule and to what cell every molecule belongs. The demultiplexer uses barcodes available in the [barcodes folder](https://github.com/BuysDB/SingleCellMultiOmics/tree/master/singlecellmultiomics/modularDemultiplexer/barcodes). The right set of barcodes is automatically detected if `-use` is not supplied.
    
Example code:
```
demux.py -merge _ *fastq.gz -use CS2C8U6NH -hd 0 -o processed_chic --y 
```
  <br>
  
**2. PolyA trimming.** <br>
Vasa reads still contain polyA stretches. Removing these improves mapping performance. We only do this for R2 because we'll be doing single-end mapping anyway.
  
Example code:
```
trim_vasa.py demultiplexedR2.fastq.gz poly_trimmed.R2.fastq.gz -min_read_len 20 
```
  <br>
  
**3. Adapter trimming.** <br>
Removes residual primers, this improves mapping performance. We only do this for R2 because we'll be doing single-end mapping anyway.
  
Example code:
```
cutadapt -o trimmed.R2.fastq.gz poly_trimmed.R2.fastq.gz -m 3 -a "IlluminaSmallAdapterConcatBait=GGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTT" -a "IlluminaIndexAdapter=GGAATTCTCGGGTGCCAAGGAACTCCAGTCACN{6}ATCTCGTATGCCGTCTTCTGCTTG"  -A "IlluminaPairedEndPCRPrimer2.0=AGATCGGAAGAGCGN{6}CAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG;min_overlap=5" -A "universalPrimer=GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT;min_overlap=5" -a  "IlluminaGEX=TTTTTAATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGACGATC;min_overlap=5" -a "IlluminaMultiplexingPCRPrimer=GGAACTCCAGTCACN{6}TCTCGTATGCCGTCTTCTGCTTG;min_overlap=5" -A "Aseq=TGGCACCCGAGAATTCCA" -a "Aseq=TGGCACCCGAGAATTCCA"  -a "illuminaSmallRNAAdapter=TCGTATGCCGTCTTCTGCTTGT"
```
  
<br>
  
**4. Mapping.** <br>
Transcriptome reads are best mapped single-end using STAR.
  
Example code:
```
STAR --runThreadN 8 --readFilesCommand zcat  --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --outSAMmultNmax 10 --outFilterMultimapNmax 10 --genomeDir /hpc/hub_oudenaarden/group_references/ensembl/97/mus_musculus/129B6masked_ribo_star_index_99" --outSAMattributes All --readFilesIn trimmed.R2.fastq.gz --outSAMunmapped Within --outFileNamePrefix STAR_mapped_R2

```
  <br>
  
**5. Tagging.** <br>
The mapped bam file is tagged using bamtagmultiome.py `-method vasa`. This places all information that was placed in the read name for alignment back into the bamfile as tags. 
  
Example code:
```
bamtagmultiome.py --multiprocess -tagthreads 8 -introns /hpc/hub_oudenaarden/group_references/ensembl/97/mus_musculus/introns.gtf.gz -exons /hpc/hub_oudenaarden/group_references/ensembl/97/mus_musculus/exons.gtf.gz -method vasa STAR_mapped_R2Aligned.sortedByCoord.out.bam -o tagged.bam

```

  <br>
  
**6. Generating library quality statistics plots.** <br>
Read the [Library-statistics-plots](https://github.com/BuysDB/SingleCellMultiOmics/wiki/Library-statistics-plots) page for information about library quality statistics.  
Example code:
```
libraryStatistics.py processed_transcriptome/{library} -tagged_bam tagged.bam
```
  
  <br>
  
**7. Reformatting bam file headers for Velocyto.** <br>
This is necessary because otherwise Velocyto crashes. 
  
Example code:
```
scmoConvert.py tagged.bam tagged_converted.bam -fmt cellranger 

```
  
  <br>
  
**8. Deduplicating.** <br>
Velocyto does also support deduplicating, but is very agressive about it. This is why we first deduplicate ourselves and then disable velocyto's deduplication (see step 9). 
  
Example code:
```
samtools view -F 1024 tagged_converted.bam -bS -@32 > tagged_samtools-dedup.bam; samtools index tagged_samtools-dedup.bam
```
  
  <br>
  
**9. Running Velocyto.** <br>
This runs velocyto and writes the output .loom file to a folder called `rna_counts`. The option `-U` disables umi deduplication: "If this flag is used the data is assumed UMI-less and reads are counted instead of molecules". 
  
Example code:
```
velocyto run -e samplename_x -@ 8 -o rna_counts tagged_samtools-dedup.bam /hpc/hub_oudenaarden/group_references/ensembl/97/mus_musculus/annotations.gtf -U "
```
  
  <br>
  
**10. Making count tables.** <br>
 The mapped and tagged bam file is converted to a count table using bamToCountTable.py. It takes `-joinedFeatureTags` as columns and `-sampleTags` as rows for this count table. For columns, we use the `SM` tag (samplename). For rows, we use the `GN` tag for exonic counts and the `IN` tag for the intronic counts
  
  
  
  
Example code exonic count table:
```
bamToCountTable.py tagged.bam -sampleTags SM -joinedFeatureTags reference_name,GN -o count_table.csv --dedup --r1only 
```
  
Example code intronic count table:
```
bamToCountTable.py tagged.bam -sampleTags SM -joinedFeatureTags reference_name,IN -o count_table.csv --dedup --r1only
```
  
  
  
  
  
  
  
