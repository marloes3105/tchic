#!/bin/bash

#SBATCH -t 24:00:00
#SBATCH --mem=200G
#SBATCH -o log.out

source activate condaP1

#take merged bam file from SCMO and split into good and bad cells, text file comes from R-script to identify high quality cells. Adjust file names
bamExtractSamples.py /hpc/hub_oudenaarden/Marloes/data/tchic/gastruloids/OUD5650/k27me3/processed_chic/merged_tagged.bam rep1H3K27me3.txt -o K27me3_.bam

#take only good cells and split bam file into cell types. txt file comes from python script for clustering and cell typing. Adjust file names.
bamExtractSamples.py K27me3_good.bam VAN5399_OUD5650_cells_clusters_celltypes.txt -o K27me3_.bam

#make count tables for all bam files in folder surrounding TSS. bed file with TSS regions has to be in same folder
for file in *.bam; do echo "submission.py -y -time 3 -m 50 -t 2 -s ./cluster -sched slurm -N bin200 --nenv \"source activate condaP1;bamToCountTable.py --filterXA -minMQ 50 $file -o $file.TSS.csv -sampleTags SM -joinedFeatureTags reference_name --dedup --r1only -bedfile TSS_mus_musculus_cleaned_10kb.bed\"";done | sh

#make cover plots for all bam files in the folder
for file in *.bam; do echo "submission.py -y -time 3 -m 50 -t 2 -s ./cluster -sched slurm -N bin200 --nenv \"source activate condaP1;bamCoverage --binSize 200 --smoothLength 1000 --ignoreDuplicates --centerReads --normalizeUsing CPM -b $file -o $file.bw\"";done | sh




