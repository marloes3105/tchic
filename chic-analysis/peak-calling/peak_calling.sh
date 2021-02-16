#!/bin/bash

#SBATCH -t 12:00:00
#SBATCH --mem=10G
#SBATCH -o log.out

source activate cutruntools

for sample in *.bam
do
        bedtools bamtobed -bedpe -i $sample > $sample.bed;
        awk '$1==$4 && $6-$2 < 1000 {print $0}' $sample.bed > $sample.clean.bed;
        cut -f 1,2,6 $sample.clean.bed | sort -k1,1 -k2,2n -k3,3n > $sample.fragments.bed;
        bedtools genomecov -bg -i $sample.fragments.bed -g sizes.genome > $sample.fragments.bedgraph;
        bash SEACR_1.3.sh $sample.fragments.bedgraph 0.01 non stringent $sample;
        rm $sample.bed;
        rm $sample.clean.bed;
        rm $sample.fragments.bed
done
