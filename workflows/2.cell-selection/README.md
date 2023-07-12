## This script takes a chic csv file and identifies good and bad cells based on enriched bin coverage in neighbour cells. <br/>
### Necessary inputs:
- merged chic csv file (we use 50kb bins) - can be the output of snakemake 1. `-chic_csv`
- transcriptome csv file with cell names and umap coordinates (this is used to calculate neighbours) `-trans_csv`
- output folder - does not have to exist yet. `-o`

### optional:<br/>
- number of neighbours to use, default: 100. `-neighbors`
- cutoff to use, default: 500. `-cutoff`

## What does it do?
This script uses the transcriptome PCAs to determine neighbours. It averages over these neighbours and finds highly covered bins. Then it goes back to the individual cell to calculate its enrichment in these pre-selected bins. If the enrichment is high, this means there is a low amount of noise. If the enrichment is low, that means that reads for this cell are present in lots of bins where they shouldn't be, in other words: noise. A cutoff is calculated, cells with a low enrichment score are assigned as 'bad' and cells with a high enrichment score are assigned as 'good'. This threshold can be manually changed after integrating the results in the jupyter noteboook. You can check the filtering by looking at the produced histogram. Based on this you can decide whether the suggested cutoof is usable, or has to be adjusted. The full QC metrics are written to 'QC.csv'. 'Qselector.csv' is used as input for the next step of this workflow (step 3 - snakemake).

### To run this script: 

```python
source activate conda;python cellSelection.py -chic_csv Gast_d3-7_H3K4me3.50kb.csv -trans_csv 20210623_rep234_day34567_cells_PCAs.csv -o cellSelection
```

### To submit and run this script: 

```python
submission.py --nenv -N cellselection -m 20 -time 80 -y "source activate conda;python cellSelection.py -chic_csv Gast_d3-7_H3K4me3.50kb.csv -trans_csv 20210623_rep234_day34567_cells_PCAs.csv -o cellSelection"
```

## If you still need to generate a chic csv table:
Example

```python
 submission.py --nenv -N count_table -y "source activate conda;bamToCountTable.py -bin 50_000 -minMQ 50 --noNames merged_bamfile.bam -sampleTags SM -joinedFeatureTags reference_name -binTag DS --r1only -o output.csv --dedup"
```


