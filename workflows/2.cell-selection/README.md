## This script takes a chic csv file and identifies good and bad cells based on enriched bin coverage in neighbour cells. <br/>
### Necessary inputs:
- merged chic csv file (we use 50kb bins) - can be the output of snakemake 1. `-chic_csv`
- transcriptome csv file with cell names and PCA coordinates (this is used to calculate neighbours) `-trans_csv`

## What does it do?
This notebook uses the transcriptome PCAs to determine neighbours. It averages over these neighbours and finds highly covered bins. Then it goes back to the individual cell to calculate its enrichment in these pre-selected bins. If the enrichment is high, this means there is a low amount of noise. If the enrichment is low, that means that reads for this cell are present in lots of bins where they shouldn't be, in other words: noise. A cutoff is calculated, cells with a low enrichment score are assigned as 'bad' and cells with a high enrichment score are assigned as 'good'. This threshold can be  changed. You can check the filtering by looking at the produced QC plots. Based on this you can decide whether the suggested cutoff is usable, or has to be adjusted. 

The output of this notebook can be used to filter cells based on ChIC signal, and remove QCfail cells before generating the transcriptome UMAP in notebook 2b.


## If you still need to generate a chic csv table:
Example

```python
 submission.py --nenv -N count_table -y "source activate conda;bamToCountTable.py -bin 50_000 -minMQ 50 --noNames merged_bamfile.bam -sampleTags SM -joinedFeatureTags reference_name -binTag DS --r1only -o output.csv --dedup"
```


