# step 2

## transcriptome modality
(1) read in plates using the script '0_TCHC_ReadInPlates'. 
(2) perform basic QC on the transcriptome fraction using scanpy, and principle components are calculated and exported in notebook 2a: '2a_TCHIC_QC_until_PCs'.
(3) filter ChIC modality using notebook 2b (see 'ChIC' header)
(4) go back to the transcriptome modality, filter out the QCfail cells and generate the UMAP using notebook 2c: '2c_TCHIC_PCA_to_UMAP'.

## ChIC modality
filter the ChIC modality using notebook 2b: '2b_cell_selection'. This notebook uses the transcriptome PCAs to determine neighbours. It averages over these neighbours and finds highly covered bins. Then it goes back to the individual cell to calculate its enrichment in these pre-selected bins. If the enrichment is high, this means there is a low amount of noise. If the enrichment is low, that means that reads for this cell are present in lots of bins where they shouldn't be, in other words: noise. A cutoff is calculated, cells with a low enrichment score are assigned as 'bad' and cells with a high enrichment score are assigned as 'good'. This threshold can be  changed. You can check the filtering by looking at the produced QC plots. Based on this you can decide whether the suggested cutoff is usable, or has to be adjusted. 

The output of this notebook can be used to filter cells based on ChIC signal, and remove QCfail cells before generating the transcriptome UMAP in notebook 2c.

This script takes a chic csv file and identifies good and bad cells based on enriched bin coverage in neighbour cells. <br/>
#### Necessary inputs:
- merged chic csv file (we use 50kb bins) - can be the output of snakemake 1. `-chic_csv`
- transcriptome csv file with cell names and PCA coordinates (this is used to calculate neighbours) `-trans_csv`


#### If you still need to generate a chic csv table:
Example

```python
 submission.py --nenv -N count_table -y "source activate conda;bamToCountTable.py -bin 50_000 -minMQ 50 --noNames merged_bamfile.bam -sampleTags SM -joinedFeatureTags reference_name -binTag DS --r1only -o output.csv --dedup"
```


