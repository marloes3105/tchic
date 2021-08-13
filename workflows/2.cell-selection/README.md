## This script takes a chic csv file and identifies good and bad cells based on enriched bin coverage in neighbour cells. <br/>
### Necessary inputs:
- merged chic csv file (we use 50kb bins) - can be the output of snakemake 1. `-chic_csv`
- transcriptome csv file with cell names and umap coordinates (this is used to calculate neighbours) `-trans_csv`
- output folder - does not have to exist yet. `-o`

### optional:<br/>
- number of neighbours to use, default: 100. `-neighbors`
- cutoff to use, default: 500. `-cutoff`

## What does it do?
This script uses the transcriptome UMAP coordinates to determine neighbours. It averages over these neighbours and finds highly covered bins. Then it goes back to the individual cell to calculate its coverage in these pre-selected bins. If the coverage is high, this means there is a low amount of noise. If the coverage is low, that means that reads for this cell are present in lots of bins where they shouldn't be, in other words: noise. A cutoff is calculated, cells with a low enrichment score are assigned as 'bad' and cells with a high enrichment score are assigned as 'good'. This whole process is repeated 3 times to filter out varying degrees of bad cells. After each round, you can check the filtering by looking at the produced plots ('coverageplot' and UMAP plot). Based on this you can decide whether one round of selection is enough, or if you need 2 or all 3 for filtering out bad cells. The full QC metrics are written to 'QC_roundx.csv'. 'Qselector_roundx.csv' is used as input for the next step of this workflow (step 3 - snakemake).

### To run this script: 

```python
source activate conda;python cellSelection.py -chic_csv Gast_d3-7_H3K4me3.50kb.csv -trans_csv 20210623_rep234_day34567_cells_umapCoordinates.csv -o cellSelection
```

### To submit and run this script: 

```python
submission.py --nenv -N cellselection -m 20 -time 80 -y "source activate conda;python cellSelection.py -chic_csv Gast_d3-7_H3K4me3.50kb.csv -trans_csv 20210623_rep234_day34567_cells_umapCoordinates.csv -o cellSelection"
```
