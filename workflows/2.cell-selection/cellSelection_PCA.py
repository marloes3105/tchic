#!/usr/bin/env python3

# import packages
import argparse
import os
import numpy as np
import scipy
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
import datetime
import warnings
from pandas.core.common import SettingWithCopyWarning
import multiprocessing
from sklearn.neighbors import NearestNeighbors
warnings.simplefilter(action='ignore', category=SettingWithCopyWarning)


###### SETUP
if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Separate good from bad cells in chic csv file""")
    argparser.add_argument('-chic_csv', metavar='chic_csv', help = "ChIC csv file",
                           type=str, required=True)
    argparser.add_argument('-trans_csv', metavar='trans_csv', help = "Transcriptome based PCA csv file",
                           type=str, required=True)
    argparser.add_argument('-cutoff', default=500, type=int)
    argparser.add_argument('-neighbors', default=100, type=int)
    argparser.add_argument( '-o', type=str, help="output data folder", default='./cellSelection/')
    argparser.add_argument('-subset', type = int, help = "How big the parallelized subsets will be. E.g. if you type -subset 500, it will parallelized into \
                            chunks of 500 cells.", default = 500)
    argparser.add_argument('-pools', type = int, help = 'number of pools to run in parallel', default = 8)
    args = argparser.parse_args()

if not os.path.exists(args.o):
    os.makedirs(args.o)


##### INITIATE variables
chic_csv_file = args.chic_csv
trans_csv_file = args.trans_csv
readcutoff = args.cutoff
neighbors = args.neighbors
subsets = args.subset
pools = args.pools

##### START script
print(str(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+ ': Reading in files')
rawcounts = pd.read_csv(chic_csv_file, index_col=(0,1,2), header= 0, low_memory=False)
coordinates = pd.read_csv(trans_csv_file, index_col=(0), header=0, low_memory=False)




# Subset - remove when working
#rawcounts = rawcounts.iloc[:, : 250]
#coordinates = coordinates.reindex(index = rawcounts.columns.tolist())

print(str(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': Done')



## definitions
def calculateCovEnr(df):
    
    
    returndf=pd.DataFrame(columns={'QC','coverage','enrichment'}, index=cellnames)

    for name in cellnames:
        #select neighbours
        subset=df.loc[cellnames[neighborslist.loc[name]]].T
        subset['sum']=(subset/subset.sum(axis=0)).sum(axis=1)
        subset = subset.sort_values('sum', ascending=False)
        subsetsums = subset['sum']
        covercounter1  =0
        covercounter2 = 0
        counteraim = sum(subsetsums)*.8
        for counter2 in subsetsums:
            if covercounter1 < counteraim:
                covercounter1 = covercounter1+counter2
                covercounter2 = covercounter2+1

        coverage = covercounter2/50000
        testlist = subset[name]

        returndf.loc[name,'coverage'] = (coverage * 100)
        returndf.loc[name,'enrichment'] = ((testlist.head(round(len(testlist)*coverage))).sum()/(testlist).sum())*100
    return returndf

def plots(data, color = 'cornflowerblue', iteration = 'x'):
    mpl.rcParams['figure.dpi'] = 200

    plt.figure(figsize=(10,6))
    fig = sns.histplot(data=data.T, x='enrichment', color = color)
    fig.axes.set_title("After round " + str(iteration),fontsize=30)
    fig.set_xlabel("Enrichment",fontsize=20)
    fig.set_ylabel("Number of cells",fontsize=20)
    fig.tick_params(labelsize=15)
    fig = plt.axvline(x=thresh, color='crimson',lw = 2)

    sns.despine()
    plt.savefig(f'{args.o}/coveragePlot_round' + str(iteration) + '.png')


####### start
print(str(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': Iteration round 1 started')


#sort based on readsum
coordinates = coordinates.T
cellsum = rawcounts.sum()
rawgood = rawcounts.loc[:,(cellsum > readcutoff)]
celltypes = coordinates.loc[:,coordinates.columns.isin(rawgood.columns)]
rawgood = rawgood.loc[:,rawgood.columns.isin(celltypes.columns)]


#generate neighbors
celltypes=celltypes.T
celltypes_sparse=scipy.sparse.csr_matrix(celltypes.values)
nbrs = NearestNeighbors(n_neighbors=neighbors, algorithm='auto').fit(celltypes_sparse)
distances, indices = nbrs.kneighbors(celltypes_sparse)
neighborslist=pd.DataFrame(indices)
neighborslist.index=celltypes.index
cellnames=neighborslist.index


#normalise data for pseudobulks
maxnorm = (rawgood / rawgood.sum(axis=0))*10000
rawgood = rawgood.fillna(0)
maxnorm = maxnorm.fillna(0)


resultstable = calculateCovEnr(maxnorm.T).T
thresh = resultstable.loc['enrichment'].median() - 2 * (resultstable.loc['enrichment'].std())

#first round result
firstround = resultstable.copy()
firstround.loc['QC']='bad'
for cellName in firstround:
    if firstround[cellName]["enrichment"] > thresh:
        firstround[cellName]["QC"] = "good"


plots(firstround, color = 'lightskyblue', iteration = '1')

#write tables
print(str(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': Writing iteration 1...')
firstround.to_csv(f'{args.o}/QC.csv', sep = '\t')
firstround.loc['QC'].to_csv(f'{args.o}/QCselector.csv', sep = '\t')


print(str(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': All done!')
