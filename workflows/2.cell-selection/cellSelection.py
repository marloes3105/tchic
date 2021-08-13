#!/usr/bin/env python3

# import packages
import argparse
import os
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
import datetime
import warnings
from pandas.core.common import SettingWithCopyWarning
warnings.simplefilter(action='ignore', category=SettingWithCopyWarning)

###### SETUP
if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Separate good from bad cells in chic csv file""")
    argparser.add_argument('-chic_csv', metavar='chic_csv', help = "ChIC csv file",
                           type=str, required=True)
    argparser.add_argument('-trans_csv', metavar='trans_csv', help = "Transcriptome csv file with umap coordinates",
                           type=str, required=True)
    argparser.add_argument('-cutoff', default=500, type=int)
    argparser.add_argument('-neighbors', default=100, type=int)
    argparser.add_argument( '-o', type=str, help="output data folder", default='./cellSelection/')
    args = argparser.parse_args()

if not os.path.exists(args.o):
    os.makedirs(args.o)


##### INITIATE variables
chic_csv_file = args.chic_csv
trans_csv_file = args.trans_csv
readcutoff = args.cutoff
neighbors = args.neighbors

##### START script
print(str(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+ ': Reading in files')
rawcounts = pd.read_csv(chic_csv_file, index_col=(0,1,2), header= 0, low_memory=False)
coordinates = pd.read_csv(trans_csv_file, index_col=(0), header=0, low_memory=False, delimiter="\t")

# Subset - remove when working
#rawcounts = rawcounts.iloc[:, : 250]
#coordinates = coordinates.reindex(index = rawcounts.columns.tolist())

print(str(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': Done')

## definitions
def calculateCovEnr(df):
    df.loc['distance'] = 100
    df.loc['coverage'] = 0
    df.loc['enrichment'] = 0

    for cellName in df:
        #distancelist = df.copy()
        refpoint = np.array(df[cellName][['V1','V2']])
        for cellName2 in df:
            #calculate Euclidean distances
            otherpoint = np.array(df[cellName2][['V1','V2']])
            df[cellName2]['distance'] = np.linalg.norm(refpoint - otherpoint)

        #select neighbours
        neighbours = df.sort_values('distance', ascending=True, axis = 1).iloc[:, : neighbors]
        neighboursfull = maxnorm.loc[:, maxnorm.columns.isin(neighbours.columns)]

        neighboursfull['sum'] = neighboursfull.sum(axis=1)
        neighboursfull = neighboursfull.sort_values('sum', ascending=False)
        neighborsums = neighboursfull['sum']
        covercounter1  =0
        covercounter2 = 0
        counteraim = sum(neighborsums)*.8

        for counter3 in neighborsums:
            if covercounter1 < counteraim:
                covercounter1 = covercounter1+counter3
                covercounter2 = covercounter2+1

        if cellName in neighboursfull.columns:
            testlist = neighboursfull[cellName]

        coverage = covercounter2/len(neighborsums)

        df[cellName]['coverage'] = (coverage * 100)
        df[cellName]['enrichment'] = (sum(testlist.head(round(len(testlist)*coverage)))/sum(testlist))*100

    df = df.drop(index='distance')

    return df

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

    plt.figure(figsize=(8,8))
    fig = sns.scatterplot(data=celltypes.T, x="V1", y="V2", hue="coverage", palette = 'viridis')
    fig.set(xticklabels=[],yticklabels=[], xlabel='', ylabel='')
    fig.axes.set_title('Coverage on UMAP, round' + str(iteration),fontsize=20)
    sns.despine()
    plt.savefig(f'{args.o}/UMAP_round' + str(iteration) + '.png')

####### First round
print(str(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': Iteration round 1 started')
coordinates.columns = ['V1','V2']
coordinates = coordinates.T
#sort based on readsum
cellsum = rawcounts.sum()
readcutoff2 = cellsum.median() + 2 * cellsum.std()
rawgood = rawcounts.loc[:,(cellsum > readcutoff)]
celltypes = coordinates.loc[:,coordinates.columns.isin(rawgood.columns)]

#normalise data for pseudobulks
maxnorm = (rawgood / rawgood.sum(axis=0))*10000
rawgood = rawgood.fillna(0)
maxnorm = maxnorm.fillna(0)
neighbors = neighbors + 1

celltypes = calculateCovEnr(celltypes)
thresh = celltypes.loc['enrichment'].median() - 2 * (celltypes.loc['enrichment'].std())

#first round result
firstround = celltypes.copy()
firstround.loc['status']='bad'
for cellName in firstround:
    if firstround[cellName]["enrichment"] > thresh:
        firstround[cellName]["status"] = "good"

plots(firstround, color = 'lightskyblue', iteration = '1')

#write tables
print(str(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': Writing iteration 1...')
firstround.to_csv(f'{args.o}/QC_round1.csv', sep = '\t')
firstround.loc['status'].to_csv(f'{args.o}/QCselector_round1.csv', sep = '\t')


####### Second round
print(str(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': Iteration round 2 started')
rawgood = rawcounts.loc[:,rawcounts.columns.isin(firstround.loc[:,firstround.loc["status"]=="good"].columns)]
celltypes = coordinates.loc[:,coordinates.columns.isin(rawgood.columns)]

#normalise data for pseudobulks
maxnorm=(rawgood/rawgood.sum(axis=0))*10000
rawgood=rawgood.fillna(0)
maxnorm=maxnorm.fillna(0)

celltypes = calculateCovEnr(celltypes)
thresh = celltypes.loc['enrichment'].median() - 2 * (celltypes.loc['enrichment'].std())

#second round result
secondround = celltypes.copy()
secondround.loc['status']='bad'
for cellName in secondround:
    if secondround[cellName]["enrichment"] > thresh:
        secondround[cellName]["status"] = "good"

plots(secondround, color = 'cornflowerblue', iteration = '2')

#write tables
print(str(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': Writing iteration 2...')
secondround.to_csv(f'{args.o}/QC_round2.csv', sep = '\t')
secondround.loc['status'].to_csv(f'{args.o}/QCselector_round2.csv', sep = '\t')


####### Third round
print(str(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': Iteration round 3 started')
rawgood = rawcounts.loc[:,rawcounts.columns.isin(secondround.loc[:,secondround.loc["status"]=="good"].columns)]
celltypes = coordinates.loc[:,coordinates.columns.isin(rawgood.columns)]

#normalise data for pseudobulks
maxnorm=(rawgood/rawgood.sum(axis=0))*10000
rawgood=rawgood.fillna(0)
maxnorm=maxnorm.fillna(0)

celltypes = calculateCovEnr(celltypes)

thresh = celltypes.loc['enrichment'].median() - 2 * (celltypes.loc['enrichment'].std())

#third round result
thirdround = celltypes.copy()
thirdround.loc['status']='bad'
for cellName in thirdround:
    if thirdround[cellName]["enrichment"] > thresh:
        thirdround[cellName]["status"] = "good"

plots(thirdround, color = 'navy', iteration = '3')

#write tables
print(str(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': Writing iteration 3...')
thirdround.to_csv(f'{args.o}/QC_round3.csv', sep = '\t')
thirdround.loc['status'].to_csv(f'{args.o}/QCselector_round3.csv', sep = '\t')
print(str(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': All done!')

