#!/usr/bin/env python3

# import packages
from pandas.io.parsers import read_csv
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import matplotlib as mpl
import sys, os
import numpy as np
import argparse
import datetime

## config
if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Normalise ChIC data using percentile normalisation. Currently does not normalise for gene length, so only suitable for TSS tables, not gene body
        tables. This will be added later.""")
    argparser.add_argument('-csv', metavar='csv', help = "Input csv file. Currently only supports TSS tables since gene body tables need to be corrected for gene length.",
                           type=str, required=True)
    argparser.add_argument('-percentile', metavar='percentile', help = "Which percentile to normalise.", default=95,
                           type=int)
    argparser.add_argument( '-o', type=str, help="output csv file.", required=True)
    argparser.add_argument( '--plots', action='store_true', help="If True, make QC plots.")
    argparser.add_argument( '-plot_output', type=str, help="If --plots is True, specify output folder here.", default = './normalisation_plots/')
    argparser.add_argument( '--genebody', action='store_true', help="If True, normalise to gene body length. This requires a -bedfile input too.")
    argparser.add_argument( '-bedfile', type=str, help="If --genebody is True, specify bedfile to use for gene body length normalisation here \
     (this is the same bed file as the one used for generating the gene body table).")
    args = argparser.parse_args()

if args.plots:
    if not os.path.exists(args.plot_output):
        os.makedirs(args.plot_output)

##### initialise files
input_file = args.csv
percentile = args.percentile

#read in table
print(str(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+ ': Reading in files')
input_table = pd.read_csv(input_file, index_col=(0,1,2,3),
                        header= 0, low_memory=False)
input_table = input_table.iloc[1:]
input_table.index.names = ['reference_name','start','end', 'bname']
raw = input_table.droplevel(['start', 'end', 'reference_name'], axis=0)
raw = raw[~raw.index.duplicated(keep='first')]

# normalise
print(str(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+ ': Normalising data')
normalised = raw.copy()
for cell in normalised:
    normalised[cell] = normalised[cell]/np.nanpercentile(normalised[cell],percentile)

# for gene body length normalisation
if args.genebody:
    #open bedfile
    print(str(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+ ': Starting gene body length normalisation')
    bedfile = args.bedfile
    bed = []
    with open(bedfile)as f:
        for line in f:
            bed.append(line.strip().split())
    bed = pd.DataFrame(bed)
    bed = bed.set_index(3)
    #calculate gene length
    bed['length'] = (bed[2].astype(int) - bed[1].astype(int))/100000
    # normalise to gene length and drop empty (gene) rows
    normalised = normalised.join(bed['length'], how='inner')
    for cell in normalised:
        normalised[cell] = normalised[cell]/normalised['length']
    normalised = normalised.drop('length', axis = 1)

# save to csv file
print(str(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+ ': Saving normalised table')
normalised.to_csv(args.o)

#figures for QC
if args.plots:
    print(str(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+ ': Making QC plots')
    fig1_output = f'{args.plot_output}/normalisationQC_perCell.png'
    fig, ax = plt.subplots(4,1,figsize = (10,16))
    plt.suptitle("ChIC mark, mean and max cuts per cell for raw and normalised data")
    sns.histplot(raw.mean(), color="lightgreen", ax=ax[0], kde = True, bins = 125)
    sns.histplot(normalised.mean(), color="mediumseagreen", ax=ax[1], kde = True, bins = 125)
    sns.histplot(raw.max(), color="limegreen", ax=ax[2],kde=True,bins=25 )
    sns.histplot(normalised.max(), color="forestgreen", ax=ax[3],kde=True,bins=25)
    ax[0].set_xlabel("Mean cuts per cell")
    ax[1].set_xlabel("Mean cuts per cell, normalised")
    ax[0].set_ylabel("Number of cells")
    ax[1].set_ylabel("Number of cells")
    ax[2].set_xlabel("Max cuts per cell")
    ax[3].set_xlabel("Max cuts per cell, normalised")
    ax[2].set_ylabel("Number of cells")
    ax[3].set_ylabel("Number of cells")
    plt.savefig(fig1_output)

    fig2_output = f'{args.plot_output}/normalisationQC_perGene.png'
    fig, ax = plt.subplots(4,1,figsize = (10,16))
    plt.suptitle("ChIC mark, mean and max cuts per gene for raw and normalised data")
    sns.histplot(raw.mean(axis=1), color="lightgreen", ax=ax[0], kde = True)
    sns.histplot(normalised.mean(axis=1), color="mediumseagreen", ax=ax[1], kde = True)
    sns.histplot(raw.max(axis=1), color="limegreen", ax=ax[2],kde=True,bins=25)
    sns.histplot(normalised.max(axis=1), color="forestgreen", ax=ax[3],kde=True,bins=25)
    ax[0].set_xlabel("Mean cuts per gene")
    ax[1].set_xlabel("Mean cuts per gene, normalised")
    ax[0].set_ylabel("Number of genes")
    ax[1].set_ylabel("Number of genes")
    ax[2].set_xlabel("Max cuts per gene")
    ax[3].set_xlabel("Max cuts per gene, normalised")
    ax[2].set_ylabel("Number of genes")
    ax[3].set_ylabel("Number of genes")
    plt.savefig(fig2_output)

print(str(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+ ': All done!')
