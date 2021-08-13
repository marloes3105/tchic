#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import seaborn as sns
import numpy as np
import pkg_resources
import pysam
from pysamiterators import MatePairIterator
from itertools import product
from collections import Counter
from singlecellmultiomics.bamProcessing import sorted_bam_file, merge_bams
from singlecellmultiomics.utils import reverse_complement, is_main_chromosome, pool_wrapper
from singlecellmultiomics.bamProcessing.bamFunctions import get_contigs_with_reads
from multiprocessing import Pool
import random
from more_itertools import chunked
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn import datasets, metrics, model_selection
import argparse
import os
import sys

from singlecellmultiomics.barcodeFileParser import barcodeFileParser
bparse = barcodeFileParser.BarcodeParser(pkg_resources.resource_filename(
            'singlecellmultiomics',
            'modularDemultiplexer/barcodes/'))
id_to_cs2_barcode = { v:k for k,v in bparse.barcodes['celseq2'].items() }

## config
if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Clean up ChIC bam file, remove transcriptome contamination in T-ChIC libraries.""")
    argparser.add_argument('-chic_bam', metavar='chic_bamfile', help = "Positive training dataset. This is a ChIC bam file used to detect ChIC-specific reads.",
                           type=str, required=True)
    argparser.add_argument('-vasa_bam', metavar='vasa_bamfile', help = "Negative training dataset. This is a VASA bam file used to detect transcriptome-specific reads.",
                           type=str, required=True)
    argparser.add_argument('-tchic_bam', metavar='tchic_bamfile', help = "Input bam file that needs to be cleaned up.",
                           type=str, required=True)
    argparser.add_argument( '-o', type=str, help="output data folder", default='./')
    argparser.add_argument('-selection', type=str, help="Select which features you want to use for training and spliting the data. Example: -selection  ['TA','TC','TG','TT','CT','CC','CA']", default=None)
    args = argparser.parse_args()

if not os.path.exists(args.o):
    os.makedirs(args.o)

##### initialise files
positive_train_path = args.chic_bam
negative_train_path = args.vasa_bam
tchic_bamfile = args.tchic_bam

twomers = { a+b:i for i,(a,b) in enumerate( product('ACGTN',repeat=2) )}

##### functions
def read_to_unmapped_vect(read):
    if read is None:
        return [0]*15
    # Clipped content:
    cnt = Counter([
        read.query_sequence[qpos] for qpos, rpos
        in read.get_aligned_pairs(matches_only=False) if rpos is None]
    )
    # Aligned content:
    cnt_al = Counter([
        read.query_sequence[qpos] for qpos, rpos
        in read.get_aligned_pairs(matches_only=True) if rpos is not None]
    )
    chic_barcode_index = int( read.get_tag('bi') )
    expected_bc = id_to_cs2_barcode[chic_barcode_index]
    # Determine if a CS2 barcode is present in the read:
    sequence = read.query_sequence
    if read.is_reverse:
        sequence = reverse_complement(sequence)

    bcode_count = sequence.count(expected_bc)
    if bcode_count>0:
        bcode_start = sequence.index(expected_bc)
    else:
        bcode_start = 0

    random_barcode_hits = 0
    for b,target in bparse.barcodes['celseq2'].items():
        if b in sequence:
            random_barcode_hits+=1

    if read.is_reverse:
        return [ cnt_al['T'], cnt_al['G'], cnt_al['C'], cnt_al['A'],  cnt['T'], cnt['G'], cnt['C'], cnt['A'], read.is_proper_pair, read.is_reverse, read.is_qcfail, read.mapping_quality, bcode_count, bcode_start, random_barcode_hits ]
    else:
        return [ cnt_al['A'], cnt_al['C'], cnt_al['G'], cnt_al['T'],  cnt['A'], cnt['C'], cnt['G'], cnt['T'], read.is_proper_pair, read.is_reverse, read.is_qcfail, read.mapping_quality, bcode_count, bcode_start, random_barcode_hits ]

read_labels = ['aligned_A','aligned_C','aligned_G','aligned_T','clipped_A','clipped_C','clipped_G','clipped_T','proper_pair','reverse','qcfail' ,'mq', 'bcode_count', 'bcode_start','any_cs2']

labels = []
for mate in ('R1','R2'):
    for label in read_labels:
        labels.append( f'{mate}_{label}')

labels.append('Fragment size')

for mer in twomers:
    labels.append(mer)

def to_feature_vector(R1,R2, read):
    x = [0]*len(twomers)
    x[twomers[read.get_tag('lh')]] = 1
    fs = 0
    if R1 is not None and R2 is not None:
        fs = max(R1.reference_end, R2.reference_end) - min(R1.reference_start, R2.reference_start)

    return read_to_unmapped_vect(R1) + read_to_unmapped_vect(R2)  + [fs] + x

def sample_bam_to_feature_matrix(bam_path, n):
    X = []
    with pysam.AlignmentFile(bam_path) as t:
        for i,(R1, R2) in enumerate(MatePairIterator(t, ignore_collisions=True)):
            if random.random()>0.05:
                continue
            if len(X)>n:
                break
            for read in (R1,R2):
                if read is not None:
                    break
            if read.is_duplicate:
                continue
            X.append( to_feature_vector(R1,R2, read) )
    return X

colormap = plt.get_cmap('RdYlBu_r')
colormap.set_bad((0,0,0))

def add_color_tag_to_read( read, value ):
    if read is None:
        return
    try:
        cfloat = colormap( (value-0.5)*2 )[:3]
    except Exception as e:
        cfloat = colormap._rgba_bad[:3]
    read.set_tag('YC', '%s,%s,%s' % tuple((int(x * 255) for x in cfloat)))

def write_classification(reads,y, chic, transcriptome):
     for (r1,r2), (chic_prob, trans_prob) in zip(reads, y):
        if r1 is not None:
            r1.set_tag('pC', chic_prob)
            r1.set_tag('pT', trans_prob)
        if r2 is not None:
            r2.set_tag('pC', chic_prob)
            r2.set_tag('pT', trans_prob)
        if chic_prob>trans_prob:
            target = chic
            add_color_tag_to_read(r1, chic_prob)
            add_color_tag_to_read(r2, chic_prob)
        else:
            target = transcriptome
            add_color_tag_to_read(r1, trans_prob)
            add_color_tag_to_read(r2, trans_prob)
        if r1 is not None:
            target.write(r1)
        if r2 is not None:
            target.write(r2)

def classify_bam(bam_in, classifier, contig, batch_size = 100_000):
    with pysam.AlignmentFile(bam_in, threads=8) as t:
        out_trans = f'{args.o}/transcriptome_filtered_{contig}.bam'
        out_chic = f'{args.o}/chic_filtered_{contig}.bam'
        with sorted_bam_file(out_trans, bam_in,fast_compression=True,) as transcriptome, \
             sorted_bam_file(out_chic, bam_in,fast_compression=True) as chic:
            X = []
            reads = []
            #5,232,644-5,559,426
            for i,(R1, R2) in enumerate(MatePairIterator(t, contig=contig, ignore_collisions=True)):
                for read in (R1,R2):
                    if read is not None:
                        break
                X.append( to_feature_vector(R1,R2, read) )
                reads.append( (R1,R2) )
                if len(X)>=batch_size:
                    y = classifier.predict_proba(pd.DataFrame(X).iloc[:,indexes].values.tolist())
                    write_classification(reads,y,chic,transcriptome)
                    X = []
                    reads = []

            X = pd.DataFrame(X).iloc[:,indexes]
            X = X.values.tolist()
            # Write remainder
            if len(X)>0:
                y = classifier.predict_proba(X)
                write_classification(reads,y,chic,transcriptome)
                X = []
                reads = []
    return out_trans, out_chic

##### start
X = []
y = []
x = sample_bam_to_feature_matrix(negative_train_path, 50_000)
_y = len(x)*['T']
X += x
y += _y
x = sample_bam_to_feature_matrix(positive_train_path, 50_000)
_y = len(x)*['C']
X += x
y += _y
X = pd.DataFrame(X,columns=labels)

y=np.array(y)

clf = RandomForestClassifier(n_jobs=16, oob_score=True)

selection = args.selection
selection  = selection.split(',')

if selection == None:
    selection = labels

indexes = [X.columns.get_loc(c) for c in selection if c in X]

X = X.iloc[:,indexes]

### PLOTS
## Plot2
fig2_output = f'{args.o}/chic_cleanup_TrueFalsePositive.png'

if False:
    {
        'all features':np.mean( cross_val_score(clf,X,y) ),
 #       'without cs2 barcode':np.mean( cross_val_score(clf,X[[col for col in X if not any( (c in col  for c in ('bcode_count', 'bcode_start','any_cs2')))]],y)),
  #      'only cs2 barcode':np.mean( cross_val_score(clf,X[[col for col in X if any( (c in col  for c in ('bcode_count', 'bcode_start','any_cs2')))]],y))
    }


fig, ax = plt.subplots()
for x_dat,y_dat,name in [

#    (X[[col for col in X if any( (c in col  for c in ('bcode_count', 'bcode_start','any_cs2')))]],y,'CS2 barcode presence'),
 #   (X[[col for col in X if any( (c in col  for c in ('aligned', 'clipped')))]],y,'Aligned and clipped bases'),
  #  (X[[col for col in X if any( (c in col  for c in twomers))]],y,'Ligation motif'),
  #  (X[[col for col in X if not any( (c in col  for c in ('bcode_count', 'bcode_start','any_cs2','aligned')))]],y,'Original features'),
    (X,y,'All features'),

]:

    X_train, X_test, y_train, y_test = model_selection.train_test_split(
        x_dat, y_dat, random_state=0)
    clf.fit(X_train, y_train)

    metrics.plot_roc_curve(clf, X_test, y_test,name=name,ax=ax,)
sns.despine()
plt.tight_layout()
plt.savefig(fig2_output)

## FeatureImportance Figure
clf.fit(X,y)
fig3_output = f'{args.o}/chic_cleanup_featureImportance.png'
mpl.rcParams['figure.dpi'] = 200
pd.DataFrame( clf.feature_importances_, index=X.columns.tolist(), columns=['Feature importance'] ).plot.barh(figsize=(5,12))
plt.grid()
sns.despine()
plt.tight_layout()
plt.savefig(fig3_output)
###

with pysam.AlignmentFile(tchic_bamfile) as a:
    contigs = [q for q in a.references if is_main_chromosome(q)]

chic_bams = []
trans_bams = []

clf.n_jobs = 1

with Pool() as workers:
    for trans_bam, chic_bam in workers.imap(pool_wrapper, ((classify_bam, {'classifier':clf, 'bam_in':tchic_bamfile,
                                                          'contig':contig })
                                                          for contig in contigs )):
        chic_bams.append(chic_bam)
        pysam.index(chic_bam)
        trans_bams.append(trans_bam)
        pysam.index(trans_bam)

merge_bams(chic_bams, f'{args.o}/chic_classified.bam')
merge_bams(trans_bams, f'{args.o}/trans_classified.bam')
