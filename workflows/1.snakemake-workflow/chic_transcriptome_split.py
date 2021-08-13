#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import pandas as pd
import pysam
from pysamiterators import MatePairIterator
from itertools import product
from collections import Counter
from singlecellmultiomics.bamProcessing import sorted_bam_file, merge_bams
import random
from more_itertools import chunked
from sklearn.ensemble import RandomForestClassifier
from singlecellmultiomics.utils import is_main_chromosome, pool_wrapper
from singlecellmultiomics.bamProcessing import merge_bams
from multiprocessing import Pool
import seaborn as sns
import matplotlib as mpl
from singlecellmultiomics.bamProcessing.bamFunctions import get_contigs_with_reads
import os

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
    args = argparser.parse_args()

if not os.path.exists(args.o):
    os.makedirs(args.o)

pos_train = args.chic_bam
neg_train = args.vasa_bam
tchic_bamfile = args.tchic_bam

twomers = { a+b:i for i,(a,b) in enumerate( product('ACGTN',repeat=2) )}

def read_to_unmapped_vect(read):
    if read is None:
        return [0]*8
    cnt = Counter([
        read.query_sequence[qpos] for qpos, rpos 
        in read.get_aligned_pairs(matches_only=False) if rpos is None]
    )
    if read.is_reverse:
        return [ cnt['T'], cnt['G'], cnt['C'], cnt['A'], read.is_proper_pair, read.is_reverse, read.is_qcfail, read.mapping_quality ]
    else:
        return [ cnt['A'], cnt['C'], cnt['G'], cnt['T'], read.is_proper_pair, read.is_reverse, read.is_qcfail, read.mapping_quality ]
    
read_labels = ['A','C','G','T','proper_pair','reverse','qcfail' ,'mq']

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
        for i,(R1, R2) in enumerate(MatePairIterator(t)):      
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
            for i,(R1, R2) in enumerate(MatePairIterator(t, contig=contig)):
                for read in (R1,R2):
                    if read is not None:
                        break
                X.append( to_feature_vector(R1,R2, read) )
                reads.append( (R1,R2) )
                if len(X)>=batch_size:
                    y = classifier.predict_proba(X)
                    write_classification(reads,y,chic,transcriptome)
                    X = []
                    reads = []
            # Write remainder
            if len(X)>0:
                y = classifier.predict_proba(X)
                write_classification(reads,y,chic,transcriptome)
                X = []
                reads = []
    return out_trans, out_chic

X = []
y = []
x = sample_bam_to_feature_matrix(neg_train, 50_000)
_y = len(x)*['T']
X += x
y += _y
x = sample_bam_to_feature_matrix(pos_train, 50_000)
_y = len(x)*['C']
X += x
y += _y
    
clf = RandomForestClassifier(n_jobs=16, oob_score=True)  
clf.fit(X,y)

fig_output = f'{args.o}/chic_cleanup_featureImportance.png'

mpl.rcParams['figure.dpi'] = 200
pd.DataFrame( clf.feature_importances_, index=labels, columns=['Feature importance'] ).plot.barh(figsize=(5,8))
plt.grid()
sns.despine()
plt.tight_layout()
plt.savefig(fig_output)

with pysam.AlignmentFile(pos_train) as a:
    contigs = [q for q in a.references if is_main_chromosome(q)]

chic_bams = []
trans_bams = []

clf.n_jobs = 1

with Pool() as workers:
    for trans_bam, chic_bam in workers.imap(pool_wrapper, ((classify_bam, {'classifier':clf, 'bam_in':tchic_bamfile,
                                                          'contig':contig })
                                                          for contig in get_contigs_with_reads(tchic_bamfile))):
        chic_bams.append(chic_bam)
        pysam.index(chic_bam)
        trans_bams.append(trans_bam)
        pysam.index(trans_bam)

merge_bams(chic_bams, f'{args.o}/chic_classified.bam')
merge_bams(trans_bams, f'{args.o}/trans_classified.bam')
    
