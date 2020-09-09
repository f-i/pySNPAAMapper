#!/usr/bin/env python
'''
Algorithm for preprocessing the gene structure to build annotation structure for
each exon

USAGE:
    python Algorithm_preprocessing_exon_annotation_regulatoryRegion.py ChrAll_knownGene.txt.exons

'''

import sys
from fileinput import FileInput
import pandas as pd

def argsCheck(numArgs):
    '''Checks if in proper number of arguments are passed; gives instructions on
    proper use'''
    if len(sys.argv) < numArgs or len(sys.argv) > numArgs:
        print("\nPlease supply ONLY ONE file to parse.\n\nCorrect Syntax is:\n    python " + sys.argv[0] + " ChrAll_knownGene.txt.exons\n\nThe 4 output files will be the input file name with \".[strings]\" appended.\n")
        exit(1)  # Aborts program. (exit(1) indicates that an error occurred)

def link(filename, who, a):
    '''filename is the name of input file as the only argument for running this
    script; who is a string from a list: 'cds', 'intron', 'utr5', 'utr3',
    'upstream' and 'downstream'; a is a Pandas DataFrame read from filename'''
    start = who + 'Start'
    end = who + 'End'
    b = a[((a[start].notnull())|(a[end].notnull()))&(~a.chrom.str.contains('_'))]
    b.loc[:,start] = b[start].astype(int)
    b.loc[:,end] = b[end].astype(int)
    new_col = 'chrom_' + who + 'Start'
    b.loc[:,new_col] = b.chrom + '_'+b[start].astype(str)
    c = b[[new_col,end,'name']]
    print('Saving to file '+filename+'.'+who+'_link')
    c.to_csv(filename+'.'+who+'_link',
             index=False,sep='\t',header=False)
    d = b.groupby('chrom')[start].apply(list).apply(sorted)
    print('Saving to file '+filename+'.'+who)
    d.to_csv(filename+'.'+who,
             index=True,sep='\t',header=False)
    with FileInput(files=[filename+'.'+who], inplace=True) as f:
        for line in f:
            print(line.replace('[', '').replace(']\n', '').replace(', ', '\t'))
    e = b.groupby(new_col, sort=False).agg({end:max,'name':'first'})
    print('Saving to file '+filename+'.'+who+'_link_shrink')
    e[end].to_csv(filename+'.'+who+'_link_shrink',
                  index=True,sep='\t',header=False)
    print('Saving to file '+filename+'.'+who+'_gene')
    e['name'].to_csv(filename+'.'+who+'_gene',
                     index=True,sep='\t',header=False)

# Uses the function defined above to check if the number of arguments is correct
argsCheck(2)

inFile = sys.argv[1]  # Stores file one for input checking; note that sys.argv[0] is this script file name

print("Reading your data file "+inFile+" ...")
a = pd.read_csv(inFile,delimiter="\t")

whos = ['cds', 'intron', 'utr5', 'utr3', 'upstream', 'downstream']
for i in whos:
    link(inFile, i ,a)
