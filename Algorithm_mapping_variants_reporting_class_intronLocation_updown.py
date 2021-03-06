#!/usr/bin/env python
'''
Algorithm for mapping identified variants onto the genomic location and
reporting the hit class

USAGE:
    python Algorithm_mapping_variants_reports_class_intronLocation.py ChrAll_knownGene.txt.exons VCF_input_file_in_tab_delimited_format.txt IntronExon_boundary_in_bp
'''

import sys
import pandas as pd


def link(who, a):
    '''who is a string from a list: 'cds', 'intron', 'utr5', 'utr3', 'upstream'
    and 'downstream'; a is a Pandas DataFrame of 'ChrAll_knownGene.txt.exons' '''
    start = who + 'Start'
    end = who + 'End'
    b = a[((a[start].notnull())|(a[end].notnull()))&(~a.chrom.str.contains('_'))]
    b.loc[:,start] = b[start].astype(int)
    b.loc[:,end] = b[end].astype(int)
    new_col = 'chrom_' + who + 'Start'
    b.loc[:,new_col] = b.chrom + '_'+b[start].astype(str)
    d = b.groupby('chrom')[start].apply(list).apply(sorted)
    e = b.groupby(new_col, sort=False).agg({end:max,'name':'first'})
    return d,e

inFile = sys.argv[1]  # Stores file one for input checking; note that sys.argv[0] is this script file name
print("Reading your first data file "+inFile+" ...")
a = pd.read_csv(inFile, delimiter="\t")

inFile0 = sys.argv[2]
print("Reading your second data file, the VCF input file in tab delimited format "+inFile+" ...")
d = pd.read_csv(inFile0, delimiter="\t", header=26)  #header numbers start from 0
d["HIT"] = ""
d["gene"] = ""

if len(sys.argv) == 4:
    intronOpt=1
    exonbuffer=sys.argv[3]
    print("The program assumes that you DO want to report how far the variant falls in the exon boundary.\nOnly variants flanking its nearby exon with <= {} bp is reported\n".format(exonbuffer))
elif len(sys.argv) == 3:
    intronOpt=-1
    print("The program assumes that you do NOT want to report how far the variant falls in the exon boundary.\n")
else:
    print("The input commands do not meet the requirement. Please see the README file and try again.\n")
    exit(1)  # Aborts program. (exit(1) indicates that an error occurred)

df=pd.DataFrame()
for i in range(len(d['POS'])):
    whos = ['cds', 'intron', 'utr5', 'utr3', 'upstream', 'downstream']
    for j in whos:
        b,c = link(j,a)
        tmp_n=d['POS'][i]
        tmpchr=d['#CHROM'][i]
        if tmpchr in b.index:
            try:
                ri=b[tmpchr].index(tmp_n)
            except ValueError:
                tmp=b[tmpchr][:]  #this is how we clone a list
                tmp.append(tmp_n)
                ri=sorted(tmp).index(tmp_n)-1
        if j=='cds':
            if tmp_n>b[tmpchr][ri] and \
            tmp_n<=c.loc[tmpchr+'_'+str(b[tmpchr][ri]),j+'End'] and \
            b[tmpchr][ri]!='NA' and \
            c.loc[tmpchr+'_'+str(b[tmpchr][ri]),j+'End']!='NA':
                d.loc[i,'HIT']="CDSHIT"
                d.loc[i,'gene']=c.loc[tmpchr+'_'+str(b[tmpchr][ri]),'name']
                df = df.append(d.iloc[[i]])
                print(d.iloc[[i]].to_string(index=False, header=False),
                      "CDSHIT",c.loc[tmpchr+'_'+str(b[tmpchr][ri]),'name'])
        if j=='intron':
            if intronOpt==1:
                if tmp_n>=b[tmpchr][ri] and \
                tmp_n<(b[tmpchr][ri]+exonbuffer) or \
                tmp_n<=(c.loc[tmpchr+'_'+str(b[tmpchr][ri]),j+'End']+1) and \
                tmp_n>(c.loc[tmpchr+'_'+str(b[tmpchr][ri]),j+'End']+1-exonbuffer) and \
                b[tmpchr][ri]!='NA' and \
                c.loc[tmpchr+'_'+str(b[tmpchr][ri]),j+'End']!='NA':
                    if (tmp_n-b[tmpchr][ri])<exonbuffer:
                        hitbuffer=tmp_n-b[tmpchr][ri]+1
                        direction="right"
                    if (c.loc[tmpchr+'_'+str(b[tmpchr][ri]),j+'End']+1-tmp_n)<exonbuffer:
                        hitbuffer=c.loc[tmpchr+'_'+str(b[tmpchr][ri]),j+'End']+1-tmp_n+1
                        direction="left"
                    d.loc[i,'HIT']="INTRONHIT."+str(hitbuffer)+"."+direction
                    d.loc[i,'gene']=c.loc[tmpchr+'_'+str(b[tmpchr][ri]),'name']
                    df = df.append(d.iloc[[i]])
                    print(d.iloc[[i]].to_string(index=False, header=False),
                          "INTRONHIT."+str(hitbuffer)+"."+direction,
                          c.loc[tmpchr+'_'+str(b[tmpchr][ri]),'name'])
            if intronOpt==-1:
                if tmp_n>=b[tmpchr][ri] and \
                tmp_n<=(c.loc[tmpchr+'_'+str(b[tmpchr][ri]),j+'End']+1) and \
                b[tmpchr][ri]!='NA' and \
                c.loc[tmpchr+'_'+str(b[tmpchr][ri]),j+'End']!='NA':
                    d.loc[i,'HIT']="INTRONHIT"
                    d.loc[i,'gene']=c.loc[tmpchr+'_'+str(b[tmpchr][ri]),'name']
                    df = df.append(d.iloc[[i]])
                    print(d.iloc[[i]].to_string(index=False, header=False),
                          "INTRONHIT",
                          c.loc[tmpchr+'_'+str(b[tmpchr][ri]),'name'])
        if j=='utr5' or j=='upstream':
            if tmp_n>b[tmpchr][ri] and \
            tmp_n<=(c.loc[tmpchr+'_'+str(b[tmpchr][ri]),j+'End']+1) and \
            b[tmpchr][ri]!='NA' and \
            c.loc[tmpchr+'_'+str(b[tmpchr][ri]),j+'End']!='NA':
                d.loc[i,'HIT']=j.upper()+"HIT"
                d.loc[i,'gene']=c.loc[tmpchr+'_'+str(b[tmpchr][ri]),'name']
                df = df.append(d.iloc[[i]])
                print(d.iloc[[i]].to_string(index=False, header=False),
                      j.upper()+"HIT",c.loc[tmpchr+'_'+str(b[tmpchr][ri]),'name'])
        if j=='utr3' or j=='downstream':
            if tmp_n>=b[tmpchr][ri] and \
            tmp_n<=c.loc[tmpchr+'_'+str(b[tmpchr][ri]),j+'End'] and \
            b[tmpchr][ri]!='NA' and \
            c.loc[tmpchr+'_'+str(b[tmpchr][ri]),j+'End']!='NA':
                d.loc[i,'HIT']=j.upper()+"HIT"
                d.loc[i,'gene']=c.loc[tmpchr+'_'+str(b[tmpchr][ri]),'name']
                df = df.append(d.iloc[[i]])
                print(d.iloc[[i]].to_string(index=False, header=False),
                      j.upper()+"HIT",c.loc[tmpchr+'_'+str(b[tmpchr][ri]),'name'])

outFile = inFile0 + ".append"
try:
    print("Saving data to file "+outFile+" ...")
    df.to_csv(outFile,index=False,sep='\t',na_rep='NA',float_format='%.9g')
except IOError:
    print("Failed to create " + outFile)
    exit(1)
