import pandas as pd

a = pd.read_csv("/home/i/Desktop/git/public/pySNPAAMapper/ChrAll_knownGene.txt.exons",delimiter="\t")

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

d = pd.read_csv("/home/i/Desktop/git/public/pySNPAAMapper/007_crop.vcf",
                delimiter="\t",header=26)  #header numbers start from 0

intronOpt=-1
exonbuffer=6

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
                print(d.iloc[[i]].to_string(index=False, header=False),
                      "CDSHIT",c.loc[tmpchr+'_'+str(b[tmpchr][ri]),'name'])
        if j=='intron':
            if intronOpt==1 and tmp_n>=b[tmpchr][ri] and \
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
                print(d.iloc[[i]].to_string(index=False, header=False),
                      "INTRONHIT."+str(hitbuffer)+"."+direction,
                      c.loc[tmpchr+'_'+str(b[tmpchr][ri]),'name'])
            if intronOpt==-1 and tmp_n>=b[tmpchr][ri] and \
            tmp_n<=(c.loc[tmpchr+'_'+str(b[tmpchr][ri]),j+'End']+1) and \
            b[tmpchr][ri]!='NA' and \
            c.loc[tmpchr+'_'+str(b[tmpchr][ri]),j+'End']!='NA':
                print(d.iloc[[i]].to_string(index=False, header=False),
                      "INTRONHIT",
                      c.loc[tmpchr+'_'+str(b[tmpchr][ri]),'name'])
        if j=='utr5' or j=='upstream':
            if tmp_n>b[tmpchr][ri] and \
            tmp_n<=(c.loc[tmpchr+'_'+str(b[tmpchr][ri]),j+'End']+1) and \
            b[tmpchr][ri]!='NA' and \
            c.loc[tmpchr+'_'+str(b[tmpchr][ri]),j+'End']!='NA':
                print(d.iloc[[i]].to_string(index=False, header=False),
                      j.upper()+"HIT",c.loc[tmpchr+'_'+str(b[tmpchr][ri]),'name'])
        if j=='utr3' or j=='downstream':
            if tmp_n>=b[tmpchr][ri] and \
            tmp_n<=c.loc[tmpchr+'_'+str(b[tmpchr][ri]),j+'End'] and \
            b[tmpchr][ri]!='NA' and \
            c.loc[tmpchr+'_'+str(b[tmpchr][ri]),j+'End']!='NA':
                print(d.iloc[[i]].to_string(index=False, header=False),
                      j.upper()+"HIT",c.loc[tmpchr+'_'+str(b[tmpchr][ri]),'name'])
