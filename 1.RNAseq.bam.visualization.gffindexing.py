# coding: utf-8
import sys
import pandas as pd


file_gff = sys.argv[1] # gff file '/ref/Cre/DroughtNet/PhytozomeV10_download/References/Alyrata/annotation/Alyrata_107_v1.0.gene_exons.gff3'
df_gff   = pd.read_csv(file_gff,sep='\t',comment='#',header=None)

dicT2G      = {} # transcript to genename
df_gff_mRNA = df_gff[(df_gff[2] == 'mRNA')]

def infoparse(x):
    key   = [i.split('=')[0] for i in x.split(';')]
    value = [i.split('=')[1] for i in x.split(';')]
    return dict(zip(key,value))

def get_genename(x):
    if x[2] == 'gene':
        return None
    key   = [y.split('=')[0] for y in x[8].split(';')]
    value = [y.split('=')[1] for y in x[8].split(';')]
    dic   = dict(zip(key,value))
    return dic['Parent']
def get_transcript(x):
    key   = [y.split('=')[0] for y in x[8].split(';')]
    value = [y.split('=')[1] for y in x[8].split(';')]
    dic   = dict(zip(key,value))
    if x[2] == 'gene':
        return None
    elif x[2] == 'mRNA':
        return dic['ID']
    else: 
        return dic['Parent']
df_gff_mRNA['genename']   = df_gff_mRNA.apply(get_genename,axis=1)
df_gff_mRNA['transcript'] = df_gff_mRNA.apply(get_transcript,axis=1)
df_gff_mRNA_ix            = df_gff_mRNA.set_index('transcript')
df_gff['transcriptname']  = df_gff.apply(get_transcript,axis=1)



def get_genename2(x):
    if x != None:
        return df_gff_mRNA_ix.loc[x]['genename']
    else:
        return None
df_gff['genename'] = df_gff['transcriptname'].apply(get_genename2)

# grep longest transcript names
df_gff_mRNA_ix['longest'] = df_gff_mRNA_ix[8].apply(lambda x : infoparse(x)['longest'])

def get_longest(x):
    if x == None:
        return None
    try:
        return df_gff_mRNA_ix.loc[x]['longest']
    except KeyError:
        return None
df_gff['longest'] = df_gff['transcriptname'].apply(get_longest)
df_gff_index      = df_gff.set_index(['genename','longest'])
df_gff_index.to_pickle(file_gff+'.pandas.df.pk')

