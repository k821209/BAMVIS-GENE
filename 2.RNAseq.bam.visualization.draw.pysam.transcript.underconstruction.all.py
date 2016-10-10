#!/usr/bin/python

from __future__ import print_function
import sys
sys.path.append('/mnt/c/ubuntu.download/pipelines')
import pandas as pd
import subprocess
import re
import numpy as np
import kang
import pysam
from tqdm import tqdm
import glob
#+ preparing references
df_gff_index = pd.read_pickle('/mnt/c/ubuntu.download/pipelines/pandas_df/Creinhardtii_281_v5.5.gene.gff3.pandas.df.pk')
file_fa      = './Creinhardtii_281_v5.0.fa'
dic_fa       = kang.Fasta2dic(file_fa)
df_gff_ix    = df_gff_index.reset_index().set_index(['transcriptname',2])
df_gff_ix.sortlevel(inplace=True)

#- preparing references done


#+ global vars.
genelist_pre = set(df_gff_ix.reset_index()['transcriptname'].dropna())
bamlist                 = ['./all.merged.bam']

donelist = set([x.replace('.bamvisgene.svg','').replace('./','') for x in glob.glob('./*.svg')])
genelist = list(genelist_pre-donelist)
print("%d total transcript, %d transcripts done, %d genes todo"%(len(genelist_pre),len(donelist),len(genelist)))
print(genelist[0])
print(list(donelist)[0])
#+ svg configure 
rect     = '<rect x="%d" y="%d" width="%d" height="%d" style="fill:%s;stroke:%s;stroke-width:1;fill-opacity:1;stroke-opacity:0.3" />'
svg_line = '<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="%s" stroke-width="%d" />'
text     = '<text x="%d"  y="%d" style="font-family:Arial;font-size:%d;stroke:black;fill:black;">%s</text>'
#- svg configure done 

#+ color conf : haruhi theme 
c_utr     = '#f50002'
c_cds     = '#53a1b5'
c_fstrand = '#5d3c2d'
c_sstrand = '#f5b024'
#- color conf done

#+ defs 
def init_svg(f):  
   print('''<?xml version="1.0" encoding="utf-8" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"
  "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<svg height="%s" width="%s" version="1.1" viewBox="0 0 %s %s" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">'''%(canvas_height,canvas_width,canvas_width,canvas_height),file=f)  
    
def end_svg(f):
    print('''</svg>
    ''',file=f)

def get_ratio(x):
    return int(canvas_width * float(x)/float(real_width))
def get_ratio_f(x):
    return (canvas_width * float(x)/float(real_width))
def draw_words(y,strWord,f):
    print(text%(15,y+8,10,strWord),file=f)
    return y+8+each_space+5
def draw_gene(start_height,left,right,genename,f):
    Outfile = f
    df = df_gff_ix.loc[(genename)]
    df = df.reset_index()
    df = df.sort(3)

    fill_CDS   = c_cds
    strock_CDS = 'black'

    fill_UTR   = c_utr
    strock_UTR = 'black'

    for i in df.index:
        if df.loc[i][2] == 'mRNA':
            strand  = df.loc[i][6]
            mx1     = get_ratio(df.loc[i][3] - left + 100)
            my1     = start_height + int(each_height/2) + genestart + 1
            mx2     = get_ratio(df.loc[i][4] - left + 100)
            my2     = start_height + int(each_height/2) + genestart + 1
            print(svg_line%(mx1,my1,mx2,my2,'black',1),file=Outfile)
    for i in df.index:
        if 'UTR' in df.loc[i][2]:
            ux1 = get_ratio(df.loc[i][3] - left + 100)
            uy1 = start_height + genestart + 1 
            uw  = get_ratio(df.loc[i][4] - df.loc[i][3])
            uh  = each_height
            print(rect%(ux1,uy1,uw,uh,fill_UTR,strock_UTR),file=Outfile)
        elif 'CDS' in df.loc[i][2]:
            cx1 = get_ratio(df.loc[i][3] - left + 100)
            cy1 = start_height + genestart + 1
            cw  = get_ratio(df.loc[i][4] - df.loc[i][3])
            ch  = each_height
            print(rect%(cx1,cy1,cw,ch,fill_CDS,strock_CDS),file=Outfile)
    print(text%(mx1,genestart+each_height*4,10,'%s (%s) %s:%d-%d'%(genename,strand,chromosome,left,right)),file=Outfile)
    return start_height + genestart+each_height*4
def match_count(a,b):
        count = 0
        if len(a)!=len(b):
            return None 
        else:
            for n,i in enumerate(a):
                if a[n] == b[n]:
                    count += 1
                if a[n] != b[n]:
                    break 
        return count
def match_list(primer):
    i = 0
    j = 0 
    match_list = []
    primer_seq = primer 
    while j < len(geneseq):
        i = 0
        inseq = geneseq[j:]
        mc_list = []
        while i != len(inseq):
            mc = match_count(inseq[i:i+len(primer_seq)],primer_seq)
            mc_list.append([i,mc])
            i += 1
        loc,mc = max(mc_list,key=lambda x:x[1])
        match_list.append([j+loc,j+loc+mc])
        j = j + loc + mc
        
        primer_seq = primer_seq[mc:]
        if primer_seq == '':
            break
    return match_list
def draw_primer(start_height,primerlist,strand):
    for strF,strR in primerlist:
        if strand == '+':
            strR = kang.rev_comp(strR)        
        else:
            strF = kang.rev_comp(strF)
        x_list = []
        for x1,x2 in match_list(strF):
            x_list.append(x1)
            x_list.append(x2)
            y1 = start_height + int(each_height/2)
            y2 = start_height + int(each_height/2)
            print(svg_line%(get_ratio(x1+100),y1,get_ratio(x2+100),y2,'blue',4),file=Outfile)
        x_list.sort()
        print(svg_line%(get_ratio(x_list[0]+100),y1,get_ratio(x_list[-1]+100),y2,'blue',1),file=Outfile)
        x_list_2 = []
        for x3,x4 in match_list(strR):
            x_list_2.append(x3)
            x_list_2.append(x4)
            y1 = start_height + int(each_height/2)
            y2 = start_height + int(each_height/2)
            print(svg_line%(get_ratio(x3+100),y1,get_ratio(x4+100),y2,'red',4),file=Outfile)
        x_list_2.sort()
        print(svg_line%(get_ratio(x_list_2[0]+100),y1,get_ratio(x_list_2[-1]+100),y2,'red',1),file=Outfile)
        if strand == '+':
            print(svg_line%(get_ratio(x_list[-1]+100),y1,get_ratio(x_list_2[0]+100),y2,'gray',1),file=Outfile)
        else: print(svg_line%(get_ratio(x_list_2[-1]+100),y1,get_ratio(x_list[0]+100),y2,'gray',1),file=Outfile)
    return y2 + each_height

def cigar_parse(cigar):
    match = re.findall(r'(\d+)(\w)', cigar)
    return match


# pysam is problematic because it cannot handle efficiently the mate read information. 
# maybe it is better to have table in advance using pysam 


def draw_alignment(start_height,chromosome,left,right,file_bam,f,total_canvas_rows=5):
    Outfile = f
    real_width        = right-left+1+100+100
    gene_space        = np.zeros([total_canvas_rows,real_width])
    samfile           = pysam.Samfile( file_bam, "rb" )
    it                = samfile.fetch(chromosome,int(left),int(right))
    dicFID2row = {} # fragment id : line.query_name
    dicFID2count = {}
    dicFID2line = {}
    for line in tqdm(it):
        if line.is_paired == True and line.is_proper_pair == False:
            continue
        if line.is_duplicate   == True:
            continue
        if line.is_qcfail      == True:
            continue
        if line.is_secondary   == True:
            continue
        try:
            dicFID2line[line.query_name].append(line)
        except KeyError:
            dicFID2line[line.query_name] = [line]
    dicFID2line_list = dicFID2line.keys()
    dicFID2line_list.sort(key=lambda x : int(dicFID2line[x][0].reference_start))
    for fragmentid in tqdm(dicFID2line_list):
        lines = dicFID2line[fragmentid]
        lines.sort(key=lambda x : int(x.reference_start))
        line  = lines[0] 
        fragstart      = line.reference_start-left+1+100
        if fragstart < 0 :
            continue
        cigar          = line.cigarstring
        try:
            cigarM       = cigar_parse(cigar)
        except TypeError:
            continue
        cigarstrings  = [x[1] for x in cigarM]
        cigarvalues   = [x[0] for x in cigarM]

        if line.is_paired == True:
            if line.is_proper_pair == False:
                continue
            if len(lines) < 2:
                continue
            try: 
                fragmentsize   = line.mpos+lines[1].reference_length - line.reference_start # reference length of fragment (intron added)
            except TypeError:
                continue 
            bridge_length  = line.mpos - line.reference_start - line.reference_length
            if fragmentsize > 0 and bridge_length < 10 :
                continue   
        else:
            fragmentsize   = line.reference_length

        if abs(fragmentsize) == 0 :
            continue
        if fragmentsize < 0 :
            print(fragmentid,len(lines),cigar,line.mpos,line.reference_start,line.reference_length)
            print(1111)
            continue
        
        #+ row selection
        try:
            srow = dicFID2row[fragmentid]
        except KeyError:
            srow = None
            for nrow, row in enumerate(gene_space):
                if 1:
                    if max(row[fragstart:fragstart+fragmentsize]) > 0 :
                        continue
                    else:
                        srow = nrow
                        gene_space[srow,fragstart:fragstart+fragmentsize] += 1
                        break
            dicFID2row[fragmentid] = srow

        if srow == None:
            continue
        #- row selection


        if line.is_paired == True:
            if fragstart < 0:
                continue
            #keepdoing      = 1
            #for i in np.arange(fragstart,fragstart+fragmentsize):
            #    try:
            #        flatten = np.sum(gene_space[:,i])
            #    except IndexError:
            #        keepdoing = 0
            #    if flatten == total_canvas_rows:
            #        keepdoing = 0
            #if keepdoing == 0 :
            #    continue
            if line.is_read1 == True and line.is_reverse == False:
                blockcolor = c_sstrand
            elif line.is_read2 == True and line.is_reverse == True:
                blockcolor = c_sstrand
            elif line.is_read1 == True and line.is_reverse == True:
                blockcolor = c_fstrand
            elif line.is_read2 == True and line.is_reverse == False:
                blockcolor = c_fstrand

            #+ draw bridge for PE read 1
            if 1:
                coverleng = 0
                for j, cigarstring in enumerate(cigarstrings):
                    if cigarstring == 'M' or cigarstring == 'N' or cigarstring == 'I':
                        coverleng += int(cigarvalues[j])
                if 1:
                    x1 = get_ratio(fragstart + coverleng)
                    y1 = start_height + srow * (each_height + each_space) + int(each_height/2)
                    y2 = start_height + srow * (each_height + each_space) + int(each_height/2) + int(each_height/2)
                    x3 = get_ratio(line.mpos - left+1 + 100)
                    x2 = (x3 + x1)/2
                    print(svg_line%(x1,y1,x2,y2,'black',1),file=Outfile)
                    print(svg_line%(x2,y2,x3,y1,'black',1),file=Outfile)
                    #print(text%(x1,y1,5,fragmentid),file=Outfile)
            #- draw bridge for PE

        ##################################### cigar check start

            #+ draw match and gap
            for eline in lines:
                w = 0
                fragstart = eline.reference_start - left + 1 + 100
                cigarstrings = eline.cigarstring
                cigar        = eline.cigarstring
                try:
                    cigarM       = cigar_parse(cigar)
                except TypeError:
                    continue
                cigarstrings  = [x[1] for x in cigarM]
                cigarvalues   = [x[0] for x in cigarM]
 
                for j, cigarstring in enumerate(cigarstrings):
                    if cigarstring == 'M':
                        x = get_ratio(w + fragstart) # get ratio have gap with actual value
                        if x < 0:
                            continue
                        y = start_height + srow * (each_height + each_space)
                        h = int(each_height)
                        ew = int(cigarvalues[j]) # each width
                        w += ew
                        print(rect%(x,y,get_ratio(ew),h,blockcolor,'black'),file=Outfile)
                    elif cigarstring == 'N' or cigarstring == 'I':
                        x1 = get_ratio(w + fragstart)
                        if x1 < 0:
                            continue
                        y1 = start_height + srow * (each_height + each_space) + int(each_height/2)
                        ew = int(cigarvalues[j]) # each width
                        w += ew
                        x2 = get_ratio(w + fragstart)
                        print(svg_line%(x1,y1,x2,y1,'pink',1),file=Outfile)
            #- draw match and gap 



        else:
            fragstart      = line.reference_start-left+1+100
            if fragstart < 0:
                continue

            if line.is_reverse == True:
                    blockcolor = c_sstrand
            else:
                    blockcolor = c_fstrand
            w = 0
            for j, cigarstring in enumerate(cigarstrings):
                if cigarstring == 'M':
                    x = get_ratio(w + fragstart)
                    if x < 0:
                        continue
                    y = start_height + srow * (each_height + each_space)
                    h = int(each_height)
                    ew = int(cigarvalues[j]) # each width
                    w += ew
                    print(rect%(x,y,get_ratio(ew),h,blockcolor,'black'),file=Outfile)
                elif cigarstring == 'N' or cigarstring == 'I':
                    x1 = get_ratio(w + fragstart)
                    if x1 < 0:
                        continue
                    y1 = start_height + srow * (each_height + each_space) + int(each_height/2)
                    ew = int(cigarvalues[j]) # each width
                    w += ew
                    x2 = get_ratio(w + fragstart)
                    print(svg_line%(x1,y1,x2,y1,'pink',1),file=Outfile)
    return start_height+(each_height + each_space) * len(gene_space)



for genename in genelist:
    outfilename             = genename+'.bamvisgene.svg'
    strand                  = df_gff_ix.loc[(genename)][6].values[0]
    chromosome, left, right = df_gff_ix.loc[(genename,'mRNA')][[0,3,4]].values[0]
    #+ alignment vars.
    each_height   = 5
    each_space    = 2
    box_height    = 3
    bridge_height = 2
    genestart     = 3
    total_canvas_rows = 10 # 2 for gene model 
    canvas_width      = 1000
    canvas_height     = (total_canvas_rows + 10)*(each_height+each_space)*len(bamlist) + 1000 # 1000 : just more blank 
    real_width        = right-left+1+100+100
    #- alignments vars. done
    #eneseq                 = dic_fa[chromosome][left-1:right]
    
    with open(outfilename,'w') as f:
        init_svg(f)
        endp = draw_gene(5,left,right,genename,f)
        for bam in bamlist:
            endp = draw_alignment(endp,chromosome,left,right,bam,f,20)
            endp = draw_words(endp,bam.split('/')[-1],f)
            print(1)
        end_svg(f)
