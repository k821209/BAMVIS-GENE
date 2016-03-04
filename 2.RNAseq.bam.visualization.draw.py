# coding: utf-8
from __future__ import print_function
import pandas as pd
import subprocess
import re
import numpy as np
import sys

genename = sys.argv[3] # genename 'Thhalv10014963m.g.v1.0'
file_bam = sys.argv[2] # input bam file ! should be with index file, run samtools index ok.bam  
#df_gene_gff = pd.read_pickle('../../genemodel_correction/tophat/df.cre.gff.gene.index.pickle')
df_gff_index = pd.read_pickle(sys.argv[1]) # index file 'Esalsugineum_173_v1.0.gene.gff3.pandas.df.pk' 
df_gff_ix = df_gff_index.reset_index().set_index(['genename','longest',2])

df_gff_ix.loc[(genename,'1','mRNA')][[0,3,4]].values
def draw_alignment(genename):
    print ('samtools view %s %s:%d-%d > temp.sam'%(file_bam,chromosome,left,right))
    subprocess.call('samtools view %s %s:%d-%d > temp.sam'%(file_bam,chromosome,left,right),shell=True)
    subprocess.call('cut -f 1-20 temp.sam > temp.sam.cut',shell=True)

    df_sam = pd.read_csv('temp.sam.cut',sep='\t',header=None)

    df_sam['SRRID']  = df_sam[0].apply(lambda x:x.split('.')[0])
    df_sam['READID'] = df_sam[0].apply(lambda x:x.split('.')[1])
    try:
        df_sam['PAIRID'] = df_sam[0].apply(lambda x:x.split('.')[2])
        bPaired = 1
    except IndexError:
        bPaired = 0

    mask = (df_sam[8] > 0)
    df_sam_forward = df_sam[mask]
    df_sam_forward = df_sam_forward.sort(3)
    READID_LIST    = df_sam_forward['READID']
    if bPaired == 1:
        df_sam_index = df_sam.set_index(['SRRID','READID','PAIRID'])
    else: 
        df_sam_index = df_sam.set_index(['SRRID','READID'])
    df_sam_index = df_sam_index.sort([3])

    def cigar_parse(cigar):
        match = re.findall(r'(\d+)(\w)', cigar)
        return match

    # SVG
    rect = '<rect x="%d" y="%d" width="%d" height="%d" style="fill:%s;stroke:%s;stroke-width:1;fill-opacity:0.8;stroke-opacity:0.3" />'
    line = '<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="%s" stroke-width="1" />'
    text = '<text x="%d"  y="%d" style="font-family:Arial;font-size:%d;stroke:black;fill:black;">%s</text>'

    Outfile = open('%s.%s.svg'%(genename,file_bam.split('/')[-1]),'w')

    each_height   = 5
    each_space    = 2
    box_height    = 3
    bridge_height = 2
    genestart     = 3

    total_canvas_rows = 100 + genestart# 2 for gene model

    canvas_width = 1000
    canvas_height = total_canvas_rows*(each_height+each_space)

    real_width = right-left+1+100+100

    print('''<svg height="%d" width="%d">'''%(canvas_height,canvas_width),file=Outfile)
    gene_space = np.zeros([total_canvas_rows,real_width]) # 100 for paired end stretching out of border

    ############################################################################ drawing space set
    def get_ratio(x):
        return int(float(x)/float(real_width) * canvas_width)


    df = df_gff_index.loc[(genename,'1')]
    df = df.reset_index()
    df = df.sort(3)

    fill_CDS   = 'yellow'
    strock_CDS = 'black'

    fill_UTR   = 'blue'
    strock_UTR = 'black'



    for i in df.index:
        if df.loc[i][2] == 'mRNA':
            x1 = get_ratio(df.loc[i][3] - left + 100)
            y1 = int(each_height/2) + genestart + 1
            x2 = get_ratio(df.loc[i][4] - left + 100)
            y2 = int(each_height/2) + genestart + 1
            print(line%(x1,y1,x2,y2,'black'),file=Outfile)
        elif 'UTR' in df.loc[i][2]:
            x1 = get_ratio(df.loc[i][3] - left + 100)
            y1 = genestart + 1 
            w  = get_ratio(df.loc[i][4] - df.loc[i][3])
            h  = each_height
            print(rect%(x1,y1,w,h,fill_UTR,strock_UTR),file=Outfile)
        elif 'CDS' in df.loc[i][2]:
            x1 = get_ratio(df.loc[i][3] - left + 100)
            y1 = genestart + 1
            w  = get_ratio(df.loc[i][4] - df.loc[i][3])
            h  = each_height
            print(rect%(x1,y1,w,h,fill_CDS,strock_CDS),file=Outfile)

    print(text%(15,genestart+each_height*4,10,genename),file=Outfile)

    ########################################################################## gene model drawing


    SRRID_list = set(df_sam_index.index.get_level_values('SRRID'))
    readnumber = 0
    bPass = 1
    for SRRID in SRRID_list:
        READID_list = df_sam_index.loc[SRRID].index.get_level_values('READID')
        for READID in READID_list:
            #print(READID)
            if bPaired == 1:  
                df = df_sam_index.loc[(SRRID,READID)]
                if len(df) == 2:
                    pass
                else: continue
            else:
                if len(df_sam_index.loc[(SRRID,READID)]) < 10: # in the case of multiple match. Just select first one ...
                    df = df_sam_index.loc[(SRRID,READID)].values[0]
                else:
                    df = df_sam_index.loc[(SRRID,READID)].values
            ##################################### cigar check start
            gap_length = 0
            match_length = 0
            if bPaired == 1:
                for cigar in df[5]:
                    cigarM       = cigar_parse(cigar)
                    #print(cigarM)
                    cigarstrings = [x[1] for x in cigarM]
                    cigarvalues  = [x[0] for x in cigarM]
                    for i, cigarstring in enumerate(cigarstrings):
                        if cigarstring == 'N' or cigarstring == 'I':
                            gap_length += int(cigarvalues[i])
                        if cigarstring == 'M':
                            match_length += int(cigarvalues[i])
            else:
                try:
                    cigarM       = cigar_parse(df[5])
                except IndexError:
                    print(df) # Sanity check 
                cigarstrings = [x[1] for x in cigarM]
                cigarvalues  = [x[0] for x in cigarM]
                for i, cigarstring in enumerate(cigarstrings):
                    if cigarstring == 'N' or cigarstring == 'I':
                        gap_length += int(cigarvalues[i])
                    if cigarstring == 'M':
                        match_length += int(cigarvalues[i])

            #print(gap_length,match_length)
            if bPaired == 1:
                fragstart      = min(df[3])-left+100
            else:
                #print (df)
                fragstart    = df[3]-left+100
            if fragstart < 0:
                continue
            if bPaired == 1:
                fragmentsize   = max(df[8])
                fragmentsize_r = min(df[8])
                insertsize     = fragmentsize - gap_length - 200
            if bPaired == 1:
                if insertsize < 10:
                    continue
                if fragmentsize_r > 0 :
                    continue
            else: pass
            ##################################################################### Error remove
            srow = None
            for nrow, row in enumerate(gene_space):
                if nrow <= 3:
                    continue
                if bPaired == 1:
                    if max(row[fragstart:fragstart+fragmentsize]) > 0 :
                        continue
                    else:
                        srow = nrow
                        gene_space[srow,fragstart:fragstart+fragmentsize] += 1
                        break
                else:
                    #print(fragstart)
                    if max(row[fragstart:fragstart+match_length+gap_length]) > 0 :
                        continue
                    else:
                        srow = nrow
                        gene_space[srow,fragstart:fragstart+match_length+gap_length] += 1
                        break

            if srow == None:
                continue
            #print(srow)
            ##################################################################### Row selection 
            readnumber = srow
            #print(readnumber)
            if bPaired == 1:
                for i in df.index:
                    fragsize          = df.loc[i][8]
                    startpoint        = df.loc[i][3]
                    secondstartpoint  = df.loc[i][7]
                    cigar        = df.loc[i][5]
                    cigarM       = cigar_parse(cigar)
                    #print(cigarM)
                    cigarstrings = [x[1] for x in cigarM]
                    cigarvalues  = [x[0] for x in cigarM] 
                    ##################################################################
                    coverleng = 0
                    for j, cigarstring in enumerate(cigarstrings):
                        if cigarstring == 'M' or cigarstring == 'N' or cigarstring == 'I':
                            coverleng += int(cigarvalues[j])
                    if fragsize > 0:
                        x1 = get_ratio(startpoint - left + 100 + coverleng)
                        y1 = readnumber * (each_height + each_space) + int(each_height/2)

                        y2 = readnumber * (each_height + each_space) + int(each_height/2) + int(each_height/2)
                        x3 = get_ratio(secondstartpoint - left + 100)
                        x2 = (x3 + x1)/2
                        print(line%(x1,y1,x2,y2,'black'),file=Outfile)
                        print(line%(x2,y2,x3,y1,'black'),file=Outfile)

                    ################################################################# paired end line
                    #print(i)
                    w = 0
                    #print (cigarstrings)
                    #print (cigarvalues)
                    for j, cigarstring in enumerate(cigarstrings):
                        if cigarstring == 'M':
                            x = get_ratio(w + startpoint - left + 100)
                            #print(x)
                            if x < 0:
                                continue
                            y = readnumber * (each_height + each_space)

                            h = int(each_height)
                            ew = int(cigarvalues[j]) # each width
                            w += ew
                            #print(x,y,w,h)
                            print(rect%(x,y,get_ratio(ew),h,'pink','black'),file=Outfile)

                        elif cigarstring == 'N' or cigarstring == 'I':
                            x1 = get_ratio(w + startpoint - left + 100)
                            if x1 < 0:
                                continue
                            y1 = readnumber * (each_height + each_space) + int(each_height/2)
                            ew = int(cigarvalues[j]) # each width
                            w += ew
                            x2 = get_ratio(w + startpoint - left + 100)

                            print(line%(x1,y1,x2,y1,'red'),file=Outfile)
            else:
                fragsize          = df[8]
                startpoint        = df[3]
                secondstartpoint  = df[7]
                #print(df)
                cigar             = df[5]
                cigarM            = cigar_parse(cigar)
                #print(cigarM)
                cigarstrings      = [x[1] for x in cigarM]
                cigarvalues       = [x[0] for x in cigarM]
                w = 0 
                for j, cigarstring in enumerate(cigarstrings):
                    if cigarstring == 'M':
                        x = get_ratio(w + startpoint - left + 100)
                        #print(x)
                        if x < 0:
                            continue
                        y = readnumber * (each_height + each_space)

                        h = int(each_height)
                        ew = int(cigarvalues[j]) # each width
                        w += ew
                        #print(x,y,w,h)
                        print(rect%(x,y,get_ratio(ew),h,'pink','black'),file=Outfile)

                    elif cigarstring == 'N' or cigarstring == 'I':
                        x1 = get_ratio(w + startpoint - left + 100)
                        if x1 < 0:
                            continue
                        y1 = readnumber * (each_height + each_space) + int(each_height/2)
                        ew = int(cigarvalues[j]) # each width
                        w += ew
                        x2 = get_ratio(w + startpoint - left + 100)

                        print(line%(x1,y1,x2,y1,'red'),file=Outfile)


    print('''</svg>
    ''',file=Outfile)

    Outfile.close()               
chromosome, left, right = df_gff_ix.loc[(genename,'1','mRNA')][[0,3,4]].values[0]
print(genename, chromosome, left, right) 
draw_alignment(genename)
#get_ipython().magic(u'save "RNAseq.bam.visualization.draw.py" _ih[193:200]')
