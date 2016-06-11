from __future__ import print_function
import numpy as np
## retrieve blocks with continuously have values that are larger than depth_cut in array.
def get_block(array,depth_cut=10,lim_len_block=100):
    block_list = []
    #print(len(np.shape(array)))
    if len(np.shape(array)) == 1:
        rows = 1
        block = []
        for n,j in enumerate(array):
            if j > depth_cut:
                block.append(n)
            else:
                if len(block) > lim_len_block:
                    block_list.append([block[0],block[-1]])
                    block = []
                else:
                    block = []
        if len(block) > 0:
            block_list.append([block[0],block[-1]])
    else: 
        rows, columns = np.shape(array)       
        for i in range(rows):
            earray = array[i]
            block = []
            for n,j in enumerate(earray):
                if j > depth_cut:
                    block.append(n)
                else:
                    if len(block) > lim_len_block:
                        block_list.append([i,block[0],block[-1]])
                        block = []
                    else:
                        block = []
        if len(block) > 0:
            block_list.append([block[0],block[-1]])
    return block_list

#0.read paired
#1.read mapped in proper pair
#2.read unmapped
#3.mate unmapped
#4.read reverse strand
#5.mate reverse strand
#6.first in pair
#7.second in pair
#8.not primary alignment
#9.read fails platform/vendor quality checks
#10.read is PCR or optical duplicate
#11.supplementary alignment


def flagparser(intFlag):
    dic_key = ['supplementary alignment','read is PCR or optical duplicate','read fails platform/vendor quality checks',\
               'not primary alignment','second in pair','first in pair','mate reverse strand','read reverse strand','mate unmapped',\
              'read unmapped','read mapped in proper pair','read paired']
    bFlag     = "{0:b}".format(intFlag)
    addzero   = '0'*(12-len(bFlag))
    dic_value = list(addzero+bFlag)
    #print (addzero+bFlag)
    return(dict(zip(dic_key,map(int,dic_value))))

def list2txt(outfilename,inlist):
    Outfile = open(outfilename,'w')
    for each in inlist:
        print(each,file=Outfile)
    Outfile.close()


gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}
def rev_comp(strSeq):
	dicComp = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
	strCseq = ''
	for i in strSeq:
		try:
			strCseq += dicComp[i.upper()]
		except KeyError:
			strCseq += 'N'
	# End of for i
	return(strCseq[::-1])
def translation(strSeq):
	strPep = ''
	for i in range(int(len(strSeq)/3)):
		try:
			strPep += gencode[strSeq[i*3:i*3+3].upper()]
		except KeyError:
			strPep += 'X'
	# End of for i
	return(strPep)
def fasta2dic(file_fasta, dic):
	bulk      = open(file_fasta).read()
	bulk_list = bulk.split('>')
	for each in bulk_list:
		if each == '':
			continue
		strHD  = each.split('\n')[0].split()[0]
		strSeq = ''.join(each.split('\n')[1:])
		dic[strHD] = strSeq
	return(dic)
def Fasta2dic(file_fasta):
	dic       = {}
	bulk      = open(file_fasta).read()
	bulk_list = bulk.split('>')
	for each in bulk_list:
		if each == '':
			continue
		strHD  = each.split('\n')[0].split()[0]
		strSeq = ''.join(each.split('\n')[1:])
		dic[strHD] = strSeq
	return(dic)
def Fasta2dic_all(file_fasta):
	dic       = {}
	bulk      = open(file_fasta).read()
	bulk_list = bulk.split('>')
	for each in bulk_list:
		if each == '':
			continue
		strHD  = each.split('\n')[0]
		strSeq = ''.join(each.split('\n')[1:])
		dic[strHD] = strSeq
	return(dic)

def dic2fa(dic,filename):
	Outfile = open(filename,'w')
	for key in dic:
		print('>'+key,file=Outfile)
		print(dic[key],file=Outfile)
	Outfile.close()

		



