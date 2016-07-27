'''f
Created on 2012. 9. 18.

@author: Akard3
'''

import string
import os
import sys
import struct
import re
import time
import fileinput








def baseconvert(n, base):
	"""convert positive decimal integer n to equivalent in another base (2-36)"""
	digits = "0123456789abcdefghijklmnopqrstuvwxyz"
	try:
		n = int(n)
		base = int(base)
	except:
		return ""


	if n < 0 or base < 2 or base > 36:
		return ""

	s = ""
	while 1:
		r = n % base
		s = digits[r] + s
		n = n / base
		if n == 0:
			break

	return s

def    formatChrNameForCElegans(reference_name_str):
	if reference_name_str == "I":
		return "chr1"
	elif reference_name_str == "II":
		return "chr2"
	elif reference_name_str == "III":
		return "chr3"
	elif reference_name_str == "IV":
		return "chr4"
	elif reference_name_str == "V":
		return "chr5"
	elif reference_name_str == "X":
		return "chrX"
	else:
		return -1

def    formatChrNameForHuman(reference_name_str):
	if reference_name_str == "1":
		return "chr1"
	elif reference_name_str == "2":
		return "chr2"
	elif reference_name_str == "3":
		return "chr3"
	elif reference_name_str == "4":
		return "chr4"
	elif reference_name_str == "5":
		return "chr5"
	elif reference_name_str == "6":
		return "chr6"
	elif reference_name_str == "7":
		return "chr7"
	elif reference_name_str == "8":
		return "chr8"
	elif reference_name_str == "9":
		return "chr9"
	elif reference_name_str == "10":
		return "chr10"
	elif reference_name_str == "11":
		return "chr11"
	elif reference_name_str == "12":
		return "chr12"
	elif reference_name_str == "13":
		return "chr13"
	elif reference_name_str == "14":
		return "chr14"
	elif reference_name_str == "15":
		return "chr15"
	elif reference_name_str == "16":
		return "chr16"
	elif reference_name_str == "17":
		return "chr17"
	elif reference_name_str == "18":
		return "chr18"
	elif reference_name_str == "19":
		return "chr19"
	elif reference_name_str == "20":
		return "chr20"
	elif reference_name_str == "21":
		return "chr21"
	elif reference_name_str == "22":
		return "chr22"
	elif reference_name_str == "23":
		return "chr23"
	elif reference_name_str == "X":
		return "chrX"
	elif reference_name_str == "Y":
		return "chrY"
	else:
		return -1


def checkParentID(part1, part2):
	if part1[0] == -1 or part2[0] == -1:
		return False
	tmp1 = part1[-1].split(";")[0].split("=")[1]
	tmp2 = part2[-1].split(";")[0].split("=")[1]
	if tmp1 == tmp2:
		return True
	else:
		return False

'''

def addData(data,strand,start_first,start_second,end_first,end_second,chr_name,frame,Parent):
    if strand == '+':
        if data.has_key(end_first):
            check = 0
            for i in range(len(data[end_first])):
                if data[end_first][i][0] == start_second:
                    temp = data[end_first][i]
                    temp2 = temp[1]
                    if temp2.has_key(SAM_file_name):
                        data[end_first][i][1][SAM_file_name] = data[end_first][i][1].get(SAM_file_name)+1 
                    else:
                        data[end_first][i][1][SAM_file_name] = 1
                    if start_first < temp[3]:
                        temp[3] = start_first
                    if end_second > temp[4]:
                        temp[4] = end_second
                    data[end_first][i] = temp
                    check = 1
            if check == 0:
                data[end_first].append([start_second,{SAM_file_name:1},strand,start_first,end_second,chr_name,frame,Parent])
        else:
            data[end_first] = []
            data[end_first].append([start_second,{SAM_file_name:1},strand,start_first,end_second,chr_name,frame,Parent])
    else:
        if data.has_key(start_first):
            check = 0
            for i in range(len(data[start_first])):
                if data[start_first][i][0] == end_second:
                    temp = data[start_first][i]
                    temp2 = temp[1]
                    if temp2.has_key(SAM_file_name):
                        data[start_first][i][1][SAM_file_name] = data[start_first][i][1].get(SAM_file_name)+1 
                    else:
                        data[start_first][i][1][SAM_file_name] = 1
                    if end_first > temp[3]:
                        temp[3] = end_first
                    if start_second < temp[4]:
                        temp[4] = start_second
                    data[start_first][i] = temp
                    check = 1
            if check == 0:
                data[start_first].append([end_second,{SAM_file_name:1},strand,end_first,start_second,chr_name,frame,Parent])
        else:
            data[start_first] = []
            data[start_first].append([end_second,{SAM_file_name:1},strand,end_first,start_second,chr_name,frame,Parent])
    return data
'''
def decomp(CIGAR,CIGAR_begin,sign,strand):
    del_location = []
    seq_location = []
    pattern = re.compile('S|N|M|D|I|[0-9]+')
    CIGAR_info = [(m.group()) for m in pattern.finditer(CIGAR)]
    CIGAR_interval = [int(n) for i,n in enumerate(CIGAR_info) if i%2 == 0]
    index_D = [i for i,n in enumerate(CIGAR_info) if n == sign]
    for i in index_D:
        del_interval = int(CIGAR_info[i-1])
        del_start = CIGAR_begin
        for j in range((i-1)/2):
            if CIGAR_info[2*j+1] != 'I' and  CIGAR_info[2*j+1] != 'S':
                del_start += CIGAR_interval[j]
        #del_start = CIGAR_begin + sum(CIGAR_interval[0:(i-1)/2])
        del_end = del_start + del_interval
        del_location.append((del_start,del_end))
        seq_finder = 0
        for j in range((i-1)/2):
            if CIGAR_info[2*j+1] != 'D' and CIGAR_info[2*j+1] != 'N':
                seq_finder += CIGAR_interval[j]
        seq_location.append((seq_finder,seq_finder+CIGAR_interval[(i-1)/2]))
    return    [del_location,seq_location]



def addBAM(ID):
	id = 1
	for line in sys.stdin: #SourceFile:
	    #print line
	    if line == "":
	        continue
	    if line.startswith(">"):
	    	continue
	    if line.startswith("@"):
	        continue
	    if line.startswith("track"):
	        continue
	    data = line.strip()
	    #data = data.split('\t')
	    data = data.split()
	    #if not data[2].startswith('chr'):
	    #	continue
	    
	    QNAME = data[0]
	    FLAG = int(data[1])
	    RNAME = data[2]
	    POS = int(data[3])-1
	    MAPQ = int(data[4])
	    CIGAR = data[5]
	    RNEXT = data[6]
	    PNEXT = int(data[7])
	    TLEN = int(data[8])
	    SEQ = data[9]
	    QUAL = data[10]
	    p = re.compile('[0-9]+M')
	    ###############################
	    #if RNAME != chr:
	    #	continue
	    ###############################
	    if check_sample[0] == 1:
	    	if 'N' in CIGAR:
	    	 	id += 1
	    		if id > check_sample[1]:
	    	  		break
	    
	    
	    tmp_flag = baseconvert(FLAG,2)
	    if tmp_flag == "0" or tmp_flag == 0:
	        strand = "+"
	    elif len(tmp_flag)<5:
	    	strand = "+"
	    elif int(tmp_flag[-5]) == 0:
	        strand = "+"
	    elif int(tmp_flag[-5]) == 1:
	        strand = "-"
	    else:
	        print QNAME, CIGAR, POS, TLEN, "can't parse strand from: ", FLAG
	        continue
	    #chr_name = formatChrNameForCElegans(RNAME)
	    if Species == 'Human':
	    	chr_name = formatChrNameForHuman(RNAME)
	    	#chr_name = RNAME
	    	
	    elif Species == 'CElegans':
    		chr_name = formatChrNameForCElegans(RNAME)
	    
	    if chr_name not in chr_ref:
	    	continue
	    #chr_name = formatChrNameForHuman(RNAME)
	    if chr_name == -1:
	    	continue
	    #############################################
	    CIGAR_begin = POS
	    #############################################
	    if 'S' in CIGAR:
	    	continue

	    '''
	    if 'D' in CIGAR:
	    	sign = 'D'
	    	[del_location,dummy] = decomp(CIGAR,CIGAR_begin,sign,strand)
	    	del_temp = del_data.get(chr_name)
	    	for coord_info in del_location:
	    		if del_temp.has_key(coord_info[0]):
	    			check = 0
	    			for i,check_end in enumerate(del_temp[coord_info[0]]):
	    				if check_end[0] == coord_info[1] and strand == check_end[2]:
	    					if check_end[1].has_key(SAM_file_name):
	    						del_temp[coord_info[0]][i][1][SAM_file_name] += 1
    						else:
    							del_temp[coord_info[0]][i][1][SAM_file_name] = 1
	    					check = 1
	    			if check == 0:
	    				del_temp[coord_info[0]].append([coord_info[1],{SAM_file_name:1},strand])
	    		else:
	    			del_temp[coord_info[0]] = [[coord_info[1],{SAM_file_name:1},strand]]
	    	del_data[chr_name] = del_temp		
	    	#for i in del_location:
	    	# i[0] is starting position of del
	    	# i[1] is ending position of del
	    	# make a database for a del and compare the del position if they are overlap
	    	
	    	'''
	    if 'N' in CIGAR:
	    	sign = 'N'
	    	[data_location, seq_location] = decomp(CIGAR,CIGAR_begin,sign,strand)
	    	data = full_data.get(chr_name)
	    	for coord_info in data_location:
	    		if data.has_key(coord_info[0]):
	    			check = 0
	    			for i,check_end in enumerate(data[coord_info[0]]):
	    				if check_end[0] == coord_info[1] and strand == check_end[2]:
	    					if check_end[1].has_key(SAM_file_name):
	    						data[coord_info[0]][i][1][SAM_file_name] += 1
    						else:
    							data[coord_info[0]][i][1][SAM_file_name] = 1
	    					#data[coord_info[0]][i][1] += 1
	    					check = 1
	    			if check == 0:
	    				data[coord_info[0]].append([coord_info[1],{SAM_file_name:1},strand])
    			else:
    				data[coord_info[0]] = [[coord_info[1],{SAM_file_name:1},strand]]
    		full_data[chr_name] = data
    	
















SAM_file_name_In_String = sys.argv[1]
Output_TMP_file_name = sys.argv[2]

Species = 'Human' #'CElegans', 'Human'
check_sample = [0,1000] #0: check off, 1: check on


full_data = {'chr1':{},'chr2':{},'chr3':{},'chr4':{},'chr5':{},'chr6':{},'chr7':{},'chr8':{},'chr9':{},'chr10':{},'chr11':{},'chr12':{},'chr13':{},'chr14':{},'chr15':{},'chr16':{},'chr17':{},'chr18':{},'chr19':{},'chr20':{},'chr21':{},'chr22':{},'chr23':{},'chrX':{},'chrY':{}}
file_name_info = {}
in_data = {'chr1':{},'chr2':{},'chr3':{},'chr4':{},'chr5':{},'chr6':{},'chr7':{},'chr8':{},'chr9':{},'chr10':{},'chr11':{},'chr12':{},'chr13':{},'chr14':{},'chr15':{},'chr16':{},'chr17':{},'chr18':{},'chr19':{},'chr20':{},'chr21':{},'chr22':{},'chr23':{},'chrX':{},'chrY':{}}
del_data = {'chr1':{},'chr2':{},'chr3':{},'chr4':{},'chr5':{},'chr6':{},'chr7':{},'chr8':{},'chr9':{},'chr10':{},'chr11':{},'chr12':{},'chr13':{},'chr14':{},'chr15':{},'chr16':{},'chr17':{},'chr18':{},'chr19':{},'chr20':{},'chr21':{},'chr22':{},'chr23':{},'chrX':{},'chrY':{}}
chr_ref = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chr23','chrX','chrY']
mu_data = {'chr1':{},'chr2':{},'chr3':{},'chr4':{},'chr5':{},'chr6':{},'chr7':{},'chr8':{},'chr9':{},'chr10':{},'chr11':{},'chr12':{},'chr13':{},'chr14':{},'chr15':{},'chr16':{},'chr17':{},'chr18':{},'chr19':{},'chr20':{},'chr21':{},'chr22':{},'chr23':{},'chrX':{},'chrY':{}}



SAM_file_name = len(file_name_info) # 0
file_name_info[SAM_file_name_In_String] = SAM_file_name

ID = 0
threshold = 2

addBAM(ID)



s = open(Output_TMP_file_name,'w')

ID = 0
temp = full_data.keys()
temp.sort()

s.write('#')
for i in file_name_info:
	s.write(i+':'+str(file_name_info[i]))
	if i != file_name_info.keys()[-1]:
		s.write(',')
s.write('\n')

s.write('#Splice'+'\n')
for chromosome in temp:
    print chromosome,
    
    data = full_data[chromosome]
    for key,value in data.items():
    	for list in value:
    		#temp_sum = sum(list[1].values())
    		#if temp_sum < threshold:
    		#	continue
    		#sum = list[1]
      		#s.write(str(key)+'\t'+str(list[0])+'\t'+str(sum)+'\n')
      		s.write(chromosome+'\t'+str(key)+'\t'+str(list[0])+'\t')
    		for i in list[1]:
    			s.write(str(i)+':'+str(list[1][i]))
    			if i != list[1].keys()[len(list[1])-1]:
    				s.write(',')
    		s.write('\t'+list[2]+'\n')
	
s.write('#Deletion'+'\n')		
for chromosome in temp:
	deletion = del_data[chromosome]
	
	for key,value in deletion.items():
		for list in value:
			#temp_sum = sum(list[1].values())
    		#if temp_sum < threshold:
    		#	continue
			s.write(chromosome+'\t'+str(key)+'\t'+str(list[0])+'\t')
			for i in list[1]:
				s.write(str(i)+':'+str(list[1][i]))
				if i != list[1].keys()[len(list[1])-1]:
					s.write(',')
			s.write('\t'+list[2]+'\n')		
      		
s.write('#Insertion'+'\n')
for chromosome in temp:
	insertion = in_data[chromosome]
	for key,value in insertion.items():
		for list in value:
			#temp_sum = sum(list[1].values())
    		#if temp_sum < threshold:
    		#	continue
			s.write(chromosome+'\t'+str(key)+'\t'+str(list[0])+'\t')
			for i in list[1]:
				s.write(str(i)+':'+str(list[1][i]))
				if i != list[1].keys()[len(list[1])-1]:
					s.write(',')
			s.write('\t'+list[2]+'\n')

s.write('#Mutation'+'\n')		
for chromosome in temp:
	Mutation = mu_data[chromosome]
	
	for key,value in Mutation.items():
		for list in value:
			s.write(chromosome+'\t'+str(key)+'\t'+str(list[0])+'\t')
			for i in list[1]:
				s.write(str(i)+':'+str(list[1][i]))
				if i != list[1].keys()[len(list[1])-1]:
					s.write(',')
			s.write('\t'+list[2]+'\n')		


s.close()

print 'END processing: ', SAM_file_name_In_String



