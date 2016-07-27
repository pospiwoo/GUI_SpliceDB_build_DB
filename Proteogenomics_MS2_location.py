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



# Function definitions from here



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

'''
def inReduce(CIGAR,CIGAR_begin):
	in_location = []
	pattern = re.compile('S|N|M|D|I|[0-9]+')
	CIGAR_info = [(m.group()) for m in pattern.finditer(CIGAR)]
	CIGAR_interval = [int(n) for i,n in enumerate(CIGAR_info) if i%2 == 0]
	index_D = [i for i,n in enumerate(CIGAR_info) if n == 'D']
	for i in index_D:
		in_interval = CIGAR_info[i-1]
		in_start = CIGAR_begin + sum(CIGAR_interval[0:(i-1)/2])
		in_end = in_start + in_interval
		in_location.append((in_start,in_end))
	return in_location

def check_mutation(M_string,DNA_string,string_begin,line):
    mu_info = []
    for i,j in enumerate(M_string):
        if j != DNA_string[i]:
            mu_info.append((string_begin + i,j))
            #mu_location =  string_begin + i - len(mu_string)
            #mu_string += j
            
    if len(mu_info) > 5:
    	print line
    	print M_string
    	print DNA_string
    	print string_begin
    	return []
    return mu_info #[mu_location-1, mu_string]
'''

def addGFF(SourceFile,ID):
	# needs to figure out how to distinguish insertion, deletion and splice information.
	return

def addSAM(SourceFile,ID):
	id = 1
	for line in SourceFile:
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
	    
	    if Species == 'Human':
	    	#chr_name = formatChrNameForHuman(RNAME)
	    	chr_name = RNAME
	    	
	    elif Species == 'CElegans':
    		chr_name = formatChrNameForCElegans(RNAME)
	    
	    if chr_name not in chr_ref:
	    	continue
	    
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
	    
	    	
	    if 'I' in CIGAR:
	    	#in_string = 'Bla Bla,' # calculate insertion string from location information and input strng
	    	sign = 'I'
	    	[in_location,seq_location] = decomp(CIGAR,CIGAR_begin,sign,strand)
	    	in_temp = in_data.get(chr_name)
	    	#print in_temp,in_location,CIGAR,chr_name
	    	for i,coord_info in enumerate(in_location):
	    		in_string = SEQ[seq_location[i][0]:seq_location[i][1]]
	    		#print in_temp
	    		if in_temp.has_key(coord_info[0]):
	    			check = 0
	    			for i,check_end in enumerate(in_temp[coord_info[0]]):
	    				#print in_temp[coord_infom[0]], check_end, coord_infom, CIGAR
	    				if check_end[0] == in_string and strand == check_end[2]:
	    					if check_end[1].has_key(SAM_file_name):
	    						in_temp[coord_info[0]][i][1][SAM_file_name] = in_temp[coord_info[0]][i][1].get(SAM_file_name)+1
    						else:
    							in_temp[coord_info[0]][i][1][SAM_file_name] = 1
	    					#in_temp[coord_infom[0]][i][1] += 1
	    					check = 1
	    			if check == 0:
	    				in_temp[coord_info[0]].append([in_string,{SAM_file_name:1},strand])
	    		else:
	    			in_temp[coord_info[0]] = [[in_string,{SAM_file_name:1},strand]]
	    						
	    	in_data[chr_name] = in_temp
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
    	
    	#mutation detection and recording
    	'''
		if 'M' in CIGAR:
			sign = 'M'
			[M_location, seq_location] = decomp(CIGAR,CIGAR_begin,sign,strand)
			mu_temp = mu_data.get(chr_name)
	    	for i,coord_info in enumerate(M_location):
				M_string = SEQ[seq_location[i][0]:seq_location[i][1]]
				DNA_string = DNA_data[M_location[i][0]:M_location[i][1]]
				mu_info = check_mutation(M_string,DNA_string,M_location[i][0],line)
				if mu_info == []:
					continue
				for j in mu_info:
					if mu_temp.has_key(j[0]):
						check = 0
						for k,check_end in enumerate(mu_temp[j[0]]):
							if check_end[0] == j[1] and strand == check_end[2]:
								if check_end[1].has_key(SAM_file_name):
									mu_temp[j[0]][k][1][SAM_file_name] += 1
								else:
									mu_temp[j[0]][k][1][SAM_file_name] = 1
								#in_temp[coord_infom[0]][i][1] += 1
								check = 1
						if check == 0:
							mu_temp[j[0]].append([j[1],{SAM_file_name:1},strand])
					else:
						mu_temp[j[0]] = [[j[1],{SAM_file_name:1},strand]]
				
	    	mu_data[chr_name] = mu_temp
	    '''


def addBAM(ID):
	id = 1
	for line in sys.stdin: #SourceFile:
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
    	


def addTMP(SourceFile,ID,file_name_info):
	ignore_list = []
	table_file_name = {}
	input_file_name_info = {}
	file_temp = SourceFile.readline()
	file_temp = file_temp[1:]
	if file_temp != '':
		file_temp = file_temp.strip()
		file_temp = file_temp.split(',')
		for i in file_temp:
			i = i.split(':')
			#print i
			input_file_name_info[i[0]] = int(i[1])
	for new_file_name in input_file_name_info:
		if file_name_info.has_key(new_file_name):
			print new_file_name + ' already in file'
			ignore_list.append(input_file_name_info[new_file_name])
			return
		else:
			file_name_number = len(file_name_info)
			file_name_info[new_file_name] = file_name_number
			table_file_name[input_file_name_info[new_file_name]] = file_name_number
	## table_file_name >> convert from input file name info to new added file name info
	## ignore_list >> contain information of overlapped file name info
	id = 0
	for line in SourceFile:
		if check_sample[0] == 1:
			id += 1
			if id > check_sample[1]:
				break
		line = line.strip()
		if re.match(line,'#Splice'):
			case = 0
			continue
		elif re.match(line,'#Deletion'):
			case = 1
			continue
		elif re.match(line,'#Insertion'):
			case = 2
			continue
		elif re.match(line,'#Mutation'):
			case = 3
			continue
		line = line.split('\t')
		chromosome = line[0]
		first = int(line[1])
		if case != 2 and case != 3:
			second = int(line[2])
		else:
			second = line[2]
		files = line[3]
		strand = line[4]
		files = files.split(',')
		file_info = {}
		for i in files:
			i = i.split(':')
			file_info[int(i[0])] = int(i[1])
		if case == 0:
			if full_data.has_key(chromosome):
				tmp_data = full_data[chromosome]
			else:
				continue
		elif case == 1:
			if del_data.has_key(chromosome):
				tmp_data = del_data[chromosome]
			else:
				continue
		elif case == 2:
			if in_data.has_key(chromosome):
				tmp_data = in_data[chromosome]
			else:
				continue
		elif case == 3:
			if mu_data.has_key(chromosome):
				tmp_data = mu_data[chromosome]
			else:
				continue
		new_file_info = {}
		for i in file_info:
			if i not in ignore_list:
				new_file_info[table_file_name[i]] = file_info[i]
		if tmp_data.has_key(first):
			check = 0
			for j,check_end in enumerate(tmp_data[first]):
				if check_end[0] == second and check_end[2] == strand:
					for i in new_file_info:
						#print ':j:',j,':new_file_info:',new_file_info,':i:',i,':tmp_data[first][j]',tmp_data[first][j]
						tmp_data[first][j][1][i] = new_file_info[i]
					check = 1
			if check == 0:
				tmp_data[first].append([second,new_file_info,strand])
		else:
			tmp_data[first] = [[second,new_file_info,strand]]	

'''
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
'''

def addVCF(SourceFile,ID):
	for line in SourceFile:
		if line.startswith('#'):
			continue
		line = line.split('\t')
		chromosome = 'chr'+str(line[0].strip())
		if chromosome not in chr_ref:
			print 'Chromosome: ', chromosome,' is not readable'
			continue
		strand = '.'
		coordi = int(line[1])
		REF = line[3]
		ALT_g = line[4].split(',')
		for ALT in ALT_g:
			if len(REF) == 1 and len(ALT) == 1:
				data = mu_data.get(chromosome)
				if data.has_key(coordi):
					check = 0
					for i, check_end in enumerate(data[coordi]):
						if check_end[0] == ALT:
							data[coordi][i][1][SAM_file_name] = 1
							check = 1
					if check == 0:
						data[coordi].append([ALT,{SAM_file_name:1},strand])
				else:
					data[coordi] = [[ALT,{SAM_file_name:1},strand]]
					
			elif len(REF) >  len(ALT):
				data = del_data.get(chromosome)
				coordi = coordi + 1
				coordi_end = coordi+len(REF)-len(ALT)
				if data.has_key(coordi):
					check = 0
					for i, check_end in enumerate(data[coordi]):
						if check_end[0] == coordi_end:
							data[coordi][i][1][SAM_file_name] = 1
							check = 1
					if check == 0:
						data[coordi].append([coordi_end,{SAM_file_name:1},strand])
				else:
					data[coordi] =  [[coordi_end,{SAM_file_name:1},strand]]
					
			elif len(REF) < len(ALT):
				data = in_data.get(chromosome)
				coordi = coordi + len(REF)
				if data.has_key(coordi):
					for i, check_end in enumerate(data[coordi]):
						if check_end[0] == ALT[len(REF):]:
							data[coordi][i][1][SAM_file_name] = 1
							check = 1
					if check == 0:
						data[coordi].append([ALT[len(REF):],{SAM_file_name:1},strand])
				else:
					data[coordi] = [[ALT[len(REF):],{SAM_file_name:1},strand]]
					
	return























print 'begins: ',time.ctime()

#BAM_SAM_TMP_file_directory = sys.argv[1] #(should be a folder) #-a
#Output_TMP_file_name = sys.argv[2] #-p
#path_to_sam_tools = sys.argv[3] #-t
#path_to_dna_fasta_files = sys.argv[4] #-f
#path_to_output_sixframe_files = sys.argv[5] #-s
#Known_protein_aa_FASTA_file_name = sys.argv[6] #-k
#path_to_MS2_spectra = sys.argv[7] #-m
#path_to_MS2_results = sys.argv[8] #-r
#path_to_location_results = sys.argv[9] #-l
#final_location_result_file_name = sys.argv[10] #-o

if len(sys.argv)%2 == 0:
	print "number of input parameters are not even"
	exit()

BAM_SAM_TMP_file_directory = ""
Output_TMP_file_name = ""
path_to_sam_tools = ""
path_to_dna_fasta_files = ""
path_to_output_sixframe_files = ""
Known_protein_aa_FASTA_file_name = ""
path_to_MS2_spectra = ""
path_to_MS2_results = ""
path_to_location_results = ""
final_location_result_file_name = ""
for i in range(0,len(sys.argv)/2):
	if sys.argv[i*2+1] == "-a":
		BAM_SAM_TMP_file_directory = sys.argv[i*2+2] #(should be a folder) #-a
		#print "-a", BAM_SAM_TMP_file_directory
	elif sys.argv[i*2+1] == "-p":
		Output_TMP_file_name = sys.argv[i*2+2] #-p
	elif sys.argv[i*2+1] == "-t":
		path_to_sam_tools = sys.argv[i*2+2] #-t
	elif sys.argv[i*2+1] == "-f":
		path_to_dna_fasta_files = sys.argv[i*2+2] #-f
	elif sys.argv[i*2+1] == "-s":
		path_to_output_sixframe_files = sys.argv[i*2+2] #-s
	elif sys.argv[i*2+1] == "-k":
		Known_protein_aa_FASTA_file_name = sys.argv[i*2+2] #-k
	elif sys.argv[i*2+1] == "-m":
		path_to_MS2_spectra = sys.argv[i*2+2] #-m
	elif sys.argv[i*2+1] == "-r":
		path_to_MS2_results = sys.argv[i*2+2] #-r
	elif sys.argv[i*2+1] == "-l":
		path_to_location_results = sys.argv[i*2+2] #-l
		#print "-l", path_to_location_results
	elif sys.argv[i*2+1] == "-o":
		final_location_result_file_name = sys.argv[i*2+2] #-o
	else:
		print "unrecognized parameter", sys.argv[i*2]


#if BAM_SAM_TMP_file_directory == "":
#	print "bam/sam/spl directory not specified"
#	exit()
if Output_TMP_file_name == "":
	print "Output spl file name or folder not specified"
	exit()
#if path_to_sam_tools == "":
#	print "path to samtools not specified"
#	exit()
#if path_to_dna_fasta_files == "":
#	print "path to dna fasta files not specified"
#	exit()
if path_to_output_sixframe_files == "":
	print "path to sixframe not specified"
	exit()
if Known_protein_aa_FASTA_file_name == "":
	print "Known protein aa FASTA file name not specified"
	exit()
if path_to_MS2_spectra == "":
	print "path to MS2 spectra not specified"
	exit()
if path_to_MS2_results == "":
	print "path to output MS2 results not specified"
	exit()
if path_to_location_results == "":
	print "path to output location results not specified"
	exit()
if final_location_result_file_name == "":
	print "final location result file name not specified"
	exit()


# parameters and setting

check_sample = [0,1000] #0: check off, 1: check on

threshold = 2
Species = 'Human' #Select from, 'CElegans', 'Human'
java_mem_in_MB = 3000 
max_genomic_location_count = 3

chr_ref = []
dna_TRIE_file_name_per_chr = {}
Sixframe_file_name_per_chr = {}



## Followings are the dna .trie files to be used (must specify the exact path and file name)
#dna_trie_files = {}
#print "  -Searching for reference dna .trie files"
#Files = os.listdir(path_to_dna_fasta_files)
#
#
#for F in Files:
#	if F.find(".") < 0 or os.path.isdir(F):
#		continue
#	file_path_and_name = os.path.join(path_to_dna_fasta_files,F)
#	(Stub, Extension) = os.path.splitext(file_path_and_name)
#	(Stub0, Extension0) = os.path.splitext(F)
#
#	if F.endswith(".fa") or F.endswith(".fasta"):
#		if F.startswith("chr") or F.startswith("Chr") or F.startswith("CHR"):
#			chr_ref_tmp = F.split("_")[0]
#			chr_ref.append(chr_ref_tmp)
#			dna_TRIE_file_name_per_chr[chr_ref_tmp] = Stub + ".trie"
#
#			tmp_Sixframe_file_name_per_chr = "Sixframe_" + Stub0 + ".fasta"
#			Sixframe_file_name_per_chr[chr_ref_tmp] = os.path.join(path_to_output_sixframe_files,tmp_Sixframe_file_name_per_chr)
#
#			if os.path.exists(Stub + ".trie") == False:
#				cmd_line_str = "python PrepDB.py FASTA " + file_path_and_name
#				os.system(cmd_line_str)
#
#			if os.path.exists(Stub + "_0.fasta") == False:
#				print " Building Sixframe DB"
#				cmd_line_str = "java -jar SixFrameBuilder.jar -r " + file_path_and_name + " -w " + Sixframe_file_name_per_chr[chr_ref_tmp]
#				os.system(cmd_line_str)
#
## end of parameters












# Using chr representation, we initiate the data structure 
full_data = {}
in_data = {}
del_data = {}
mu_data = {}
for i in range(0,len(chr_ref)):
	full_data[chr_ref[i]] = {}
	in_data[chr_ref[i]] = {}
	del_data[chr_ref[i]] = {}
	mu_data[chr_ref[i]] = {}




###Retreive folder
#
#FileNames = []
#Files = os.listdir(BAM_SAM_TMP_file_directory)
#file_name_info = {}
#
## find .spl file first
#Files = os.listdir(BAM_SAM_TMP_file_directory)
#for F in Files:
#	if F.find(".") < 0 or os.path.isdir(F):
#		continue
#	if F.endswith(".bam") or F.endswith(".BAM"):
#		continue
#	
#	curr_file_name = F
#	if file_name_info.has_key(curr_file_name):
#		print 'Current file already in memory :', curr_file_name
#		continue
#	SAM_file_name = len(file_name_info)
#
#	if F.endswith(".spl") :
#		file_path_and_name = os.path.join(BAM_SAM_TMP_file_directory,F)
#		print " -Processing .spl file :",F
#		SourceFile = open(file_path_and_name,'r')
#		ID = 0
#		addTMP(SourceFile,ID,file_name_info)
#		continue
#	#else:
#	#	print 'Wrong input file format : ' + F
#	#	continue
#
#
## find BAM file next
#for F in Files:
#	if F.find(".") < 0 or os.path.isdir(F):
#		continue
#	elif F.endswith(".bam") or F.endswith(".BAM"):
#		file_path_and_name = os.path.join(BAM_SAM_TMP_file_directory,F)
#		(Stub, Extension) = os.path.splitext(F)
#		Single_BAM_Output_TMP_file_name = Stub + ".spl"
#		
#		#print Single_BAM_Output_TMP_file_name
#		if os.path.exists(BAM_SAM_TMP_file_directory + "/" + Single_BAM_Output_TMP_file_name):
#			print "BAM file :", Single_BAM_Output_TMP_file_name, "already exists"
#			continue
#		elif file_name_info.has_key(F):
#			print 'Current file already in memory :', F
#			continue
#		elif file_name_info.has_key(F.replace(".bam",".sam")):
#			print 'Current file already in memory :', F
#			continue
#		elif file_name_info.has_key(F.replace(".sam",".bam")):
#			print 'Current file already in memory :', F
#			continue
#		elif file_name_info.has_key(F.replace(".spl",".bam")):
#			print 'Current file already in memory :', F
#			continue
#		elif file_name_info.has_key(F.replace(".spl",".sam")):
#			print 'Current file already in memory :', F
#			continue
#
#		#print file_name_info
#		print " -Processing .bam file :",F
#		cmd_line_str = path_to_sam_tools + "/samtools view " + file_path_and_name +" | "+ "python SAMtoGFF_indel_vcf_single_BAM_file_suwoo_ver_20130724.py " + F + " " + Single_BAM_Output_TMP_file_name
#
#		#cmd_line_str = os.path.join(path_to_sam_tools,samtools)
#		#cmd_line_str = cmd_line_str + " view " + BAM_file_name + " | python SAMtoGFF_indel_vcf_single_BAM_file_suwoo_ver_20130724.py " + Output_TMP_file_name
#		
#		os.system(cmd_line_str)
#
#		#ID = 0
#		#addBAM(ID)
#		continue
#
#
#
#
## find other files
#Files = os.listdir(BAM_SAM_TMP_file_directory)
#for F in Files:
#	if F.find(".") < 0 or os.path.isdir(F):
#		continue
#	if F.endswith(".bam") or F.endswith(".BAM"):
#		continue
#	
#	curr_file_name = F
#	if file_name_info.has_key(curr_file_name):
#		print 'Current file already in memory :', curr_file_name
#		continue
#	SAM_file_name = len(file_name_info)
#
#
#	if F.endswith(".spl") :
#		file_path_and_name = os.path.join(BAM_SAM_TMP_file_directory,F)
#		print " -Processing .spl file :",F
#		SourceFile = open(file_path_and_name,'r')
#		ID = 0
#		addTMP(SourceFile,ID,file_name_info)
#		continue
#	elif F.endswith(".sam") or F.endswith(".SAM"):
#		file_name_info[curr_file_name] = SAM_file_name
#		file_path_and_name = os.path.join(BAM_SAM_TMP_file_directory,F)
#		print " -Processing .sam file :",F
#		SourceFile = open(file_path_and_name,'r')
#		ID = 0
#		addSAM(SourceFile,ID)
#		continue
#	elif F.endswith(".vcf") :
#		file_name_info[curr_file_name] = SAM_file_name
#		file_path_and_name = os.path.join(BAM_SAM_TMP_file_directory,F)
#		print " -Processing .vcf file :",F
#		SourceFile = open(file_path_and_name,'r')
#		ID = 0
#		addVCF(SourceFile,ID)
#		continue
#	#else:
#	#	print 'Wrong input file format : ' + F
#	#	continue
#
#
#
#
#
#
#
#
#temp = full_data.keys()
#temp.sort()
#for chromosome in temp:	
#	(Stub, Extension) = os.path.splitext(Output_TMP_file_name)
#	Output_TMP_file_name_per_chr = Stub+"_"+chromosome+".spl"
#	Output_ms2db_file_name_per_chr = Stub+"_"+chromosome+".ms2db"
#	Output_fasta_file_name_per_chr = Stub+"_"+chromosome+".fa"
#	#Output_TMP_file_name_per_chr = Output_TMP_file_name.replace(".spl","_"+chromosome+".spl")
#
#	print " -Writting output .spl file :", Output_TMP_file_name_per_chr
#	s = open(Output_TMP_file_name_per_chr,'w')
#
#	s.write('#')
#	for i in file_name_info:
#		s.write(i+':'+str(file_name_info[i]))
#		if i != file_name_info.keys()[-1]:
#			s.write(',')
#	s.write('\n')
#
#
#	s.write('#Splice'+'\n')
#	data = full_data[chromosome]
#	for key,value in data.items():
#		for list in value:
#			s.write(chromosome+'\t'+str(key)+'\t'+str(list[0])+'\t')
#			for i in list[1]:
#				s.write(str(i)+':'+str(list[1][i]))
#				if i != list[1].keys()[len(list[1])-1]:
#					s.write(',')
#			s.write('\t'+list[2]+'\n')
#
#
#
#
#
#
#	
#	s.write('#Deletion'+'\n')
#	deletion = del_data[chromosome]	
#	for key,value in deletion.items():
#		for list in value:
#			#temp_sum = sum(list[1].values())
#    		#if temp_sum < threshold:
#    		#	continue
#			s.write(chromosome+'\t'+str(key)+'\t'+str(list[0])+'\t')
#			for i in list[1]:
#				s.write(str(i)+':'+str(list[1][i]))
#				if i != list[1].keys()[len(list[1])-1]:
#					s.write(',')
#			s.write('\t'+list[2]+'\n')		
#      		
#	s.write('#Insertion'+'\n')
#	insertion = in_data[chromosome]
#	for key,value in insertion.items():
#		for list in value:
#			#temp_sum = sum(list[1].values())
#    		#if temp_sum < threshold:
#    		#	continue
#			s.write(chromosome+'\t'+str(key)+'\t'+str(list[0])+'\t')
#			for i in list[1]:
#				s.write(str(i)+':'+str(list[1][i]))
#				if i != list[1].keys()[len(list[1])-1]:
#					s.write(',')
#			s.write('\t'+list[2]+'\n')
#
#	s.write('#Mutation'+'\n')
#	Mutation = mu_data[chromosome]
#	
#	for key,value in Mutation.items():
#		for list in value:
#			s.write(chromosome+'\t'+str(key)+'\t'+str(list[0])+'\t')
#			for i in list[1]:
#				s.write(str(i)+':'+str(list[1][i]))
#				if i != list[1].keys()[len(list[1])-1]:
#					s.write(',')
#			s.write('\t'+list[2]+'\n')
#
#	s.close()
#
#
#
#
#
#
#	print "  -Constructing splice graph from .spl file"
#
#	cmd_line_str = "java -Xmx" + str(java_mem_in_MB) + "M -jar ConstructSpliceGraphFromTmp.jar"
#	cmd_line_str = cmd_line_str + " -p " + Output_TMP_file_name_per_chr
#	cmd_line_str = cmd_line_str + " -w " + Output_ms2db_file_name_per_chr
#	cmd_line_str = cmd_line_str + " -s " + dna_TRIE_file_name_per_chr[chromosome]
#	
#	#print "    ",cmd_line_str
#	os.system(cmd_line_str)
#
#
#
#
#
#
#	print "  -Converting splice graph(.spl) to .fa file"
#	cmd_line_str = "python ACGT03102013_split.py " + Output_ms2db_file_name_per_chr + " " + Output_fasta_file_name_per_chr + " 30"
#
#	#print "    ",cmd_line_str
#	os.system(cmd_line_str)













print "  -Check knwon protein .fa file"
(Stub, Extension) = os.path.splitext(Known_protein_aa_FASTA_file_name)
known_protein_file_name = Stub + ".trie"
if os.path.exists(known_protein_file_name) == False:
	cmd_line_str = "python PrepDB.py FASTA " + Known_protein_aa_FASTA_file_name
	print "    - test : ",cmd_line_str
	os.system(cmd_line_str)







print "  -Check and run MSGF+ on splice graph .fa files"
(Stub, Extension) = os.path.splitext(Output_TMP_file_name)
splice_dir = os.path.dirname(Output_TMP_file_name)

Files = os.listdir(splice_dir)

for F in Files:
	if F.find(".") < 0 or os.path.isdir(F):
		continue
	file_path_and_name = os.path.join(splice_dir,F)
	(Stub, Extension) = os.path.splitext(F)

	if F.endswith(".fa") or F.endswith(".fasta"):
		if F.find(".revCat") > 0 or F.endswith(".revCat") > 0:
			continue
		sFiles = os.listdir(path_to_MS2_spectra)
		for sF in sFiles:
			if sF.find(".") < 0 or os.path.isdir(sF):
				continue

			(sStub, sExtension) = os.path.splitext(sF)
			out_MS2_file_name = sStub + "_" + Stub + ".mzid"
			cmd_line_str = "java -Xmx"+str(java_mem_in_MB)+"M -jar MSGFPlus.jar"
			cmd_line_str = cmd_line_str + " -s " + os.path.join(path_to_MS2_spectra,sF)
			cmd_line_str = cmd_line_str + " -o " + os.path.join(path_to_MS2_results,out_MS2_file_name)
			cmd_line_str = cmd_line_str + " -d " + file_path_and_name
			cmd_line_str = cmd_line_str + " -tda 1 "
			cmd_line_str = cmd_line_str + " -mod Mods.txt"

			#print "    ",cmd_line_str
			print "    - test : ",cmd_line_str
			os.system(cmd_line_str)





print "  -Check and run MSGF+ on sixframe .fa files"
Files = os.listdir(path_to_output_sixframe_files)

for F in Files:
	if F.find(".") < 0 or os.path.isdir(F):
		continue
	file_path_and_name = os.path.join(path_to_output_sixframe_files,F)
	(Stub, Extension) = os.path.splitext(F)

	if F.endswith(".fa") or F.endswith(".fasta"):
		if F.find(".revCat") > 0 or F.endswith(".revCat") > 0:
			continue
		sFiles = os.listdir(path_to_MS2_spectra)
		for sF in sFiles:
			if sF.find(".") < 0 or os.path.isdir(sF):
				continue

			(sStub, sExtension) = os.path.splitext(sF)
			out_MS2_file_name = sStub + "_" + Stub + ".mzid"
			cmd_line_str = "java -Xmx"+str(java_mem_in_MB)+"M -jar MSGFPlus.jar"
			cmd_line_str = cmd_line_str + " -s " + os.path.join(path_to_MS2_spectra,sF)
			cmd_line_str = cmd_line_str + " -o " + os.path.join(path_to_MS2_results,out_MS2_file_name)
			cmd_line_str = cmd_line_str + " -d " + file_path_and_name
			cmd_line_str = cmd_line_str + " -tda 1"
			cmd_line_str = cmd_line_str + " -mod Mods.txt"

			#print "    ",cmd_line_str
			print "    - test : ",cmd_line_str
			os.system(cmd_line_str)









print "  -Merge result files"
out_Merged_MS2_file_name = "Merged_MS2_results.tsv"
cmd_line_str = "python ./MzidToTsv.py " + path_to_MS2_results + " "
cmd_line_str = cmd_line_str + os.path.join(path_to_MS2_results,out_Merged_MS2_file_name)

print "    - test : ",cmd_line_str
os.system(cmd_line_str)





print "  -Find novel peptide psm"
(Stub, Extension) = os.path.splitext(out_Merged_MS2_file_name)
novel_psm_file_name = Stub + "_novel.txt"
known_psm_file_name = Stub + "_known.txt"

cmd_line_str = "java -Xmx"+str(java_mem_in_MB)+"M -jar DetermineNovelty_minus_two_peptide.jar"
cmd_line_str = cmd_line_str + " -r " + os.path.join(path_to_MS2_results,out_Merged_MS2_file_name)
cmd_line_str = cmd_line_str + " -w " + os.path.join(path_to_MS2_results,known_psm_file_name)
cmd_line_str = cmd_line_str + " -n " + os.path.join(path_to_MS2_results,novel_psm_file_name)
cmd_line_str = cmd_line_str + " -t " + known_protein_file_name

print "    - test : ",cmd_line_str
os.system(cmd_line_str)





print "FDR calculation and cut-off"
(Stub, Extension) = os.path.splitext(novel_psm_file_name)
novel_FDR_psm_file_name = Stub + "_FDR.txt"

cmd_line_str = "java -Xmx"+str(java_mem_in_MB)+"M -cp MSGFDB.jar fdr.ComputeFDR"
cmd_line_str = cmd_line_str + " -f " + os.path.join(path_to_MS2_results,out_Merged_MS2_file_name)
cmd_line_str = cmd_line_str + " 9 XXX -i 0 -n 2 -p 8 -s 12 0 -pepfdr 0.01 "
cmd_line_str = cmd_line_str + "-o " + os.path.join(path_to_MS2_results,novel_FDR_psm_file_name)

print "    - test : ",cmd_line_str
os.system(cmd_line_str)




print "Convert to .p file"
(Stub, Extension) = os.path.splitext(novel_FDR_psm_file_name)
p_file_novel_FDR_psm_file_name = Stub + ".p"

cmd_line_str = "python ConvertTSVToPickleFile.py " + os.path.join(path_to_MS2_results,novel_FDR_psm_file_name) + " "
cmd_line_str = cmd_line_str + os.path.join(path_to_MS2_results,p_file_novel_FDR_psm_file_name)

print "    - test : ",cmd_line_str
os.system(cmd_line_str)





print "Search location"



(Stub, Extension) = os.path.splitext(Output_TMP_file_name)
splice_dir = os.path.dirname(Output_TMP_file_name)

Files = os.listdir(splice_dir)

for F in Files:
	if F.find(".") < 0 or os.path.isdir(F):
		continue
	file_path_and_name = os.path.join(splice_dir,F)
	(Stub, Extension) = os.path.splitext(F)
	(pStub, pExtension) = os.path.splitext(p_file_novel_FDR_psm_file_name)
	location_file_name = "Location_" + Stub + "_" + pStub + ".txt"

	if F.endswith(".fa") or F.endswith(".fasta"):
		if F.find(".revCat") > 0 or F.endswith(".revCat") > 0:
			continue
		cmd_line_str = "python Location05242013.py "
		cmd_line_str = cmd_line_str + file_path_and_name + " "
		cmd_line_str = cmd_line_str + os.path.join(path_to_location_results,location_file_name) + " "
		cmd_line_str = cmd_line_str + os.path.join(path_to_MS2_results,p_file_novel_FDR_psm_file_name)
		os.system(cmd_line_str)

		print "    - test : ",cmd_line_str
		os.system(cmd_line_str)





Files = os.listdir(path_to_output_sixframe_files)

for F in Files:
	if F.find(".") < 0 or os.path.isdir(F):
		continue
	file_path_and_name = os.path.join(path_to_output_sixframe_files,F)
	(Stub, Extension) = os.path.splitext(F)
	(pStub, pExtension) = os.path.splitext(p_file_novel_FDR_psm_file_name)
	location_file_name = "Location_" + Stub + "_" + pStub + ".txt"

	if F.endswith(".fa") or F.endswith(".fasta"):
		if F.find(".revCat") > 0 or F.endswith(".revCat") > 0:
			continue
		cmd_line_str = "python Location05242013.py "
		cmd_line_str = cmd_line_str + file_path_and_name + " "
		cmd_line_str = cmd_line_str + os.path.join(path_to_location_results,location_file_name) + " "
		cmd_line_str = cmd_line_str + os.path.join(path_to_MS2_results,p_file_novel_FDR_psm_file_name)
		os.system(cmd_line_str)

		print "    - test : ",cmd_line_str
		os.system(cmd_line_str)
#








print "Recalculate genomic location count and filter"

cmd_line_str = "python RecalculateLocationCountAndMerge.py " + path_to_location_results + " "
cmd_line_str = cmd_line_str + final_location_result_file_name + " "
cmd_line_str = cmd_line_str + str(max_genomic_location_count)
os.system(cmd_line_str)

print "    - test : ",cmd_line_str
os.system(cmd_line_str)














print ''
print 'END: ',time.ctime()



















