'''
Created on 2013. 4. 8.

@author: Akard3
'''
import sys
import os
import pickle

'''
'''

#chr1	antisense	gene	173604912	173606273	.	+	.	ID=ENSG00000232113;Name=RP11-360D2.1
#chr1	antisense	transcript	173604912	173606273	.	+	.	ID=ENST00000417563;Name=RP11-360D2.1-001;Parent=ENSG00000232113
#chr1	antisense	exon	173604912	173604984	.	+	.	ID=exon:ENST00000417563:1;Parent=ENST00000417563
#chr1	antisense	exon	173605966	173606273	.	+	.	ID=exon:ENST00000417563:2;Parent=ENST00000417563
#chr1	miRNA	gene	25732722	25732780	.	+	.	ID=ENSG00000266468;Name=AL031284.1
#chr1	miRNA	transcript	25732722	25732780	.	+	.	ID=ENST00000577655;Name=AL031284.1-201;Parent=ENSG00000266468
#chr1	miRNA	exon	25732722	25732780	.	+	.	ID=exon:ENST00000577655:1;Parent=ENST00000577655
#chr1	processed_transcript	gene	233749750	233808258	.	+	.	ID=ENSG00000135750;Name=KCNK1
#chr1	processed_transcript	transcript	233791441	233802432	.	+	.	ID=ENST00000472869;Name=KCNK1-005;Parent=ENSG00000135750
#chr1	processed_transcript	exon	233791441	233791547	.	+	.	ID=exon:ENST00000472869:1;Parent=ENST00000472869
#chr1	processed_transcript	exon	233796220	233796412	.	+	.	ID=exon:ENST00000472869:2;Parent=ENST00000472869
#chr1	processed_transcript	exon	233797602	233797682	.	+	.	ID=exon:ENST00000472869:3;Parent=ENST00000472869
#chr1	processed_transcript	exon	233802341	233802432	.	+	.	ID=exon:ENST00000472869:4;Parent=ENST00000472869
#chr1	processed_transcript	transcript	233772270	233802736	.	+	.	ID=ENST00000487728;Name=KCNK1-004;Parent=ENSG00000135750
#chr1	processed_transcript	exon	233772270	233772327	.	+	.	ID=exon:ENST00000487728:1;Parent=ENST00000487728
#chr1	processed_transcript	exon	233773034	233773180	.	+	.	ID=exon:ENST00000487728:2;Parent=ENST00000487728
#chr1	processed_transcript	exon	233802341	233802736	.	+	.	ID=exon:ENST00000487728:3;Parent=ENST00000487728
#chr1	protein_coding	mRNA	233749750	233808258	.	+	.	ID=ENST00000366621;Name=KCNK1-001;Parent=ENSG00000135750
#chr1	protein_coding	five_prime_UTR	233749750	233749917	.	+	.	ID=five_prime_UTR:ENST00000366621:1;Parent=ENST00000366621
#chr1	protein_coding	exon	233749750	233750272	.	+	.	ID=exon:ENST00000366621:1;Parent=ENST00000366621
#chr1	protein_coding	exon	233802341	233802736	.	+	.	ID=exon:ENST00000366621:2;Parent=ENST00000366621
#chr1	protein_coding	exon	233807017	233808258	.	+	.	ID=exon:ENST00000366621:3;Parent=ENST00000366621
#chr1	protein_coding	start_codon	233749918	233749920	.	+	0	ID=start_codon:ENST00000366621:1;Parent=ENST00000366621
#chr1	protein_coding	CDS	233749918	233750272	.	+	0	ID=CDS:ENST00000366621:1;Parent=ENST00000366621
#chr1	protein_coding	CDS	233802341	233802736	.	+	2	ID=CDS:ENST00000366621:2;Parent=ENST00000366621
#chr1	protein_coding	CDS	233807017	233807276	.	+	2	ID=CDS:ENST00000366621:3;Parent=ENST00000366621
#chr1	protein_coding	stop_codon	233807274	233807276	.	+	0	ID=stop_codon:ENST00000366621:1;Parent=ENST00000366621
#chr1	protein_coding	three_prime_UTR	233807277	233808258	.	+	.	ID=three_prime_UTR:ENST00000366621:1;Parent=ENST00000366621


def	CleanPeptideString(pep_str):
	#print pep_str
	pep_str = pep_str.replace("0","")
	pep_str = pep_str.replace("1","")
	pep_str = pep_str.replace("2","")
	pep_str = pep_str.replace("3","")
	pep_str = pep_str.replace("4","")
	pep_str = pep_str.replace("5","")
	pep_str = pep_str.replace("6","")
	pep_str = pep_str.replace("7","")
	pep_str = pep_str.replace("8","")
	pep_str = pep_str.replace("9","")
	pep_str = pep_str.replace("+","")
	pep_str = pep_str.replace(".","")
	pep_str = pep_str.replace("?","_")
	pep_str = pep_str.replace("_","")
	#pep_str = inserts(pep_str,".",1)
	#pep_str = inserts(pep_str,".",-1)
	#print pep_str
	return pep_str





in_folder = sys.argv[1] # folder
output_file_name = sys.argv[2]
location_count_threshold = int(sys.argv[3])

if len(sys.argv) >= 5:
	min_pep_len = int(sys.argv[4])
else:
	min_pep_len = 11


(Stub, Extension) = os.path.splitext(output_file_name)
tmp_out_file_name = Stub + "_unfiltered.tmp"


cmd_line_str = "cat " + os.path.join(in_folder,"*.txt") + " > " + tmp_out_file_name
os.system(cmd_line_str)





location_file_name = tmp_out_file_name


location_count_dic = {}
unique_location_dic = {}


location_file = open(location_file_name,'r')
for line in location_file:
    if line == "":
	continue
    if line.startswith("#"):
	continue
    if line.startswith("chr\tstart-end\tPEP\tspec_count\tlocation_count\tFDR\tSprob"):
	continue
    curr_part = line.strip()
    curr_part = curr_part.split('\t')
    peptide = curr_part[2]
    c_peptide = CleanPeptideString(peptide)

    if len(c_peptide) <= min_pep_len:
	#print c_peptide, "aa len too short"
	continue

    unique_key = c_peptide + curr_part[1]
    if unique_location_dic.has_key(unique_key):
	continue
    else:
	unique_location_dic[unique_key] = 1

    if location_count_dic.has_key(c_peptide):
	location_count_dic[c_peptide] += 1
    else:
	location_count_dic[c_peptide] = 1

location_file.close()

#chr	start-end	PEP	spec_count	location_count	FDR	Sprob



out_file = open(output_file_name,'w')
out_file.write("#chr\tstart-end\tPEP\tspec_count\tlocation_count\tFDR\tSprob\tStrand\n")

location_file = open(location_file_name,'r')

unique_location_dic = {}
for line in location_file:
    if line == "":
	continue
    if line.startswith("#"):
	continue
    if line.startswith("chr\tstart-end\tPEP\tspec_count\tlocation_count\tFDR\tSprob"):
	continue
    curr_part = line.strip()
    curr_part = curr_part.split('\t')
    peptide = curr_part[2]
    c_peptide = CleanPeptideString(peptide)
    if location_count_dic.has_key(c_peptide):
	    if location_count_threshold < location_count_dic[c_peptide]:
		continue
    else:
	continue

    unique_key = c_peptide + curr_part[1]
    if unique_location_dic.has_key(unique_key):
	continue
    else:
	unique_location_dic[unique_key] = 1

    out_file.write(curr_part[0]+"\t")
    out_file.write(curr_part[1]+"\t")
    out_file.write(curr_part[2]+"\t")
    out_file.write(curr_part[3]+"\t")
    out_file.write(str(location_count_dic[c_peptide])+"\t")
    out_file.write(curr_part[5]+"\t")
    out_file.write(curr_part[6]+"\t")
    out_file.write(curr_part[7]+"\n")

location_file.close()

out_file.close()

