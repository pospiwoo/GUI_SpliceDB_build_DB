import string
import os
import sys
import struct
import re

python_script_path_name = sys.argv[1]
db_file_names = sys.argv[2]
p_file_name = sys.argv[3]
output_sh_file_name = sys.argv[4]


#out_file_path = Stub + "_formatted.fa"
SourceFile = open(db_file_names, "r")
OutFile  = open(output_sh_file_name, "w")

FileNames = []
FileNames_prefix = []
for line in SourceFile:
	data = line.strip()
	if data == "":
		break
	if line == "":
		break
	FileNames.append(data)
	
	data = data.split("/")[-1]
	(Stub, Extension) = os.path.splitext(data)
	(Stub1, Extension) = os.path.splitext(output_sh_file_name)
	out_file_name = "Location_" + Stub + ".txt"
	FileNames_prefix.append(out_file_name)

#	if data.find(".fasta") > 0:
#		FileNames_prefix.append(data.replace(".fasta","_Location.txt"))
#		FileNames_prefix.append(data.replace("/data/repository/CPTAC/Database/SpliceGraph/PNNL_OV_Match/","/data/home/s3cha/Postprocess/debug/PNNL_Incoming_20130505/"))
#	elif data.find(".fa") > 0:
#		FileNames_prefix.append(data.replace(".fa","_Location.txt"))
#		FileNames_prefix.append(data.replace("/data/repository/CPTAC/Database/SpliceGraph/PNNL_OV_Match/","/data/home/s3cha/Postprocess/debug/PNNL_Incoming_20130505/"))
#	else:
#		print data



for i in range(0,len(FileNames)):
	output_string = "qsub -l h_vmem=9G /data/home/s3cha/Postprocess/run_cmd.sh python "
	output_string = output_string + python_script_path_name + " "
	output_string = output_string + FileNames[i] + " "
	output_string = output_string + FileNames_prefix[i] + " "
	output_string = output_string + p_file_name + "\n"
	OutFile.write(output_string)


OutFile.close()


