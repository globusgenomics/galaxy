import os, sys
import fnmatch
import csv
	
def get_headers(inputfile):	
	columnList=[]
	#line=inputfile.readlines()[0]
	filename=inputfile.get_file_name()
	try:
		f = open(filename)
		line=f.readline()	
		while(line[0]=='#' or (not line.strip())):	#remove header (starting with hash sign and empty lines to get to headerline
			line=f.readline()
		line = line.strip()		
		i=1;
		for col in line.split("\t"):
			label=str(i)+': '+str(col)
			columnList.append([label,label,False])		
			i+=1		
		
	except IOError as e:	
		pass	
	
	return columnList
	
	




