#!/usr/bin/env python

import sys, os, commands, string, time, subprocess


'''
n = len(sys.argv)
for i in range(0,n):
    print sys.argv[i]
'''

input_GFF = sys.argv[1]
input_BAM = sys.argv[2]
event_type = sys.argv[3]

if sys.argv[4]=="none": #$events_file == null
    events_file = ""    
else:
    events_file = sys.argv[4]   

if sys.argv[5]=="false": #$settings_file.InputSource == "default"
    settings_file = ""    
else:
    settings_file = sys.argv[6]   #Input_data1_source=uploaded input

read_len = sys.argv[7]

if sys.argv[8]=="none": #$job_name == null
    job_name = ""    
else:
    job_name = sys.argv[8]   

htmlfile = sys.argv[9]
galaxyhome = sys.argv[10]
input_GFF_FOLDER = input_GFF.split('.')[0]+'_files'

OUTPUT_FOLDER = htmlfile.split('.')[0]+'_files'   #/media/Work/galaxy-globus-crdata/database/files/000/dataset_66_files

if not os.path.exists(OUTPUT_FOLDER):
	os.makedirs(OUTPUT_FOLDER)
else:
	pass


#command = "python ./misopy/run_events_analysis.py --compute-genes-psi "+input_GFF_FOLDER+" "+input_BAM+" --output-dir "+OUTPUT_FOLDER+" --read-len "+read_len
#print command

parameter = " --run "+input_GFF_FOLDER+" "+input_BAM+" --output-dir "+OUTPUT_FOLDER+" --read-len "+read_len

if not event_type == "default":
    parameter = parameter + " --event-type "+ event_type

if not events_file == "":
    parameter = parameter + " --events-file "+ events_file

if not settings_file == "":
    parameter = parameter + " --settings-filename "+ settings_file

if not job_name == "":
    parameter = parameter + " --job-name "+ job_name




galhtmlprefix = """<?xml version="1.0" encoding="utf-8" ?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<meta name="generator" content="Galaxy %s tool output - see http://getgalaxy.org/" />
<title></title>
<link rel="stylesheet" href="/static/style/base.css" type="text/css" />
</head>
<body>
<div class="document">
"""
galhtmlattr = """Galaxy tool %s run at %s</b><br/>"""
galhtmlpostfix = """</div></body></html>\n"""


def timenow():
    """return current time as a string
    """
    return time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))

def htmloutput(htmlout,outputfolder):
	rstyle="""<style type="text/css">
        tr.d0 td {background-color: oldlace; color: black;}
        tr.d1 td {background-color: aliceblue; color: black;}
        </style>"""    
        res = [rstyle,]
        res.append(galhtmlprefix % os.path.basename(sys.argv[0]))   
        res.append(galhtmlattr % ('Compute Psi values for genes',timenow()))
        flist = [x for x in os.listdir(outputfolder) if not x.startswith('.')] 

	for root, dirs, files in os.walk(outputfolder): 
		for name in files: 
			print os.path.join(root, name)  

        if len(flist) > 0:
            res.append('<b>The following output files were created (click the filename to view/download a copy):</b><hr/>')
            res.append('<table>\n')

	    for root, dirs, files in os.walk(outputfolder): 
		for name in files: 
		    if cmp(root,outputfolder) == 0:  #if this is a file at outputfolder
                        res.append('<tr><td><a href="%s">%s</a></td></tr>\n' % (name,name))
                    else:  #if this is a dir at outputfolder
			subdir = root.replace(outputfolder,'')[1:]
			res.append('<tr><td><a href="%s">%s</a></td></tr>\n' % (os.path.join(subdir, name),name))

            res.append('</table><p/>\n') 

        res.append(galhtmlpostfix) 
        outf = open(htmlout,'w')
        outf.write(''.join(res))   
        outf.write('\n')
        outf.close()



##uncomment when running on VM

os.chdir("%s/tools/miso/" % galaxyhome )     #uncomment when running on cluster
#os.chdir("/media/Work/galaxy-proteonics/tools/miso/")    #for local running

command="python2.7 ./miso.py "+parameter

#command = "python2.7 ./run_events_analysis.py  --compute-genes-psi /media/Work/galaxy-proteonics/database/files/000/dataset_150_files /media/Work/galaxy-proteonics/database/files/000/dataset_144_files/miso_BamBai.bam --output-dir /media/Work/galaxy-proteonics/database/files/000/dataset_154_files --read-len 36 --job-name job"

print command

subprocess.call(command,shell=True)

htmloutput(htmlfile,OUTPUT_FOLDER)


