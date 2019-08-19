#!/usr/bin/env python

import sys, os, commands, string, time, subprocess

galaxyhome=os.environ.get('GALAXY_HOME')


n = len(sys.argv)
for i in range(0,n):
    #print sys.argv[i]
    if sys.argv[i]=="flag":   #flag is used to judge the number of countfile (may have multiple countfile)
	m=i


labelname = sys.argv[1]

#compose multiple countfile to one string
countfile=""
for j in range(2,m):
    print sys.argv[j]
    countfile += sys.argv[j]+","

#remove last comma
countfile=countfile[0:-1]

event_type = sys.argv[m+1]

read_len = sys.argv[m+2]
overhang_len = sys.argv[m+3]


if sys.argv[m+4]=="none": #$events_file == null
    events_file = ""    
else:
    events_file = sys.argv[m+4]   

if sys.argv[m+5]=="false": #$settings_file.InputSource == "default"
    settings_file = ""    
else:
    settings_file = sys.argv[m+6]   #Input_data1_source=uploaded input


if sys.argv[m+7]=="none": #$job_name == null
    job_name = ""    
else:
    job_name = sys.argv[m+7]   

htmlfile = sys.argv[m+8]


OUTPUT_FOLDER = htmlfile.split('.')[0]+'_files'   #/media/Work/galaxy-globus-crdata/database/files/000/dataset_66_files

if not os.path.exists(OUTPUT_FOLDER):
	os.makedirs(OUTPUT_FOLDER)
else:
	pass


parameter = " --compute-events-psi "+labelname+" "+countfile+" --event-type "+event_type+" --read-len "+read_len+" --overhang-len "+overhang_len+" --output-dir "+OUTPUT_FOLDER

if not events_file == "":
    parameter = parameter + " --events-file "+ events_file

if not settings_file == "":
    parameter = parameter + " --settings-filename "+ settings_file

if not job_name == "":
    parameter = parameter + " --job-name "+ job_name

#print parameter


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
        res.append(galhtmlattr % ('Compute Psi values for events',timenow()))
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

os.chdir(galaxyhome + "/tools/miso/")     #uncomment when running on cluster
#os.chdir("/media/Work/galaxy-proteonics/tools/miso/")    #for local running

command="python2.7 ./run_events_analysis.py "+parameter

#command = "python2.7 /media/Work/galaxy-proteonics/tools/miso/run_events_analysis.py  --compute-events-psi se-sample /media/Work/galaxy-proteonics/tools/miso/test-data/se-counts/se_test.counts --output-dir /media/Work/galaxy-proteonics/tools/miso/test-output/SE-output --read-len 35 --overhang-len 4 --event-type SE"

print command

subprocess.call(command,shell=True)

htmloutput(htmlfile,OUTPUT_FOLDER)


