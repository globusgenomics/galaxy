#!/usr/bin/env python

import sys, os, commands, string, time, subprocess

galaxyhome=os.environ.get('GALAXY_HOME')

'''
n = len(sys.argv)
for i in range(0,n):
    print sys.argv[i]
'''

input_gene = sys.argv[1]

htmlfile = sys.argv[2]

input_gene_FOLDER = input_gene.split('.')[0]+'_files'

OUTPUT_FOLDER = htmlfile.split('.')[0]+'_files'   #/media/Work/galaxy-globus-crdata/database/files/000/dataset_66_files

if not os.path.exists(OUTPUT_FOLDER):
	os.makedirs(OUTPUT_FOLDER)
else:
	pass


#command = "python ./misopy/run_events_analysis.py --compute-genes-psi "+input_GFF_FOLDER+" "+input_BAM+" --output-dir "+OUTPUT_FOLDER+" --read-len "+read_len
#print command

parameter = " --summarize-samples "+input_gene_FOLDER+" "+OUTPUT_FOLDER
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
        res.append(galhtmlattr % ('Summarize Samples',timenow()))
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

command="python2.7 ./run_miso.py "+parameter

#command = "python2.7 ./run_events_analysis.py  --compute-genes-psi /media/Work/galaxy-proteonics/database/files/000/dataset_150_files /media/Work/galaxy-proteonics/database/files/000/dataset_144_files/miso_BamBai.bam --output-dir /media/Work/galaxy-proteonics/database/files/000/dataset_154_files --read-len 36 --job-name job"

print command

subprocess.call(command,shell=True)

htmloutput(htmlfile,OUTPUT_FOLDER)


