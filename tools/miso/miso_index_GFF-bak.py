#!/usr/bin/env python

import sys, os, commands, string, time, subprocess

galaxyhome=os.environ.get('GALAXY_HOME')

INDEX_GFF = sys.argv[1]
htmlout = sys.argv[2]

outputdir = htmlout.split('.')[0]+'_files'   #/media/Work/galaxy-globus-crdata/database/files/000/dataset_66_files

if not os.path.exists(outputdir):
	os.makedirs(outputdir)
else:
	pass

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
        res.append(galhtmlattr % ('Index GFF',timenow()))
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


os.chdir(galaxyhome + "/tools/miso/")     #uncomment when running on cluster
#os.chdir("/media/Work/galaxy-proteonics/tools/miso/")    #for local running

if INDEX_GFF != '':
	subprocess.call(["python","./misopy/index_gff.py","--index",INDEX_GFF, outputdir])
else:
	print "lack GFF file"


htmloutput(htmlout,outputdir)


