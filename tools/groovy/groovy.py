#!/usr/bin/env python

import os, sys, commands, string, subprocess, shutil

galaxyhome=os.environ.get('GALAXY_HOME')

output = sys.argv[1]

##uncomment when running on VM
os.chdir(galaxyhome + "/tools/groovy/")     #uncomment when running on cluster
#os.chdir("/media/Work/galaxy-proteonics/tools/groovy/")    #for local running

command="groovy -cp GLSGeneusRestApiUtils.groovy UpdateProcessUDF_from_Galaxy.groovy"

file = open(output, 'w')

proc = subprocess.Popen( args=command, shell=True, stdout=file, stderr=subprocess.PIPE )

