#!/usr/bin/env python
import sys, os, commands, string, subprocess, shutil

galaxyhome=os.environ.get('GALAXY_HOME')


n = len(sys.argv)
for i in range(0,n):
    print sys.argv[i]


inputFASTA = sys.argv[1]

if sys.argv[2]!="none":
    seedSize = " " + sys.argv[2]
else:
    n = ""

outputIndex = sys.argv[3]

outputIndex_FOLDER = outputIndex.split('.')[0]+'_files'

if not os.path.exists(outputIndex_FOLDER):
	os.makedirs(outputIndex_FOLDER)
else:
	pass


##uncomment when running on VM
os.chdir(galaxyhome + "/tools/snap/")     #uncomment when running on cluster
#os.chdir("/media/Work/galaxy-proteonics/tools/snap/")    #for local running

command="./snap index " + inputFASTA + " " + outputIndex_FOLDER

print command

#proc = subprocess.Popen( args=command, shell=True, stderr=subprocess.PIPE )
try:
    proc = subprocess.Popen( args=command, shell=True, stderr=subprocess.PIPE )
    returncode = proc.wait()
    stderr = ''
    buffsize = 1048576
    try:
        while True:
            stderr += proc.stderr.read(buffsize)
            if not stderr or len(stderr)% buffsize != 0:
                break
    except OverflowError:
            pass
    if returncode != 0:
        raise Exception, stderr
except Exception,e:
    print str(e)
