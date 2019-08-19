#!/usr/bin/env python
import gzip, sys, os, commands, string, subprocess, shutil

galaxyhome=os.environ.get('GALAXY_HOME')


n = len(sys.argv)
for i in range(0,n):
    print sys.argv[i]


IndexDIR = sys.argv[1]
Index_FOLDER = IndexDIR.split('.')[0]

output = " -o -bam " + sys.argv[2]

# get read length for inputFASTA
handle = gzip.open(sys.argv[3], 'rb')
r_name = handle.readline()
read = handle.readline().rstrip("\n")
length = len(read)
handle.close()

min_size = 20
if length >= 50:
    min_size = 40
elif length >= 33:
    min_size=26
elif length >= 27:
    min_size = 22
elif length >= 20:
    min_size = 20

if len(sys.argv)>4:
    inputFASTA = "paired " + Index_FOLDER + " -mrl " + str(min_size) + " -compressedFastq " + sys.argv[3] + "  " + sys.argv[4]
else:
    inputFASTA = "single " + Index_FOLDER + " -mrl " + str(min_size) + " -compressedFastq " + sys.argv[3]


##uncomment when running on VM
#os.chdir(galaxyhome + "/tools/snap/")     #uncomment when running on cluster
#os.chdir("/media/Work/galaxy-proteonics/tools/snap/")    #for local running

command="snap " + inputFASTA + output

print command

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

