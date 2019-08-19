#!/usr/bin/env python
import gzip, sys, os, commands, string, subprocess, shutil, glob

galaxyhome=os.environ.get('GALAXY_HOME')
#parser = optparse.OptionParser(description="")

n = len(sys.argv)
for i in range(0,n):
    print sys.argv[i]


IndexDIR = sys.argv[1]
Index_FOLDER = IndexDIR.split('.')[0]

output = " -o -bam " + sys.argv[2]

if "bag" not in sys.argv[3]:

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
        inputFASTA = "paired " + Index_FOLDER + " -mrl " + str(min_size) + " -compressedFastq " + sys.argv[3] + " " + sys.argv[4]
    else:
        inputFASTA = "single " + Index_FOLDER + " -mrl " + str(min_size) + " -compressedFastq " + sys.argv[3]
elif sys.argv[3] == "bag":
    # this is a bag, figure out where things are.
    files = glob.glob("%s/*.fastq.gz" % sys.argv[4])
    print files
    # get read length for inputFASTA
    handle = gzip.open(files[0], 'rb')
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

    if len(files) == 1:
        inputFASTA = "single " + Index_FOLDER + " -xf 2.0 -mrl " + str(min_size) + " -compressedFastq " + files[0]
    else:
        if "_1" in files[0]:
            inputFASTA = "paired " + Index_FOLDER + " -mrl " + str(min_size) + " -compressedFastq " + files[0] + " " + files[1]
        else:
            inputFASTA = "paired " + Index_FOLDER + " -mrl " + str(min_size) + " -compressedFastq " + files[1] + " " + files[0]

elif "bdbag" in sys.argv[3]:
    # This was and ENCODE BDBAG
    # get the metadata.tsv file and figure out if it's paired or single by using the ID given
    metadata_path = glob.glob("%s/*/*/data/metadata.tsv" % sys.argv[4])[0]
    fh = open(metadata_path, "r")
    sample_lines = []
    sample = sys.argv[3].split("~")[1]
    for line in fh:
        if sample in line:
            sample_lines.append(line)
    files = []
    if len(sample_lines) == 1:
        # single end
        fastq = glob.glob("%s/*/*/data/%s.fastq.gz" % (sys.argv[4], sample))[0]
        files.append(fastq)
    elif len(sample_lines) == 2:
        values1 = sample_lines[0].split("\t")
        values2 = sample_lines[1].split("\t")
        if values1[33] == "1":
            files.append(glob.glob("%s/*/*/data/%s.fastq.gz" % (sys.argv[4],values1[0]))[0])
            files.append(glob.glob("%s/*/*/data/%s.fastq.gz" % (sys.argv[4],values2[0]))[0])
        elif values1[33] == "2":
            files.append(glob.glob("%s/*/*/data/%s.fastq.gz" % (sys.argv[4],values2[0]))[0])
            files.append(glob.glob("%s/*/*/data/%s.fastq.gz" % (sys.argv[4],values1[0]))[0])

    # get read length for inputFASTA
    handle = gzip.open(files[0], 'rb')
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

    if len(files) == 1:
        inputFASTA = "single " + Index_FOLDER + " -t 32 -xf 2.0 -mrl " + str(min_size) + " -compressedFastq " + files[0]
    else:
        inputFASTA = "paired " + Index_FOLDER + " -t 32 -mrl " + str(min_size) + " -compressedFastq " + files[0] + " " + files[1]

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

