#!/usr/bin/env python
#By, Guruprasad Ananda.

from galaxy import eggs
import sys, os, glob

def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit()

def _sort_input_fastq(input_dir):
    forward = []
    reverse = []
    for infile in sorted(os.listdir(input_dir)):
        # assume files are fastq
        if infile.endswith(".fastq.gz"):
            if "_R1_" in infile or "_1_" in infile or "_1.fastq.gz" in infile:
                forward.append("%s/%s" % (input_dir, infile))
            elif "_R2_" in infile or "_2_" in infile or "_2.fastq.gz" in infile:
                reverse.append("%s/%s" % (input_dir, infile))
    return (forward, reverse)
    
def main():
    outfile1 = sys.argv[1]
    outfile2 = sys.argv[2]
    input_dir = sys.argv[3]
    
    input_files = sorted(glob.glob("%s/*.fastq.gz" % input_dir))
    try:
        fout = open(sys.argv[1],'w')
    except:
        stop_err("Output file cannot be opened for writing.")

    forward, reverse = _sort_input_fastq(input_dir)        
    try:
        fin = open(forward[0],'r')
        rin = open(reverse[0],'r')
    except:
        stop_err("Input file cannot be opened for reading.")
    
    cmdline = "cat %s " % (" ".join(forward))
    cmdline = cmdline + ">" + outfile1
    try:
        os.system(cmdline)
    except:
        stop_err("Error encountered with forward fastq cat.")
        
    cmdline = "cat %s " % (" ".join(reverse))
    cmdline = cmdline + ">" + outfile2
    try:
        os.system(cmdline)
    except:
        stop_err("Error encountered with reverse fastq cat.")

if __name__ == "__main__": main()
