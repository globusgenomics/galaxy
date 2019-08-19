#!/usr/bin/env python
#By, Guruprasad Ananda.

import sys, os, glob

def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit()
    
def get_paired_reads(input_dir):
    # the format of the file names must be:
    # sampleName_L00#_R#_001.fastq.gz
    forward_reads = sorted(glob.glob("%s/*_R1_*" % input_dir))
    reverse_reads = sorted(glob.glob("%s/*_R2_*" % input_dir))
    return forward_reads, reverse_reads

def main():
    outfile1 = sys.argv[1]
    outfile2 = sys.argv[2]
    indir = sys.argv[3]
    
    (f_reads, r_reads) = get_paired_reads(indir)

    cmdline1 = "cat %s > %s" %(" ".join(f_reads), outfile1)
    print cmdline1
    cmdline2 = "cat %s > %s" %(" ".join(r_reads), outfile2)
    print cmdline2

    try:
        os.system(cmdline1)
        os.system(cmdline2)
    except:
        stop_err("Error encountered with cat.")
        
if __name__ == "__main__": main()
