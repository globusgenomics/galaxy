import sys
import pyDNase
import pyDNase.footprinting as fp
import subprocess
import os

if __name__ == "__main__":

    #There should be two input files:
    #1. original .bam file
    #2. fseq output .bed file
    #3. folder needs to contain .bam.bai file for reference

    input_bed_name = sys.argv[1] 
    input_bam_name = sys.argv[2]
    input_bam_index = sys.argv[3]
    output_file_name = sys.argv[4]
    pvalue_cutoff = sys.argv[5]

    # create links for input bam and indexes
    os.symlink(input_bam_name, "input.bam")
    os.symlink(input_bam_index, "input.bam.bai")

    #wellington 
    regions = pyDNase.GenomicIntervalSet(input_bed_name)
    reads = pyDNase.BAMHandler("input.bam")
    footprinter = fp.wellington(regions[0], reads)
    footprints = footprinter.footprints(withCutoff=pvalue_cutoff)
    with open(output_file_name, "w") as resultout:
        resultout.write(str(footprints))


