#!/usr/bin/python

import optparse, os, shutil, sys, tempfile, glob, pysam
from subprocess import *
import subprocess

def run_cmd ( cmd , descriptor):
    stderr_name = tempfile.NamedTemporaryFile( prefix = "cmd_stderr" ).name
    proc = Popen( args=cmd, shell=True, stderr=open( stderr_name, 'wb' ) )
    exit_code = proc.wait()
    if exit_code:
        for line in open( stderr_name ):
            print >> sys.stderr, line
            os.unlink( stderr_name ) #clean up
            raise Exception( "Error running command: %s " % descriptor )
        os.unlink( stderr_name ) #clean up

def __main__():
    parser = optparse.OptionParser()
    parser.add_option( '', '--bams', dest='bam_list', type="string", help='File list of BAMS' )
    parser.add_option( '-g', '--group-field', dest='groupField', help='RG field to use for grouping of BAMs' )
    parser.add_option( '-s', '--keyword-separator', dest='separator', help='optional: if value in RG field contains more than just the group ID, then we need a separator' )
    parser.add_option( '-f', '--keyword-field-order', dest='fieldOrder', help='optional: if value in RG field contains more than just the group ID, then we need the order where group name appears in' )
    parser.add_option( '-o', '--output', dest='outputfile', type="string", help='output File of groups' )
    parser.add_option( '-d', '--output-dir', dest='outputdirectory', type="string", help='output directory of groups' )
    (options, args) = parser.parse_args()

    output_dir = options.outputdirectory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    ## get list of input BAM files
    directory = None
    if os.path.isdir(options.bam_list):
        directory = options.bam_list
    else:
        fh = open(options.bam_list)
        while line in fh:
            directory = line.rstrip("\n")
            break

    # group BAM list by selected field and store in dictionary
    bam_groups = {}
    for bam in glob.glob("%s/*.bam" % directory):
        samfile = pysam.AlignmentFile(bam,'rb')
        value = samfile.header['RG'][0][options.groupField]
        if options.separator is not None:
            value = value.split(options.separator)[int(options.fieldOrder)-1]
        if value in bam_groups:
            bam_groups[value].append(bam)
        else:
            bam_groups[value] = [bam]

    # write each group into a separate file in the output directory 
    # and count them and write group name and count to the output file
    fh_file = open(options.outputfile, "w")
    fh_file.write("Groups generated on BAM input from directory: %s\n" % directory)
    fh_file.write("Separating on RG field %s\n" % options.groupField)
    if options.separator is not None:
        fh_file.write("Which includes separator %s and using field %d" % (options.separator, int(options.fieldOrder)-1))
    fh_file.write("\n\nGROUP_NAME\tBAM_COUNT\n")
    

    for group in bam_groups:
        fh = open("%s/%s.txt" % (options.outputdirectory, group), "w")
        fh.write("\n".join(sorted(bam_groups[group])))
        fh.close()
        fh_file.write("%s\t%s\n" % (group, len(bam_groups[group])))
    fh_file.close()

if __name__=="__main__":
	__main__()
