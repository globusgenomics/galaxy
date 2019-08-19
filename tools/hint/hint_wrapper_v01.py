#!/usr/bin/python

import optparse, os, shutil, sys, tempfile, glob, shlex, vcf, pysam, tarfile
from subprocess import *
import subprocess

CHUNK_SIZE = 2**20 #1mb

def cleanup_before_exit( tmp_dir ):
    if tmp_dir and os.path.exists( tmp_dir ):
        shutil.rmtree( tmp_dir )

def __main__():
    parser = optparse.OptionParser()
    parser.add_option('','--input', dest="inputF", action='store', type="string", default=None, help='')
    parser.add_option('','--region', dest="regionF", action='store', type="string", default=None, help='')
    parser.add_option( '', '--output', dest='outputF', action='store', type="string", default=None, help='')
    parser.add_option( '', '--out-dir', dest='output_dir', action='store', type="string", default=None, help='If specified, the output directory for extra files.' )

    (options, args) = parser.parse_args()

    if not os.path.exists(options.output_dir): 
        os.mkdir(options.output_dir)
    input_dir = "%s/input" % options.output_dir
    output_dir = "%s/output" % options.output_dir
    config_dir = "%s/config" % options.output_dir
           
    if not os.path.exists(output_dir):
        os.mkdir(input_dir)
        os.mkdir(output_dir)
        os.mkdir(config_dir)
    # region input file
    linked_bed_name = "%s/regions.bed" % input_dir
    if not os.path.exists(linked_bed_name):
        os.symlink(options.regionF, linked_bed_name)
#    shutil.copy(options.regionF, linked_bed_name)
#    samfile = pysam.AlignmentFile(options.inputF, 'rb')
#    sn = samfile.header['RG'][0]['SM']
    linked_bam_name="%s/sample.bam" % input_dir
    if not os.path.exists(linked_bam_name):  
        os.symlink(options.inputF, linked_bam_name)
        pysam.index(linked_bam_name)
    
    input_config = open("%s/input.txt" % config_dir, 'w')
    input_config.write("name\ttype\tfile\tdata\tgroup\n")
    input_config.write("HS2\tregions\t%s\tHS\tDU_K562_HINT\n" % linked_bed_name) 
    input_config.write("DNase\treads\t%s\tDNASE\tDU_K562_HINT\n" % linked_bam_name)
    input_config.close()

    #hint_cmd = "rgt-hint --output-location %s/ %s/input.txt"  % (output_dir, config_dir)
    hint_cmd = "rgt-hint --output-location %s/ %s/input.txt"  % (output_dir, config_dir)
    print "hint cmd:%s" % hint_cmd
    stdout = tempfile.NamedTemporaryFile( prefix="hint-stdout-", dir=options.output_dir)
    stderr = tempfile.NamedTemporaryFile( prefix="hint-stderr-", dir=options.output_dir)
    proc = subprocess.Popen( args=hint_cmd, stdout=stdout, stderr=stderr, shell=True, cwd=options.output_dir )
    return_code = proc.wait()

    if return_code:
        stderr_target = sys.stderr
    else:
        stderr_target = sys.stdout
        stderr.flush()
        stderr.seek(0)
    while True:
        chunk = stderr.read( CHUNK_SIZE )
        if chunk:
            stderr_target.write( chunk )
        else:
            break
    stderr.close()
    stdout.close()
                                                                                    
    # copy files to final output locations
    shutil.copy('%s/DU_K562_HINT.bed' % output_dir, options.outputF)
    cleanup_before_exit( options.output_dir )

if __name__=="__main__": __main__()

