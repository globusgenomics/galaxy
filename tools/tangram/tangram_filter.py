#!/usr/bin/env python
"""
Run tangram_filter.pl
usage: tangram_filter.py [options]
"""

import subprocess
import glob, fnmatch, time, re, optparse, os, sys, tempfile, shutil

CHUNK_SIZE = 2**20 #1mb

"""
   tangram_filter.py
      --vcf $vcf_input
      --msk $msk_file
      --out $output
      --type $sv_type
      --window $window_size
      #if $advanced.advanced_select == "yes":
          #if $advanced.rpfp!= 2:
              --rpf $advanced.rpf
          #end if
          #if $advanced.srf != 2:
              --srf $advanced.srf
          #end if
      #end if

"""

def run_cmd ( cmd , wd_tmpdir, descriptor):
    print cmd
    stderr = tempfile.NamedTemporaryFile( prefix = "cmd_stderr" )
    proc = subprocess.Popen( args=cmd, stderr=stderr, shell=True, cwd=wd_tmpdir )

    exit_code = proc.wait()

    if exit_code:
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


def __main__():
    print "\n\nStart time:"
    print time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))

    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '', '--vcf', dest='vcf', help="Input VCF" )
    parser.add_option( '', '--msk', dest='maskFile', help='Mask file' )
    parser.add_option( '', '--type', dest='type', help='SV type' )
    parser.add_option( '', '--window', dest='window', help='Window size around the variant' )
    parser.add_option( '', '--rpf', dest='rpf', help='RPF number' )
    parser.add_option( '', '--srf', dest='srf', help="SRF number" )
    parser.add_option( '', '--out', dest='output', help='The output.' )
    ( options, args ) = parser.parse_args()

    
    ## Create a temporary directory where output will be written
    wd_tmpdir = tempfile.mkdtemp( prefix='tmp-tangram-filter-' )

    ## create the mask input file containing a list of masked files as input
    maskInput_path = "%s/mask.txt" % (wd_tmpdir)
    maskInput = open(maskInput_path, "w")
    maskInput.write("%s\t%s\t%s\n" % ("AL", options.window, options.maskFile) )
    maskInput.close()

    ## create the command line
    cmd = "tangram_filter.pl --vcf %s --msk %s " % (options.vcf, maskInput_path )

    if options.rpf:
        cmd += "-rpf %s " % (options.rpf)
    if options.srf:
        cmd += "-srf %s " % (options.srf)

    cmd += " > %s" % options.output
    log_msg = "Running Tangram filter"

    ## run the command
    run_cmd(cmd, wd_tmpdir, log_msg)

    ## clean up tmp directories
    shutil.rmtree(wd_tmpdir)

    print "\n\nEnd time:"
    print time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))

if __name__=="__main__": __main__()

