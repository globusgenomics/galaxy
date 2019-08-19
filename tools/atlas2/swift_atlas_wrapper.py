#!/usr/bin/env python
#arodri7

"""
A wrapper script for running swift on globus-galaxy
"""

import sys, optparse, os, tempfile, subprocess, shutil
    
#Parse Command Line     
parser = optparse.OptionParser()
parser.add_option( '-s', '--sites.file', dest='sites_file', type="string", help='Swift site file' )
parser.add_option( '-t', '--tc.file', dest='tc_file', type="string", help='Swift TC file' )
parser.add_option( '-b', '--bamdir', dest='bamdir', type="string", help='Input directory where all BAM file are located' )
parser.add_option( '-o', '--outputdir', dest='outputdir', type="string", help='Output directory to store VCF files' )
parser.add_option( '-f', '--fasta', dest='fasta', type="string", help='Fasta file' )
parser.add_option( '-d', '--dummy', dest='dummy_out', type="string")
(options, args) = parser.parse_args()

swift_file = args[0]

#cmd = 'sleep 10'
cmd = "swift -sites.file %s -tc.file %s %s -bamdir=%s -outputdir=%s -fasta=%s 2>&1" % (options.sites_file, options.tc_file, swift_file, options.bamdir, options.outputdir, options.fasta)
print cmd
proc = subprocess.Popen( args=cmd, shell=True )
returncode = proc.wait()

stdout = open( options.dummy_out, 'w' )
#stdout.write(options.outputdir)
stdout.write("BAM\t%s\nVCF\t%s\nREF\t%s\n") % (options.bamdir, options.outputdir, options.fasta)
stdout.close

