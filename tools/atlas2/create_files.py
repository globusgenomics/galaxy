#!/usr/bin/env python
#arodri7

"""
A wrapper script for running swift on globus-galaxy
"""

import sys, optparse, os, tempfile, subprocess, shutil
    
#Parse Command Line     
parser = optparse.OptionParser()
parser.add_option( '-q', '--qty', dest='number_files', type="int", help='Input directory where all BAM file are located' )
parser.add_option( '-o', '--outputdir', dest='outputdir', type="string", help='Output directory to store VCF files' )
parser.add_option( '-d', '--dummy', dest='dummy_out', type="string")
(options, args) = parser.parse_args()
os.makedirs(options.outputdir)

for i in range(1,options.number_files):
    filename = "%s/file%s.tmp" % (options.outputdir, i)
    open(filename, 'a').close()

stdout = open( options.dummy_out, 'w' )
stdout.write("Wrote files")
stdout.close
