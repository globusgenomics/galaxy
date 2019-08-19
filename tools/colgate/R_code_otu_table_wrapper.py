#!/usr/bin/python

import optparse, os, re, shutil, sys, tempfile, glob
from subprocess import *
import subprocess
import multiprocessing
from joblib import Parallel, delayed
from itertools import izip
from os.path import isfile, join; from os import listdir

CHUNK_SIZE = 2**20 #1mb

def __main__():
    Colgate_Dir = "/mnt/galaxyTools/tools/colgate/0.1/"
    parser = optparse.OptionParser()
    parser.add_option('-i', '--shared_table', dest="shared_table", action='store', type="string", default=None, help='')
    parser.add_option('-t', '--taxonomy', dest="taxonomy", action='store', type="string", default=None, help='')
    parser.add_option( '-o', '--output', dest='output', action='store', type="string", default=None, help='')
    parser.add_option( '', '--out-dir', dest='output_dir', action='store', type="string", default=None, help='If specified, the output directory for extra files.' )
    (options, args) = parser.parse_args()

    
    if not os.path.exists(options.output_dir):
        output_dir=options.output_dir
        os.mkdir(output_dir)

    cmd = "Rscript %s/R_code_otu_table.R -i %s -t %s -o %s/otu.txt" % (Colgate_Dir, options.shared_table, options.taxonomy, output_dir)

    print "cmd:%s" % cmd

    stdout = tempfile.NamedTemporaryFile( prefix="cmd-stdout-", dir= output_dir)
    stderr = tempfile.NamedTemporaryFile( prefix="cmd-stderr-", dir= output_dir)
    proc = subprocess.Popen( args=cmd, stdout=stdout, stderr=stderr, shell=True, cwd=output_dir)
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
    shutil.copy("%s/otu.txt" % output_dir, options.output)

if __name__=="__main__": __main__()
