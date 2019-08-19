#!/usr/bin/python

import optparse, os, re, shutil, sys, tempfile, glob
from subprocess import *
import subprocess
import multiprocessing
from joblib import Parallel, delayed
from itertools import izip
from os.path import isfile, join; from os import listdir

CHUNK_SIZE = 2**20 #1mb

def cleanup_before_exit( tmp_dir ):
    if tmp_dir and os.path.exists( tmp_dir ):
        shutil.rmtree( tmp_dir )

def __main__():
    Colgate_Dir = "/mnt/galaxyTools/tools/colgate/0.1/" 
    parser = optparse.OptionParser()
    parser.add_option('', '--otu_table', dest="otutable", action='store', type="string", default=None, help='')
    parser.add_option('', '--map_file', dest="mapfile", action='store', type="string", default=None, help='')
    parser.add_option( '', '--plot1', dest='plot1', action='store', type="string", default=None, help='')
    parser.add_option( '', '--plot2', dest='plot2', action='store', type="string", default=None, help='')
    parser.add_option( '', '--plot3', dest='plot3', action='store', type="string", default=None, help='')
    parser.add_option( '', '--plot4', dest='plot4', action='store', type="string", default=None, help='')
    parser.add_option( '', '--plot5', dest='plot5', action='store', type="string", default=None, help='')
    parser.add_option( '', '--plot6', dest='plot6', action='store', type="string", default=None, help='')
    parser.add_option( '', '--plot11', dest='plot11', action='store', type="string", default=None, help='')
    parser.add_option( '', '--plot7', dest='plot7', action='store', type="string", default=None, help='')
    parser.add_option( '', '--plot8', dest='plot8', action='store', type="string", default=None, help='')
    parser.add_option( '', '--plot9', dest='plot9', action='store', type="string", default=None, help='')
    parser.add_option( '', '--plot10', dest='plot10', action='store', type="string", default=None, help='')
    parser.add_option( '', '--plot12', dest='plot12', action='store', type="string", default=None, help='')
    parser.add_option( '', '--table1', dest='table1', action='store', type="string", default=None, help='')
    parser.add_option( '', '--table2', dest='table2', action='store', type="string", default=None, help='')
    parser.add_option( '', '--out-dir', dest='output_dir', action='store', type="string", default=None, help='If specified, the output directory for extra files.' )

    (options, args) = parser.parse_args()
#    ephemeral_dir = tempfile.mkdtemp(prefix="colgate-")
#    ephemeral_tmp_dir = "%s/tmp" % ephemeral_dir;
#    if not os.path.exists(ephemeral_tmp_dir):
#        os.mkdir(ephemeral_tmp_dir)

    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)

    PI_cmd = "Rscript %s/R_code_PI.R --otu_table %s --map_file %s --plot1 %s --plot2 %s --plot3 %s --plot4 %s --plot5 %s --plot6 %s --plot7 %s --plot8 %s --plot9 %s --plot10 %s --plot11 %s --plot12 %s --table1 %s --table2 %s" % (Colgate_Dir, options.otutable, options.mapfile, options.plot1, options.plot2, options.plot3, options.plot4, options.plot5, options.plot6, options.plot7, options.plot8, options.plot9, options.plot10, options.plot11, options.plot12, options.table1, options.table2)

    print "PI_cmd:%s" % PI_cmd
    #stdout = tempfile.NamedTemporaryFile( prefix="PI_cmd-stdout-", dir=ephemeral_tmp_dir )
    #stderr = tempfile.NamedTemporaryFile( prefix="PI_cmd-stderr-", dir=ephemeral_tmp_dir )
    #proc = subprocess.Popen( args=PI_cmd, stdout=stdout, stderr=stderr, shell=True, cwd=ephemeral_dir)
    stdout = tempfile.NamedTemporaryFile( prefix="PI_cmd-stdout-", dir= options.output_dir)
    stderr = tempfile.NamedTemporaryFile( prefix="PI_cmd-stderr-", dir= options.output_dir)
    proc = subprocess.Popen( args=PI_cmd, stdout=stdout, stderr=stderr, shell=True, cwd=options.output_dir)
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

#    shutil.copy('%s/plot1.png' % ephemeral_tmp_dir, options.plot1)
#    shutil.copy('%s/plot1.png' % options.output_dir, options.plot1)
if __name__=="__main__": __main__()
