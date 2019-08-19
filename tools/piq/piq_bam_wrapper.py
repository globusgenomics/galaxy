#!/usr/bin/python

import optparse, os, re, shutil, sys, tempfile, glob, shlex, vcf, pysam, tarfile, csv, operator
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
    piq_dir_path="/mnt/galaxyTools/tools/piq/1.3/thashim"
    parser = optparse.OptionParser()
    parser.add_option('', '--input', dest="inputF", action='store', type="string", default=None, help='')
    parser.add_option( '', '--output', dest='outputF', action='store', type="string", default=None, help='')

    (options, args) = parser.parse_args()

    ephemeral_dir = tempfile.mkdtemp(prefix="piq-")
    ephemeral_tmp_dir = "%s/tmp" % ephemeral_dir;
    if not os.path.exists(ephemeral_tmp_dir):
        os.mkdir(ephemeral_tmp_dir)

    piq_cmd2 = "Rscript %s/bam2rdata.r %s/common.r %s/bam.RData %s; " % (piq_dir_path, piq_dir_path, ephemeral_tmp_dir, options.inputF)

    print "piq_cmd2:%s" % piq_cmd2
    stdout = tempfile.NamedTemporaryFile( prefix="piq_cmd2-stdout-", dir=ephemeral_tmp_dir )
    stderr = tempfile.NamedTemporaryFile( prefix="piq_cmd2-stderr-", dir=ephemeral_tmp_dir )
    proc = subprocess.Popen( args=piq_cmd2, stdout=stdout, stderr=stderr, shell=True, cwd=ephemeral_dir)
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

    shutil.copy('%s/bam.RData' % ephemeral_tmp_dir, options.outputF)

if __name__=="__main__": __main__()

