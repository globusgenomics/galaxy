#!/usr/bin/python

import optparse, os, re, shutil, sys, tempfile, glob, shlex, vcf, pysam, tarfile, csv, operator
from subprocess import *
import subprocess
import multiprocessing
from joblib import Parallel, delayed
from itertools import izip

NCPU = multiprocessing.cpu_count()
CHUNK_SIZE = 2**20 #1mb

def parallel_jobs(motif, piq_dir_path, tmp_dir, output_dir, motif_file_name):
    cmd = "Rscript %s/pwmmatch.exact.r %s/common.r %s %d %s/; " % (piq_dir_path, piq_dir_path, motif_file_name, motif, tmp_dir)
    print "\ncmd: %s\n" % cmd
    stdout = tempfile.NamedTemporaryFile( prefix="piq-stdout-", dir=output_dir )
    stderr = tempfile.NamedTemporaryFile( prefix="piq-stderr-", dir=output_dir )
    proc = subprocess.Popen( args=cmd, stdout=stdout, stderr=stderr, shell=True, cwd=output_dir )
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
    
def __main__():
    piq_dir_path="/mnt/galaxyTools/tools/piq/1.3/thashim"
    parser = optparse.OptionParser()
    parser.add_option('', '--input', dest="inputF", action='store', type="string", default=None, help='')
    parser.add_option( '', '--output', dest='outputF', action='store', type="string", default=None, help='')
    parser.add_option( '', '--out-dir', dest='output_dir', action='store', type="string", default=None, help='If specified, the output directory for extra files.' )

    (options, args) = parser.parse_args()
    
    ##if not os.path.exists(options.output_dir):
    ##    os.mkdir(options.output_dir)

    ephemeral_dir = tempfile.mkdtemp(prefix="piq-")
    ephemeral_tmp_dir = "%s/tmp" % ephemeral_dir;
    if not os.path.exists(ephemeral_tmp_dir):
        os.mkdir(ephemeral_tmp_dir)

#    maxMotif=1
#    maxMotif=519
#    Parallel(n_jobs=NCPU)(delayed(parallel_jobs)(i, piq_dir_path, tmp_dir, options.inputF, out_dir, options.output_dir) for i in range(1,maxMotif+1))
    F = open(options.inputF,"r")
    maxMotif=0
    for line in F:
        if re.search(">",line):
            maxMotif=maxMotif+1
        else:
            continue
    F.close()
    print "max number of motif:%s" % maxMotif
    Parallel(n_jobs=NCPU)(delayed(parallel_jobs)(i, piq_dir_path, ephemeral_tmp_dir, ephemeral_dir, options.inputF) for i in range(1,maxMotif+1))
    
    shutil.copytree(ephemeral_tmp_dir, options.output_dir)

    # concatenerate calls.all and RC.calls.all for csv and bed files
    #cleanup_before_exit( options.output_dir )

if __name__=="__main__": __main__()

