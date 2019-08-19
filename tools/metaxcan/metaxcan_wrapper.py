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
    parser.add_option( '', '--covariance', dest='covar', action='store', type="string", help='' )
    parser.add_option( '', '--beta_column', dest='beta', action='store', type="string", help='' )
    parser.add_option( '', '--pvalue_column', dest='pval', action='store', type="string", help='' )
    parser.add_option( '', '--gwas_folder', dest='gwasF', action='store', type="string", help='' )
    parser.add_option( '', '--weight_db_path', dest='weight_db_path', action='store', type="string", help='' )
    parser.add_option( '', '--output_file', dest='outputF', action='store', type="string", default=None, help='')
    parser.add_option( '', '--out-dir', dest='output_dir', action='store', type="string", default=None, help='If specified, the output directory for extra files.' )

    (options, args) = parser.parse_args()

    # worker node tmp directory where all will be processed
    beta_dir= "%s/beta" % options.output_dir
    if os.path.exists(beta_dir):
            shutil.rmtree(beta_dir)
    os.mkdir(beta_dir)
    # create input files from the input_dir
    pattern=".*gz"
    metaxcan_cmd = "python /mnt/galaxyTools/tools/metaxcan/09072016/software/MetaXcan.py --beta_folder %s --weight_db_path %s --covariance %s --gwas_folder %s --gwas_file_pattern %s --compressed --beta_column %s --pvalue_column %s --output_file %s/results.csv"  % (beta_dir, options.weight_db_path, options.covar, options.gwasF, pattern, options.beta, options.pval, options.output_dir)
    print metaxcan_cmd

    stdout = tempfile.NamedTemporaryFile( prefix="metaxcan-stdout-", dir=options.output_dir )
    stderr = tempfile.NamedTemporaryFile( prefix="metaxcan-stderr-", dir=options.output_dir )
    proc = subprocess.Popen( args=metaxcan_cmd, stdout=stdout, stderr=stderr, shell=True, cwd=options.output_dir )
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
    shutil.copy('%s/results.csv' % options.output_dir, options.outputF) 
    cleanup_before_exit( options.output_dir )

if __name__=="__main__": __main__()
