#!/usr/bin/python

CHUNK_SIZE = 2**20 #1mb

import argparse, os, shutil, subprocess, sys, tempfile, shlex, vcf, pysam

parser = argparse.ArgumentParser(description='')
parser.add_argument ( '-t', dest='target_dir', help='target directory')
parser.add_argument ( '-i', dest='input_dir', help='input directory')
parser.add_argument ( '-o', dest='output', help='output log file' )

def __main__():
    args = parser.parse_args()
    output_dir = tempfile.mkdtemp(prefix="homer-");
    
    target_dirs = " ".join(args.target_dir)
    input_dirs = " ".join(args.input_dir)

    homer="/mnt/galaxyTools/tools/homer/4.9"
    homer_bin="/mnt/galaxyTools/tools/homer/4.9/bin"
    env="export PATH=%s:%s:$PATH;" % (homer, homer_bin)
    cmd="%s perl %s/getDifferentialPeaksReplicates.pl -t %s -i %s > %s/outputPeaks.txt" % (env, homer_bin, args.target_dir, args.input_dir, output_dir)

    print("cmd: %s" % cmd)

    stdout = tempfile.NamedTemporaryFile( prefix="getdifferentialpeaksreplicates-stdout-", dir=output_dir )
    stderr = tempfile.NamedTemporaryFile( prefix="getdifferentialpeaksreplicates-stderr-", dir=output_dir )  

    proc = subprocess.Popen( args=cmd, stdout=stdout, stderr=stderr, shell=True)
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
    ##shutil.copy("%s/
    stderr.close()
    stdout.close() 
    shutil.copy('%s/outputPeaks.txt' % output_dir, args.output)

if __name__=="__main__":
	__main__()
