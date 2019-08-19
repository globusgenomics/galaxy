#!/usr/bin/python

CHUNK_SIZE = 2**20 #1mb

import argparse, os, shutil, subprocess, sys, tempfile, shlex, vcf, pysam

parser = argparse.ArgumentParser(description='')
parser.add_argument ( '--input', dest='input', help='targ directory')
parser.add_argument ( '--input2', dest='control_input', help='control tag directory')
parser.add_argument ( '--mode', dest='mode', help='mode of operation')
parser.add_argument ( '-o', dest='output', help='output log file' )

def __main__():
    args = parser.parse_args()

    if args.input is not None:
        tag_dir=args.input

    if args.control_input is not None:
        cmd="findPeaks %s -style %s -o auto -i %s" % (tag_dir, args.mode, args.control_input)
    else:
        cmd="findPeaks %s -style %s -o auto" % (tag_dir, args.mode)

    print("cmd: %s" % cmd)

    tmp_dir = tempfile.mkdtemp()
    stdout = tempfile.NamedTemporaryFile( prefix="findpeaks-stdout-", dir=tmp_dir )
    stderr = tempfile.NamedTemporaryFile( prefix="findpeaks-stderr-", dir=tmp_dir )  

    proc = subprocess.Popen( args=cmd, stdout=stdout, stderr=stderr, shell=True, cwd=tmp_dir)
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
    if args.mode == "factor":
        shutil.copy('%s/peaks.txt' % tag_dir, args.output)
    elif args.mode == "histone":
        shutil.copy('%s/regions.txt' % tag_dir, args.output)
    elif args.mode == "groseq":
        shutil.copy('%s/transcripts.txt' % tag_dir, args.output)
    elif args.mode == "super":
        shutil.copy('%s/superEnhancers.txt' % tag_dir, args.output)    
    elif args.mode == "tss":
        shutil.copy('%s/tss.txt' % tag_dir, args.output)
    elif args.mode == "dnase":
        shutil.copy('%s/peaks.txt' % tag_dir, args.output)
    elif args.mode == "mC":
        shutil.copy('%s/regions.txt' % tag_dir, args.output)
    else:
        print("output file not found")
if __name__=="__main__":
	__main__()
