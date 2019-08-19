#!/usr/bin/python

CHUNK_SIZE = 2**20 #1mb

import argparse, os, shutil, subprocess, sys, tempfile, shlex, vcf, pysam

parser = argparse.ArgumentParser(description='')
parser.add_argument ( '--input', dest='input', help='the bed file')
parser.add_argument ( '-o', dest='output', help='output log file' )

def execute( cmd, output="" ):
    tmp_dir = tempfile.mkdtemp()
    try:
            err = open(tmp_dir+"/errorLog", 'a')
            if output != "":
                    out = open(output, 'w')
            else:
                    out = subprocess.PIPE
            process = subprocess.Popen( args=shlex.split(cmd), stdout=out, stderr=err )
            process.wait()
            err.close()
            if out != subprocess.PIPE:
                    out.close()
    except Exception, e:
            sys.stderr.write("problem doing : %s\n" %(cmd))
            sys.stderr.write( '%s\n\n' % str(e) )

def __main__():
    args = parser.parse_args()

    if args.input is not None:
        cmd="makeUCSCfile %s -o bedgraph" % args.input
#        cmd="makeUCSCfile %s -o %s" % (args.input, args.output) 
        
    print("cmd: %s" % cmd)

    tmp_dir = tempfile.mkdtemp()
    stdout = tempfile.NamedTemporaryFile( prefix="makeucscfile-stdout-", dir=tmp_dir )
    stderr = tempfile.NamedTemporaryFile( prefix="makeucscfile-stderr-", dir=tmp_dir )  

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

    shutil.copy('%s/bedgraph.gz' % tmp_dir, args.output)

if __name__=="__main__":
	__main__()
