#!/usr/bin/python

CHUNK_SIZE = 2**20 #1mb

import argparse, os, shutil, subprocess, sys, tempfile, shlex, vcf, pysam

parser = argparse.ArgumentParser(description='')
parser.add_argument ( '--input', dest='input_bed', help='input bed file')
parser.add_argument ( '--genome', dest='genome', help='genome')
parser.add_argument ( '--size', dest='size', help='genome')
#parser.add_argument ( '--output-log', dest='output_log', help='output log file' )
#parser.add_argument ( '--output-dir', dest='output_dir', help='output directory')
parser.add_argument ( '-p', action='append', dest='pass_through_options', help='These options are passed through directly to contra, without any modification')
parser.add_argument ( '--known-results-html', dest='known_results_html', help='known results html' )
parser.add_argument ( '--homer-results-html', dest='homer_results_html', help='homer results html' )
parser.add_argument ( '--out-concat-motifs', dest='out_concat_motifs', help='out concat motifs' )
parser.add_argument ( '--known-results-tabular', dest='known_results_tabular', help='known results tabular' )

def __main__():
    args = parser.parse_args()
    print("ptc: %s" % args.pass_through_options)
#    if args.output_dir is not None:
#        output_dir = args.output_dir
#        if not os.path.exists(output_dir):
#            os.mkdir(output_dir)
#    else:
#        output_dir = tempfile.mkdtemp();
    output_dir = tempfile.mkdtemp(prefix="homer-");
    if args.pass_through_options is not None:
        ptc = ' '.join(args.pass_through_options )
    print("ptc: %s" % ptc)
    ncpu=4
    homer="/mnt/galaxyTools/tools/homer/4.9"
    homer_bin="/mnt/galaxyTools/tools/homer/4.9/bin"
    env="export PATH=%s:%s:$PATH;" % (homer, homer_bin)
    cmd="%s perl %s/findMotifsGenome.pl %s %s %s -size %s %s -p %d" % (env, homer_bin, args.input_bed, args.genome, output_dir, args.size, ptc, ncpu)

    print("cmd: %s" % cmd)

    stdout = tempfile.NamedTemporaryFile( prefix="findmotifsgenome-stdout-", dir=output_dir )
    stderr = tempfile.NamedTemporaryFile( prefix="findmotifsgenome-stderr-", dir=output_dir )  

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

    stderr.close()
    stdout.close() 

    shutil.copy('%s/knownResults.html' % output_dir, args.known_results_html)
    shutil.copy('%s/homerResults.html' % output_dir, args.homer_results_html)    
    shutil.copy('%s/homerMotifs.all.motifs' % output_dir, args.out_concat_motifs)
    shutil.copy('%s/knownResults.txt' % output_dir, args.known_results_tabular)

if __name__=="__main__":
	__main__()
