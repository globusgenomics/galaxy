#!/usr/bin/python

#!/usr/bin/python

import optparse, os, shutil, sys, tempfile, glob, shlex, tarfile, pysam
from subprocess import *
import subprocess

CHUNK_SIZE = 2**20 #1mb

def cleanup_before_exit( tmp_dir ):
    if tmp_dir and os.path.exists( tmp_dir ):
        shutil.rmtree( tmp_dir )

def __main__():
    parser = optparse.OptionParser()
    parser.add_option( '', '--reference_name', dest='refName', action='store', type="string", help='' )
    parser.add_option( '', '--reference', dest='ref', action='store', type="string", help='' )
    parser.add_option( '', '--bam', dest='inputbam', action='store', type="string", help='' )
    #parser.add_option( '', '--summary', dest='summaryF', action='store', type="string", default=None, help='')
    parser.add_option( '', '--log', dest='logF', action='store', type="string", default=None, help='')
    parser.add_option( '', '--output', dest='outputF', action='store', type="string", default=None, help='')
    parser.add_option( '', '--output-dir', dest='output_dir', action='store', type="string", default=None, help='If specified, the output directory for extra files.' )
    parser.add_option( '', '--checkChr', dest='checkchr', action='store', type="string", default=None, help='If specified, the output directory for extra files.' )
    (options, args) = parser.parse_args()

    # worker node tmp directory where all will be processed
    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)
        #subprocess.call(['chmod', '-R', '0777', options.output_dir])
        #os.chmod(options.output_dir, 0o777)

    # create input files from the input_dir
    shutil.copy(options.inputbam, "%s/input1.bam" % options.output_dir)
    pysam.index("%s/input1.bam" % options.output_dir)
    subprocess.call(['chmod', '-R', '0777', options.output_dir])
    if options.checkchr == "false":
        chrList="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"
    else:
        chrList="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"

    cnvnator_cmd = "docker run -v %s:%s -v /tmp:/tmp mustxyk/ubuntu-cnvnator /root/run_cnvnator.sh %s %s %s/cnvnator.root \"%s\" 500 %s/cnvnator.log %s/cnvnator.out %s/input1.bam" % (options.output_dir, options.output_dir, options.refName, options.ref, options.output_dir, chrList, options.output_dir, options.output_dir, options.output_dir)

    print cnvnator_cmd

    stdout = tempfile.NamedTemporaryFile( prefix="cnvnator-stdout-", dir=options.output_dir )
    stderr = tempfile.NamedTemporaryFile( prefix="cnvnator-stderr-", dir=options.output_dir )
    proc = subprocess.Popen( args=cnvnator_cmd, stdout=stdout, stderr=stderr, shell=True, cwd=options.output_dir )
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
    #shutil.copy('%s/cnvnator.root' % options.output_dir, options.summaryF)
    shutil.copy('%s/cnvnator.log' % options.output_dir, options.logF)
    shutil.copy('%s/cnvnator.out' % options.output_dir, options.outputF)
    cleanup_before_exit( options.output_dir )

if __name__=="__main__": __main__()
