#!/usr/bin/python

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
    parser.add_option('','--counts', dest="counts", action='store', type="string", default=None, help='')
    parser.add_option('','--expr', dest="expr", action='store', type="string", default=None, help='')
    parser.add_option('','--method', dest="method", action='store', type="string", default=None, help='')
    parser.add_option('','--alpha', dest="alpha", action='store', type="string", default=None, help='')
    parser.add_option('','--regulatormethod', dest="regulatormethod", action='store', type="string", default=None, help='')
    parser.add_option('','--tfbsthreshold', dest="tfbsthreshold", action='store', type="string", default=None, help='')
    parser.add_option( '', '--output', dest='outputF', action='store', type="string", default=None, help='')
    parser.add_option( '', '--out-dir', dest='output_dir', action='store', type="string", default=None, help='If specified, the output directory for extra files.' )

    (options, args) = parser.parse_args()

    if not os.path.exists(options.output_dir): 
        os.mkdir(options.output_dir)

    cmd = "Rscript /opt/galaxy/tools/trena/makeTRN-proximalonly-demo.R %s %s %s %s"  % (options.output_dir, options.counts, options.expr, options.method)
    print cmd
#    print "hint cmd:" + hint_cmd
    stdout = tempfile.NamedTemporaryFile( prefix="trena-maketrn-stdout-", dir=options.output_dir )
    stderr = tempfile.NamedTemporaryFile( prefix="trena-maketrn-stderr-", dir=options.output_dir )
    proc = subprocess.Popen( args=cmd, stdout=stdout, stderr=stderr, shell=True, cwd=options.output_dir )
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
    shutil.copy('%s/trn.RData' % options.output_dir, options.outputF)
    #cleanup_before_exit( options.output_dir )

if __name__=="__main__": __main__()

