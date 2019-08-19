#!/usr/bin/python

#!/usr/bin/python

import optparse, os, shutil, sys, tempfile, glob, shlex, vcf, pysam, tarfile
from subprocess import *
import subprocess

CHUNK_SIZE = 2**20 #1mb

# global db definition
genome_db={}; project_db={}
genome_db['hg38']="postgres://bdds.globusgenomics.org/hg38"
genome_db['hg19']="postgres://bdds.globusgenomics.org/hg19"
genome_db['mm10']="postgres://bdds.globusgenomics.org/mm10"
genome_db['mm9']="postgres://bdds.globusgenomics.org/mm9"

project_db['wholebrain']="postgres://bdds.globusgenomics.org/wholebrain"
project_db['lymphoblast']="postgres://bdds.globusgenomics.org/lymphoblast"


def cleanup_before_exit( tmp_dir ):
    if tmp_dir and os.path.exists( tmp_dir ):
        shutil.rmtree( tmp_dir )

def __main__():
    parser = optparse.OptionParser()
    parser.add_option('','--genome', dest="genomedb", action='store', type="string", default=None, help='')
    parser.add_option('','--tissue', dest="projectdb", action='store', type="string", default=None, help='')
    parser.add_option('','--genelist', dest="genelist", action='store', type="string", default=None, help='')
    parser.add_option('','--tflist', dest="tflist", action='store', type="string", default=None, help='')
    parser.add_option('','--geneID', dest="geneID", action='store', type="string", default=None, help='')
    parser.add_option('','--upstreamsize', dest="upstreamsize", action='store', type="string", default=None, help='')
    parser.add_option('','--downstreamsize', dest="downstreamsize", action='store', type="string", default=None, help='')
    parser.add_option( '', '--output', dest='outputF', action='store', type="string", default=None, help='')
    parser.add_option( '', '--out-dir', dest='output_dir', action='store', type="string", default=None, help='If specified, the output directory for extra files.' )

    (options, args) = parser.parse_args()

    if not os.path.exists(options.output_dir): 
        os.mkdir(options.output_dir)

    cmd = "Rscript /opt/galaxy/tools/trena/getTFBSCountsInPromoters.R %s %s %s "  % (options.output_dir, genome_db[options.genomedb], project_db[options.projectdb])
    print cmd
#    print "hint cmd:" + hint_cmd
    stdout = tempfile.NamedTemporaryFile( prefix="trena-getcount-stdout-", dir=options.output_dir )
    stderr = tempfile.NamedTemporaryFile( prefix="trena-getcount-stderr-", dir=options.output_dir )
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
    shutil.copy('%s/promoter_tfbs_counts.RData' % options.output_dir, options.outputF)
    #cleanup_before_exit( options.output_dir )

if __name__=="__main__": __main__()

