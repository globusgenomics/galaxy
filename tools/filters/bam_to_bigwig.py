#!/usr/bin/env python
"""
Convert BAM files to BigWig file format
"""
import os
import sys
import subprocess
import tempfile, shutil
from optparse import OptionParser
from contextlib import contextmanager, closing

import pysam

CHUNK_SIZE = 2**20 #1mb

@contextmanager
def indexed_bam(bam_file, config):
    indexFile="%s.bai" % bam_file
    if not os.path.exists(indexFile):
        pysam.index(bam_file)
    sam_reader = pysam.Samfile(bam_file, "rb")
    yield sam_reader
    sam_reader.close()


def get_sizes(bam_file, config):
    with indexed_bam(bam_file, config) as work_bam:
        sizes = zip(work_bam.references, work_bam.lengths)
    return sizes

def cleanup_before_exit( tmp_dir ):
    if tmp_dir and os.path.exists( tmp_dir ):
        shutil.rmtree( tmp_dir, ignore_errors=True )

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("", "--input", dest="inputfile")
    parser.add_option("-o", "--output", dest="outfile")
    parser.add_option( '', '--out-dir', dest='output_dir', action='store', type="string", default=None, help='If specified, the output directory for extra files.' )
    parser.add_option( '-p', '--pass_through', dest='pass_through_options', action='append', type="string", help='These options are passed through directly to contra, without any modification.' )
    (options, args) = parser.parse_args()

    config = {"program": {"ucsc_bedGraphToBigWig": ["bedGraphToBigWig"],
                          "bedtools_genomeCoverageBed":
                          ["bedtools", "genomecov"]}}

    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)

    output_dir = "%s/output" % options.output_dir
    inputDirectory = "%s/input" % options.output_dir

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
        os.mkdir(inputDirectory)

    tmp_dir = tempfile.mkdtemp( dir=options.output_dir , prefix='tmp-TOOL-' )
    samfile = pysam.AlignmentFile(options.inputfile, 'rb')
    if "RG" in samfile.header:
        sn = samfile.header['RG'][0]['SM']
    else:
        sn = "tmpname"

    # create symlink for input file
    os.symlink(options.inputfile, "%s/%s.bam" % (inputDirectory, sn))
    inputLinkedFile=("%s/%s.bam" % (inputDirectory, sn))

    # pass through options
    ptc = ""
    if options.pass_through_options:
        ptc = ' '.join( options.pass_through_options )

    sizes = get_sizes(inputLinkedFile, config)
    print "Have %i references" % len(sizes)
    if not sizes:
        sys.stderr.write("Problem reading BAM header.\n")
        sys.exit(1)

    # Use a temp file to avoid any possiblity of not having write permission
    print "Calculating coverage..."
    outbase="%s/%s" % (tmp_dir,sn)
    bedgraph_tmp="%s.bed" % outbase
    bedgraph_file="%s_graph.bed" % outbase
    bam2graph_cmd = "bedtools genomecov -ibam  %s -bg %s > %s; sort -k1,1 -k2,2n %s > %s" % (inputLinkedFile, ptc, bedgraph_tmp, bedgraph_tmp, bedgraph_file)

    try:
        stdout = tempfile.NamedTemporaryFile( prefix="bam2bigwig-stdout-", dir=tmp_dir )
        stderr = tempfile.NamedTemporaryFile( prefix="bam2bigwig-stderr-", dir=tmp_dir )
        proc = subprocess.Popen( args=bam2graph_cmd, stdout=stdout, stderr=stderr, shell=True, cwd=tmp_dir )
        proc.communicate()
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
        stdout.close()
        stderr.close()
    except Exception, e:
        sys.stderr.write("problem doing : %s\n" %(bam2graph_cmd))
        sys.stderr.write( '%s\n\n' % str(e) )

    size_file = "%s-sizes.txt" % outbase
    with open(size_file, "w") as out_handle:
        for chrom, size in sizes:
            out_handle.write("%s\t%s\n" % (chrom, size))
    print("Converting %i MB graph file to bigwig..." % (os.path.getsize(bedgraph_file) // (1024 * 1024)))

    cl = config["program"]["ucsc_bedGraphToBigWig"] + \
        [bedgraph_file, size_file, options.outfile]
    subprocess.check_call(cl)
    print options.output_dir
    cleanup_before_exit(tmp_dir)   
