#!/usr/bin/env python

"""
Runs Emirge on single-end or paired-end data.
usage: emirge_wrapper.py [options]

See below for options
"""

import optparse, os, shutil, subprocess, sys, tempfile
import time, glob, operator

def stop_err( msg ):
    sys.stderr.write( '%s\n' % msg )
    sys.exit()

def findNewestDir(directory):
    os.chdir(directory)
    dirs = {}
    for dir in glob.glob('*'):
        if os.path.isdir(dir):
            dirs[dir] = os.path.getctime(dir)

    lister = sorted(dirs.iteritems(), key=operator.itemgetter(1))
    return lister[-1][0]

def __main__():
    print "Start time:"
    print time.strftime('%d/%m/%Y %H:%M:%S\n', time.localtime(time.time()))

    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '-f', '--read1', dest='forward_read', help='forward fastq read' )
    parser.add_option( '-r', '--read2', dest='reverse_read', help='reverse fastq read' )
    parser.add_option( '', '--single', action='store_true', dest='single', help='single-ended fastq read' )
    parser.add_option( '', '--paired', action='store_true', dest='paired', help='paired-ended fastq read' )
    parser.add_option( '', '--insert-mean', dest='insert_mean', help='Insert mean' )
    parser.add_option( '', '--insert-stdev', dest='insert_stdev', help='insert size distribution standard deviation' )
    parser.add_option( '', '--amplicon', action='store_true', dest='amplicon_lib', help='Use amplicon lib script' )
    parser.add_option( '', '--max-read-length', dest='max_read_length', help='Max read length' )
    parser.add_option( '', '--db', dest='db', help='bowtie db' )
    parser.add_option( '', '--fasta', dest='fasta', help='fasta reference file' )
    parser.add_option( '', '--indexed', action='store_true', dest='indexed', help='Use indexed reference database files' )
    parser.add_option( '', '--phred33', action='store_true', dest='phred33', help='fastqsanger format' )
    parser.add_option( '', '--out', dest='outfile', help='output file' )
    parser.add_option( '', '--outdir', dest='outdir', help='output dir' )
    parser.add_option( '', '--num_threads', dest='nt', help='number of processors' )
    parser.add_option( '', '--iterations', dest='iterations', help='iterations' )
    parser.add_option( '', '--snp_fraction', dest='snp_fraction', help='snp_fraction' )
    parser.add_option( '', '--variant_fraction', dest='variant_fraction', help='variant_fraction' )
    parser.add_option( '', '--join_threshold', dest='join_threshold', help='join_threshold' )
    parser.add_option( '', '--read_depth', dest='read_depth', help='read_depth' )
    (options, args) = parser.parse_args()

    # check that output file is empty
    #if os.path.exists(options.outfile):
    #    os.remove(options.outfile)

    # make temp directory for placement of indices
    tmp_index_dir = tempfile.mkdtemp()
    ##tmp_dir = tempfile.mkdtemp()
    tmp_dir = options.outdir
    
    if not os.path.exists(options.outdir):
        os.mkdir(options.outdir)

    # index if necessary
    if options.fasta and not options.indexed:
        ref_file = tempfile.NamedTemporaryFile( dir=tmp_index_dir )
        ref_file_name = ref_file.name
        index_name = "%s_indexed" % ref_file_name
        ref_file.close()
        os.symlink( options.fasta, ref_file_name )
        cmd1 = 'bowtie-build %s %s' % ( ref_file_name, index_name )
        try:
            tmp = tempfile.NamedTemporaryFile( dir=tmp_index_dir ).name
            tmp_stderr = open( tmp, 'wb' )

	    print cmd1  

            proc = subprocess.Popen( args=cmd1, shell=True, cwd=tmp_index_dir, stderr=tmp_stderr.fileno() )
            returncode = proc.wait()
            tmp_stderr.close()
            # get stderr, allowing for case where it's very large
            tmp_stderr = open( tmp, 'rb' )
            stderr = ''
            buffsize = 1048576
            try:
                while True:
                    stderr += tmp_stderr.read( buffsize )
                    if not stderr or len( stderr ) % buffsize != 0:
                        break
            except OverflowError:
                pass
            tmp_stderr.close()
            if returncode != 0:
                raise Exception, stderr
        except Exception, e:
            # clean up temp dirs
            if os.path.exists( tmp_index_dir ):
                shutil.rmtree( tmp_index_dir )
            if os.path.exists( tmp_dir ):
                shutil.rmtree( tmp_dir )
            stop_err( 'Error indexing reference sequence. ' + str( e ) )
    else:
        ref_file_name = options.fasta
        index_name = options.db

    # set up command line for emirge 
    cmd = list()
    if options.amplicon_lib:
        cmd.append("emirge_amplicon.py")
    else:
        cmd.append("emirge.py")

    # add the output directory
    cmd.append(tmp_dir)

    # add the reads
    cmd.append("-1 %s" % options.forward_read)
    if options.paired:
        cmd.append("-2 %s" % options.reverse_read)
        cmd.append("-i %s" % options.insert_mean)
        cmd.append("-s %s" % options.insert_stdev)

    # fasta reference file
    cmd.append("-f %s" % ref_file_name)
    cmd.append("-b %s" % index_name)

    # read length
    cmd.append("-l %s" % options.max_read_length)

    # fastq format
    if options.phred33:
        cmd.append("--phred33")

    # advanced options
    if options.nt:
        cmd.append("-a %s" % options.nt)
    if options.iterations:
        cmd.append("-n %s" % options.iterations)
    if options.snp_fraction:
        cmd.append("-p %s" % options.snp_fraction)
    if options.variant_fraction:
        cmd.append("-v %s" % options.variant_fraction)
    if options.join_threshold:
        cmd.append("-j %s" % options.join_threshold)
    if options.read_depth:
        cmd.append("-c %s" % options.read_depth)

    cmd_line = " ".join(cmd)

    # perform cmd
    try:
        tmp = tempfile.NamedTemporaryFile( dir=tmp_dir ).name
        tmp_stderr = open( tmp, 'wb' )
        print cmd_line

        proc = subprocess.Popen( args=cmd_line, shell=True, cwd=tmp_dir, stderr=tmp_stderr.fileno() )
        returncode = proc.wait()
        tmp_stderr.close()
        # get stderr, allowing for case where it's very large
        tmp_stderr = open( tmp, 'rb' )
        stderr = ''
        buffsize = 1048576
        try:
            while True:
                stderr += tmp_stderr.read( buffsize )
                if not stderr or len( stderr ) % buffsize != 0:
                    break
        except OverflowError:
            pass
        tmp_stderr.close()
        if returncode != 0:
            raise Exception, stderr

        ### get the fasta file from the last iteration
        iteration_dir = findNewestDir(tmp_dir)
        get_fasta_cmd = "emirge_rename_fasta.py %s > %s" % (iteration_dir, options.outfile)

        tmp_stderr = open( tmp, 'wb' )
        proc = subprocess.Popen( args=get_fasta_cmd, shell=True, cwd=tmp_dir, stderr=tmp_stderr.fileno() )
        returncode = proc.wait()
        tmp_stderr.close()
        # get stderr, allowing for case where it's very large
        tmp_stderr = open( tmp, 'rb' )
        stderr = ''
        buffsize = 1048576
        try:
            while True:
                stderr += tmp_stderr.read( buffsize )
                if not stderr or len( stderr ) % buffsize != 0:
                    break
        except OverflowError:
            pass
        tmp_stderr.close()
        if returncode != 0:
            raise Exception, stderr

    finally:
        # clean up temp dir
        if os.path.exists( tmp_index_dir ):
            shutil.rmtree( tmp_index_dir )
        #if os.path.exists( tmp_dir ):
        #    shutil.rmtree( tmp_dir )

    print "\n\nEnd time:"
    print time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))

if __name__=="__main__": __main__()
