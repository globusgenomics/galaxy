#!/usr/bin/env python

"""
A wrapper script for running the isaac commands.
"""

import sys, optparse, os, tempfile, subprocess, shutil
from binascii import unhexlify
from string import Template
import time

DEFAULT_ISAAC_PREFIX = "isaac_file"
CHUNK_SIZE = 2**20 #1mb


def cleanup_before_exit( tmp_dir ):
    if tmp_dir and os.path.exists( tmp_dir ):
        shutil.rmtree( tmp_dir )

def index_bam_files( bam_filename, tmp_dir ):
    bam_index_filename = "%s.bai" % bam_filename
    if not os.path.exists( bam_index_filename ):
        #need to index this bam file
        stderr_name = tempfile.NamedTemporaryFile( prefix = "bam_index_stderr" ).name
        command = 'samtools index %s %s' % ( bam_filename, bam_index_filename )
        proc = subprocess.Popen( args=command, shell=True, stderr=open( stderr_name, 'wb' ) )
        return_code = proc.wait()
        if return_code:
            for line in open( stderr_name ):
                print >> sys.stderr, line
            os.unlink( stderr_name ) #clean up
            cleanup_before_exit( tmp_dir )
            raise Exception( "Error indexing BAM file" )
        os.unlink( stderr_name ) #clean up


def open_file_from_option( filename, mode = 'rb' ):
    if filename:
        return open( filename, mode = mode )
    return None

def __main__():
    print "Start time:"
    print time.strftime('%d/%m/%Y %H:%M:%S\n', time.localtime(time.time()))

    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '-i', '--input', dest='input_bam', action='store', type="string", default=None, help='input BAM file' )
    parser.add_option( '-x', '--index', dest='input_bam_index', action='store', type="string", default=None, help='input BAM index file' )
    parser.add_option( '-R', '--db', dest='reference', action='store', type="string", default=None, help='Reference file path' )
    parser.add_option( '-t', '--type', dest='analysis_type', action='store', type="string", default=None, help='analysis type wes or wgs' )
    parser.add_option( '-o', '--output-gvcf', dest='output_gvcf', action='store', type="string", default=None, help='output path' )
    parser.add_option( '-c', '--ncpu', dest='ncpu', action='store', type="string", default=None, help='CPU qty' )
    ##parser.add_option( '-v', '--output-vcf', dest='output_vcf', action='store', type="string", default=None, help='output path' )
    (options, args) = parser.parse_args()

    # make tmpdir and output dir
    output_dir_tmp = tempfile.mkdtemp( prefix='tmp-isaac-output')
    if not os.path.isdir(output_dir_tmp):
        os.mkdir(output_dir_tmp)
    #tmp_dir = tempfile.mkdtemp( prefix='tmp-variant-isaac', dir=output_dir_tmp )
    tmp_dir = "%s/work_dir" % output_dir_tmp
    input_tmp_dir = tempfile.mkdtemp( prefix='tmp-inputs', dir=output_dir_tmp )
    input_bam_link = "%s/input.bam" % input_tmp_dir
    input_bam_index_link = "%s/input.bam.bai" % input_tmp_dir
    os.symlink(options.input_bam, input_bam_link)
    if options.input_bam_index and os.path.exists(options.input_bam_index):
        os.symlink(options.input_bam_index, input_bam_index_link)
    else:
        index_bam_files(input_bam_link, input_tmp_dir)

    # get the config file path
    config_path = None
    if options.analysis_type == "wgs":
        config_path = os.environ['ISAAC_VARIANT_CONFIG_WHOLE_PATH']
    elif options.analysis_type == "wes":
        config_path = os.environ['ISAAC_VARIANT_CONFIG_EXOME_PATH']
    else:
        sys.exit("Configuration file NOT_AVAILABLE")

    # copy the config file to the working directory
    job_config_file = "%s/config.txt" % output_dir_tmp
    shutil.copyfile(config_path,job_config_file)
    os.chdir(output_dir_tmp)

    cmd = "configureWorkflow.pl --bam=%s --ref=%s --config=%s --output-dir=%s" % (input_bam_link, options.reference, job_config_file, tmp_dir)
    cmd1 = "make -j 30 -C %s" % tmp_dir
    #cmd1 = "make -j %s -C %s" % (options.ncpu, tmp_dir)


    #set up stdout and stderr output options
    stdout_name = tempfile.NamedTemporaryFile( prefix = "stdout" ).name
    stderr_name = tempfile.NamedTemporaryFile( prefix = "stderr" ).name
   
    print cmd
    buffsize = 1048576
 
    proc = subprocess.Popen( args=cmd, stdout=open( stdout_name, 'wb' ), stderr=open( stderr_name, 'wb' ), shell=True, cwd=output_dir_tmp )

    return_code = proc.wait()
    if return_code:
        tmp_stderr = open( stderr_name, 'rb' )
        stderr = ''
        try:
            while True:
                stderr += tmp_stderr.read( buffsize )
                if not stderr or len( stderr ) % buffsize != 0:
                    break
        except OverflowError:
            pass
        tmp_stderr.close()
        if return_code != 0:
            raise Exception, stderr

    os.chdir(tmp_dir)
    #set up stdout and stderr output options for CMD1 run
    stdout_name = tempfile.NamedTemporaryFile( prefix = "stdout" ).name
    stderr_name = tempfile.NamedTemporaryFile( prefix = "stderr" ).name
    
    print cmd1
    buffsize = 1048576

    proc = subprocess.Popen( args=cmd1, stdout=open( stdout_name, 'wb' ), stderr=open( stderr_name, 'wb' ), shell=True, cwd=tmp_dir )

    return_code = proc.wait()
    if return_code:
        tmp_stderr = open( stderr_name, 'rb' )
        stderr = ''
        try:
            while True:
                stderr += tmp_stderr.read( buffsize )
                if not stderr or len( stderr ) % buffsize != 0:
                    break
        except OverflowError:
            pass
        tmp_stderr.close()
        if return_code != 0:
            raise Exception, stderr

    # place output GVCF to output file
    gvcf_file_out = "%s/results/%s.genome.vcf.gz" % (tmp_dir, "input")
    shutil.move(gvcf_file_out, options.output_gvcf)


##    # Gunzip output file
#    import gzip
#    gvcf_file_out = "%s/results/%s.genome.vcf.gz" % (tmp_dir, "input")
#    #vcf_file_out = "%s/results/%s.genome.vcf" % (tmp_dir, "input")
#    inF = gzip.GzipFile(gvcf_file_out, 'rb')
#    s = inF.read()
#    inF.close()
#
#    outF = file(options.output_gvcf, 'wb')
#    outF.write(s)
#    outF.close()
#
#    # Extract the variant regions from the GVCF (create the VCF)
#    cmd2 = "cat %s | extract_variants > %s" % (options.output_gvcf, options.output_vcf)
#    #set up stdout and stderr output options
#    stdout_name = tempfile.NamedTemporaryFile( prefix = "stdout" ).name
#    stderr_name = tempfile.NamedTemporaryFile( prefix = "stderr" ).name
#
#    print cmd2
#    buffsize = 1048576
#
#    proc = subprocess.Popen( args=cmd2, stdout=open( stdout_name, 'wb' ), stderr=open( stderr_name, 'wb' ), shell=True, cwd=output_dir_tmp )
#
#    return_code = proc.wait()
#    if return_code:
#        tmp_stderr = open( stderr_name, 'rb' )
#        stderr = ''
#        try:
#            while True:
#                stderr += tmp_stderr.read( buffsize )
#                if not stderr or len( stderr ) % buffsize != 0:
#                    break
#        except OverflowError:
#            pass
#        tmp_stderr.close()
#        if return_code != 0:
#            raise Exception, stderr

    cleanup_before_exit( tmp_dir )




    print "\n\nEnd time:"
    print time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))

if __name__=="__main__": __main__()
