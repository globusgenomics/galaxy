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


#isaac_filename = isaac_filename_from_galaxy( filename, galaxy_ext, target_dir = tmp_dir, lane_id = lane_id, direction = direction)
def isaac_filename_from_galaxy( galaxy_filename, galaxy_ext, target_dir = None, lane_id = None, direction = None):
    if direction == "forward":
        read_direction = "1"
    elif direction == "reverse":
        read_direction = "2"
 
    if galaxy_ext == "fastq":
        extension = "fastq"
    elif galaxy_ext == "fastq_gz":
        extension = "fastq.gz"
    elif galaxy_ext == "bam":
        extension = "bam"
    elif galaxy_ext == "bcl":
        extension = "bcl"

    isaac_filename = os.path.join( target_dir, "lane%s_read%s.%s" % ( lane_id, read_direction, extension ) )
    #print "ln -s %s %s" % ( galaxy_filename, isaac_filename )
    os.symlink( galaxy_filename, isaac_filename )
    return isaac_filename

def open_file_from_option( filename, mode = 'rb' ):
    if filename:
        return open( filename, mode = mode )
    return None

def prepare_flowcell_link_dir( input_tmpdir, filename, lane_id):
    # create the necessar links and directory structure for a flowcell
    source_path = os.path.dirname(filename)
    os.symlink(filename, "%s/%s" % (input_tmpdir, "RunInfo.xml"))
    os.mkdir("%s/%s" % (input_tmpdir, "Data"))
    os.mkdir("%s/%s" % (input_tmpdir, "Data/Intensities"))
    if lane_id == "all":
       os.symlink("%s/Data/Intensities/BaseCalls" % source_path)
    else:
       os.mkdir("%s/%s" % (input_tmpdir, "Data/Intensities/BaseCalls"))
       os.symlink("%s/Data/Intensities/BaseCalls/%s" % (source_path, lane_id), "%s/Data/Intensities/BaseCalls/%s" % (input_tmpdir, lane_id))


def __main__():
    print "Start time:"
    print time.strftime('%d/%m/%Y %H:%M:%S\n', time.localtime(time.time()))

    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '-p', '--pass_through', dest='pass_through_options', action='append', type="string", help='These options are passed through directly to ISAAC, without any modification.' )
    parser.add_option( '-d', '--dataset', dest='datasets', action='append', type="string", nargs=6, help='"-argument" "original_filename" "galaxy_filetype" "name_prefix" "lane_id" "forward_or_reverse"' )
    parser.add_option( '', '--extra-output-dir', dest='output_dir', action='store', type="string", default=None, help='output_directory' )
    parser.add_option( '', '--outbam', dest='output_bam', action='store', type="string", default=None, help='output_bam_path' )
    parser.add_option( '-R', '--reference', dest='reference_xml', action='store', type="string", default=None, help='Reference XML file path' )
    parser.add_option( '', '--dbname', dest='reference_name', action='store', type="string", default=None, help='Reference name' )
    parser.add_option( '', '--format', dest='input_format', action='store', type="string", default=None, help='input format' )
    (options, args) = parser.parse_args()

    # make tmpdir and output dir
    output_dir_tmp = tempfile.mkdtemp( prefix='tmp-isaac-output')
    #output_dir_tmp = options.output_dir
    if not os.path.isdir(output_dir_tmp):
        os.mkdir(output_dir_tmp)
    tmp_dir = tempfile.mkdtemp( prefix='tmp--isaac', dir=output_dir_tmp )
    input_link_dir = "%s/inputs" % tmp_dir
    isaac_tmp_dir = "%s/tmp" % tmp_dir
    os.mkdir( input_link_dir )
    os.mkdir( isaac_tmp_dir )

    if options.input_format != "bcl":
        input_filenames = []
        input_format = None
        if options.datasets:
            if options.input_format != "input_dir":
                for ( dataset_arg, filename, galaxy_ext, prefix, lane_id, direction ) in options.datasets:
                    isaac_filename = isaac_filename_from_galaxy( filename, galaxy_ext, target_dir = input_link_dir, lane_id = lane_id, direction = direction)
                    input_filenames.append( isaac_filename )
                    input_format = galaxy_ext
                    if input_format == 'fastq_gz':
                        input_format = 'fastq-gz'
            else:
                # input_dir from BCL2fastq tool
                for ( dataset_arg, filename, galaxy_ext, prefix, sample_name, direction ) in options.datasets:
                    input_link_dir = "%s/%s" % (filename, sample_name)
                    input_format = 'fastq-gz'
        # prepare command
        cmd = 'isaac-align -r %s -n %s -o %s -t %s -b %s --base-calls-format %s %s' % ( options.reference_xml, options.reference_name, output_dir_tmp, isaac_tmp_dir, input_link_dir, input_format, ' '.join( options.pass_through_options ) )
    else:
        for ( dataset_arg, filename, galaxy_ext, prefix, lane_id, direction ) in options.datasets:
        
            # prepare directory in temporary space. Create necessary directory structure and links to run file
            input_tmpdir = tempfile.mkdtemp( prefix='tmp-isaac-input-bcl' )
            prepare_flowcell_link_dir(input_tmpdir, filename, lane_id)
            xmlfile = "%s/RunInfo.xml" % input_tmpdir

        # assume bcl-gz format for all bcl files
        input_format = "bcl-gz"
        cmd = 'isaac-align -r %s -n %s -o %s -t %s -b %s --base-calls-format %s %s' % ( options.reference_xml, options.reference_name, output_dir_tmp, isaac_tmp_dir, xmlfile, input_format, ' '.join( options.pass_through_options ) )


    #set up stdout and stderr output options
    stdout_name = tempfile.NamedTemporaryFile( prefix = "stdout" ).name
    stderr_name = tempfile.NamedTemporaryFile( prefix = "stderr" ).name
    #stdout = open_file_from_option( stdout_name, mode = 'wb' )
    #stderr = open_file_from_option( stderr_name, mode = 'wb' )
    #if no stderr file is specified, we'll use our own
    #if stderr is None:
    #    stderr = tempfile.NamedTemporaryFile( prefix="isaac-stderr-", dir=tmp_dir )
   
    print cmd
    buffsize = 1048576
 
    proc = subprocess.Popen( args=cmd, stdout=open( stdout_name, 'wb' ), stderr=open( stderr_name, 'wb' ), shell=True, cwd=isaac_tmp_dir )

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

    # place output BAM to output file
    bam_file_out = "%s/Projects/%s/%s/sorted.bam" % (output_dir_tmp, options.reference_name, options.reference_name)
    shutil.move(bam_file_out, options.output_bam)

    cleanup_before_exit( tmp_dir )

    print "\n\nEnd time:"
    print time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))

if __name__=="__main__": __main__()
