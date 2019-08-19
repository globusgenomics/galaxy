#!/usr/bin/env python

"""
A wrapper script for running the fusioncatcher.
"""

import sys, optparse, os, tempfile, subprocess, shutil
from binascii import unhexlify
from string import Template
import time

DEFAULT_ISAAC_PREFIX = "fusioncatcher_file"
CHUNK_SIZE = 2**20 #1mb


def cleanup_before_exit( tmp_dir ):
    if tmp_dir and os.path.exists( tmp_dir ):
        shutil.rmtree( tmp_dir )


def filename_from_galaxy( galaxy_filename, galaxy_ext, prefix, normal_target_dir = None, target_dir = None, lane_id = None, direction = None):
    if direction == "forward":
        read_direction = "1"
    elif direction == "reverse":
        read_direction = "2"
 
    if galaxy_ext == "fastq":
        extension = "fastq"
    elif galaxy_ext == "sra":
        extension = "sra"

    if "normal" in prefix:
        fusioncatcher_filename = os.path.join( normal_target_dir, "lane%s_read%s.%s" % ( lane_id, read_direction, extension ) )
    else:
        fusioncatcher_filename = os.path.join( target_dir, "lane%s_read%s.%s" % ( lane_id, read_direction, extension ) )
    #print "ln -s %s %s" % ( galaxy_filename, isaac_filename )
    os.symlink( galaxy_filename, fusioncatcher_filename )
    return fusioncatcher_filename

def open_file_from_option( filename, mode = 'rb' ):
    if filename:
        return open( filename, mode = mode )
    return None

def __main__():
    print "Start time:"
    print time.strftime('%d/%m/%Y %H:%M:%S\n', time.localtime(time.time()))

    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '-p', '--pass_through', dest='pass_through_options', action='append', type="string", help='These options are passed through directly to ISAAC, without any modification.' )
    parser.add_option( '-d', '--dataset', dest='datasets', action='append', type="string", nargs=6, help='"-argument" "original_filename" "galaxy_filetype" "name_prefix" "lane_id" "forward_or_reverse"' )
    parser.add_option( '', '--extra-output-dir', dest='output_dir', action='store', type="string", default=None, help='output_directory' )
    parser.add_option( '-o', '--fusion_genes_output', dest='output_fusion_genes_list', action='store', type="string", default=None, help='output_bam_path' )
    parser.add_option( '-R', '--reference', dest='reference', action='store', type="string", default=None, help='Reference file path' )
    (options, args) = parser.parse_args()

    # make tmpdir and output dir   
    os.mkdir(options.output_dir)
    tmp_dir = tempfile.mkdtemp( prefix='tmp--fusioncatcher', dir=options.output_dir )
    input_link_dir = "%s/inputs" % tmp_dir
    normal_link_dir = "%s/normal_inputs" % tmp_dir
    fusioncatcher_tmp_dir = "%s/tmp" % tmp_dir
    if not os.path.exists(input_link_dir):
        os.mkdir( input_link_dir )
    if not os.path.exists(normal_link_dir):
        os.mkdir( normal_link_dir )
    if not os.path.exists(fusioncatcher_tmp_dir):
        os.mkdir( fusioncatcher_tmp_dir )

    normal_sample_flag = ""
    if options.datasets:
        for ( dataset_arg, filename, galaxy_ext, prefix, lane_id, direction ) in options.datasets:
            filename_from_galaxy( filename, galaxy_ext, prefix, normal_target_dir = normal_link_dir, target_dir = input_link_dir, lane_id = lane_id, direction = direction)
            if "normal" in prefix:
                normal_sample_flag = "--I %s" % normal_link_dir

    # prepare command
    cmd = 'fusioncatcher -i %s %s -d %s -o %s %s' % ( input_link_dir, normal_sample_flag, options.reference, options.output_dir, ' '.join( options.pass_through_options ) )


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
 
    proc = subprocess.Popen( args=cmd, stdout=open( stdout_name, 'wb' ), stderr=open( stderr_name, 'wb' ), shell=True, cwd=fusioncatcher_tmp_dir )

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
    fusion_candidates_file_out = "%s/final-list_candidate-fusion-genes.txt" % (options.output_dir)
    shutil.move(fusion_candidates_file_out, options.output_fusion_genes_list)

    cleanup_before_exit( tmp_dir )

    print "\n\nEnd time:"
    print time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))

if __name__=="__main__": __main__()
