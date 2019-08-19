#!/usr/bin/env python

"""
A wrapper script for running freebayes (take directory as input).
"""

import sys, optparse, os, tempfile, subprocess, shutil
from binascii import unhexlify
from string import Template
import time  

CHUNK_SIZE = 2**20 #1mb

def cleanup_before_exit( tmp_dir ):
    if tmp_dir and os.path.exists( tmp_dir ):
        shutil.rmtree( tmp_dir )

def open_file_from_option( filename, mode = 'rb' ):
    if filename:
        return open( filename, mode = mode )
    return None

def index_bam_files( bam_filenames, tmp_dir ):
    for bam_filename in bam_filenames:
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

def __main__():
    #liubo added, print current time for calculating execution time of GATK
    print "Start time:"
    print time.strftime('%d/%m/%Y %H:%M:%S\n', time.localtime(time.time()))

    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '-p', '--pass_through1', dest='pass_through_options1', action='append', type="string", help='These options are passed through directly, without any modification.' )
    parser.add_option( '-q', '--pass_through2', dest='pass_through_options2', action='append', type="string", help='These options are passed through directly, without any modification.' )
    parser.add_option( '', '--input_dir', dest='input_dir', action='store', type="string", help='Input directory path of BAM files' )
    parser.add_option( '', '--input_dir_file', dest='input_dir_file', action='store', type="string", help='Input directory File containing path of BAM files' )

    parser.add_option( '', '--stdout', dest='stdout', action='store', type="string", default=None, help='If specified, the output of stdout will be written to this file.' )
    parser.add_option( '', '--stderr', dest='stderr', action='store', type="string", default=None, help='If specified, the output of stderr will be written to this file.' )   
    (options, args) = parser.parse_args()

    tmp_dir = tempfile.mkdtemp( prefix='tmp-freebayes-' )
    #set up stdout and stderr output options
    stdout = open_file_from_option( options.stdout, mode = 'wb' )
    stderr = open_file_from_option( options.stderr, mode = 'wb' )
    #if no stderr file is specified, we'll use our own
    if stderr is None:
        stderr = tempfile.NamedTemporaryFile( prefix="freebayes-stderr-", dir=tmp_dir )

#run the command for setting up input files
    if options.pass_through_options1:
        cmd1 = ' '.join( options.pass_through_options1 )
        proc = subprocess.Popen( args=cmd1, stdout=stdout, stderr=stderr, shell=True, cwd=tmp_dir )


#traverse all the bam files in the input directory and format them to "-I file1.bam -I file2.bam..."
#options.input_dir is the input directory's full path
    if options.input_dir_file:
	infile = open(options.input_dir_file, 'r')
        inputDirLine = infile.readline()
        inputDirectory = inputDirLine.rstrip('\n\r ')
    else:
        inputDirectory = options.input_dir

    if inputDirectory:  
        flist = [x for x in os.listdir(inputDirectory) if not x.startswith('.') and x.endswith('.bam')] 
	inputBAM = ""
        if len(flist) > 0:
	    for root, dirs, files in os.walk(inputDirectory): 
		for name in files: 
		    if not name.startswith('.') and name.endswith('.bam'):
                        inputBAM += " --bam " + os.path.join(root, name)


    if options.pass_through_options2:
        pass_through_options2 = ' '.join( options.pass_through_options2 )

    cmd2 = "freebayes " + inputBAM + pass_through_options2

    print cmd2
 
    proc = subprocess.Popen( args=cmd2, stdout=stdout, stderr=stderr, shell=True, cwd=tmp_dir )
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
    
    cleanup_before_exit( tmp_dir )


    #print current time for calculating execution time of GATK
    print "\n\nEnd time:"
    print time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))

if __name__=="__main__": __main__()
