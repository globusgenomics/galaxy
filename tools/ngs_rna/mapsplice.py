#!/usr/bin/env python

"""  
python mapsplice_segments.py -Q fq -o 1M_36bp_output_path -c chr20_sequence_index_path
 -u reads_path/1M_36bp_fastq.txt -B chr20_sequence_index _path/index -L 18 2>36bp_time.log
"""

import sys, optparse, os, tempfile, subprocess, shutil, re


def stop_err( msg ):
    sys.stderr.write( '%s\n' % msg )
    sys.exit()


def __main__():
    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '', '--input1', dest='input1', help=' ')
    parser.add_option( '', '--readsformat', dest='readsformat', help=' ')
    parser.add_option( '', '--own_file', dest='own_file', action='store', help='' )
    parser.add_option( '', '--indexes_path', dest='indexes_path', action='store', help='' )
    parser.add_option( '', '--seglen', dest='seglen', action='store', help='' )
    parser.add_option( '', '--chrfiles_dir', dest='chrfiles_dir', action='store', help='' )
    parser.add_option( '-o', '--output', dest='output', action='store', help='write output to file' )
    parser.add_option('-d', '--outputdir', dest='outputdir', action='store', default="")
    ( options, args ) = parser.parse_args()
    #print options
    #print args
    outdir=options.outputdir
    #print outdir
    outfilename=outdir+'/'+'alignments.sam'
    
    cmd = 'python /opt/rnaseq/MapSplice_1.15.2/bin/mapsplice_segments.py'
    cmd+=' -Q '+options.readsformat
    cmd+=' -o '+options.outputdir
    cmd+=' -c '+options.chrfiles_dir
    """ 
    A comma separated (no blank space) list of FASTA or FASTQ read files (including path)
                Notes:
                For paired-end reads, the order should be as follows: reads1_end1,reads1_end2,reads2_end1,read2_end2...
                For two ends from the same read, the read names should be in the following format: read_base_name/1 and read_base_name/2
                -The read_base_name should be the same for two ends
                Format constraint: Reads names after @ or > should not contain a blank space or tab 
    """
    cmd+=' -u '+options.input1
    for i, arg in enumerate( args ):
        print i,arg
        #arg is input_file_namess
        cmd += ",%s" % arg
    if options.indexes_path:
        cmd+=' -B '+options.indexes_path
    elif options.own_file:
        cmd+=' -B '+options.own_file
    cmd+=' -L '+options.seglen

    
    # Debugging.
    #mv can overide current files with the same name
    #Moved /opt/galaxy-dist/database/job_working_directory/000/331/galaxy_dataset_554.dat to /opt/galaxy-dist/database/files/000/dataset_554.dat
    if options.output.find('job_working_directory')!=-1:
        outpath=options.output.replace('job_working_directory','files')
        outpath=outpath.replace('galaxy_dataset','dataset')
        #re.sub(pattern, repl, string, count=0, flags=0)
        #Return the string obtained by replacing the leftmost non-overlapping occurrences of pattern in string by the replacement repl
        pattern=r'(\/\d+)(\/\d+)'
        outpath=re.sub(pattern,'\g<1>',outpath)
    else:
        outpath=options.output
    
    cmd +=';  wait; cp '+outfilename+' '+options.output 
    #galaxy will mv the output files from job_working_directory to files after finishing the run
    cmd +=';  wait; cp '+outfilename+' '+outpath

    print cmd
                
    try:
        #
        # Run command.
        #
        tmp_name = tempfile.NamedTemporaryFile( dir="." ).name
        tmp_stderr = open( tmp_name, 'wb' )
        #OK to use multiple job handlers
        proc = subprocess.Popen( args=cmd, shell=True, stderr=tmp_stderr.fileno() )
        #This will deadlock when using stdout=PIPE and/or stderr=PIPE and 
        #the child process generates enough output to a pipe 
        #such that it blocks waiting for the OS pipe buffer to accept more data. Use communicate() to avoid that
        returncode = proc.wait()
        #returncode = subprocess.call( args=cmd, shell=True )
        #streamdata = proc.communicate()
        #print streamdata
        #returncode = proc.returncode
        tmp_stderr.close()
        
        # Error checking.
        if returncode != 0:
            raise Exception, "return code = %i" % returncode
           
    except Exception, e:
        # Read stderr so that it can be reported:
        tmp_stderr = open( tmp_name, 'rb' )
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
        
        stop_err( 'Error running MapSplice.\n%s\n' % ( str( e ) ) )
             
    
    
if __name__=="__main__": __main__()
