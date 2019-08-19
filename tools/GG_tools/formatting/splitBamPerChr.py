#!/usr/bin/env python
"""
Splits BAM data per chromosome.
usage: splitBamPerChr.py [options]
   -i BAM file
   -o output File name
   --files_path output for chr consensus files
"""

import time, optparse, os, sys, subprocess, tempfile, shutil

galhtmlprefix = """<?xml version="1.0" encoding="utf-8" ?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<meta name="generator" content="Galaxy tool output - see http://getgalaxy.org/" />
<title></title>
<link rel="stylesheet" href="/static/style/base.css" type="text/css" />
</head>
<body>
<div class="document">
"""
galhtmlpostfix = """</div></body></html>\n"""

def stop_err( msg ):
    sys.stderr.write( '%s\n' % msg )
    sys.exit()

def cleanup_before_exit( tmp_dir ):
    if tmp_dir and os.path.exists( tmp_dir ):
        shutil.rmtree( tmp_dir )

def create_bam_idx(bamFile, target_dir = None):
    if target_dir is None:
        target_dir = os.getcwd()

    index_name = "%s/%s.bai" % (target_dir, os.path.basename(bamFile))
    stderr_name = tempfile.NamedTemporaryFile( prefix = "bam_stderr" ).name
    command = 'samtools index %s %s' % ( bamFile, index_name )
    print command
    proc = subprocess.Popen( args=command, shell=True, stderr=open( stderr_name, 'wb' ) )
    exit_code = proc.wait()
    if exit_code:
        for line in open( stderr_name ):
            print >> sys.stderr, line
        os.unlink( stderr_name ) #clean up
        #cleanup_before_exit( tmp_dir )
        raise Exception( "Error indexing BAM file" )
    os.unlink( stderr_name ) #clean up

def __main__():
    print "\n\nStart time:"
    print time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))

    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '-i', '--input', dest='bamFile', help='The input BAM dataset' )
    parser.add_option( '', '--indexBam', dest='bamFileIndex', help='The BAM index file' )
    parser.add_option( '-o', '--output', dest='outputFile', help='Output file will be a log file pointing to the location of each of the file.' )
    parser.add_option( '', '--files_path', dest='extra_files_path', help="Location for chr consensus output" )
    ( options, args ) = parser.parse_args()

    # exit if input file empty
    if os.path.getsize( options.bamFile ) == 0:
        raise Exception, 'Initial BAM file empty'

    #print "Creating directory: " + args.files_path
    if not os.path.exists(options.extra_files_path):
        os.makedirs(options.extra_files_path)

    tmp_dir = tempfile.mkdtemp(prefix="bamIndex", dir=options.extra_files_path)

    # create link of input BAM file to the output directory
    basename = os.path.basename(options.bamFile)
    link_filename = os.path.join( tmp_dir, "%s.bam" % ( basename ) )
    os.symlink( options.bamFile, link_filename )

    # create a BAM index file
    if options.bamFileIndex:
        # create a link to the index file in the same directory as the bam file
        bai_link_filename = os.path.join(tmp_dir, "%s.bam.bai" % ( basename ) )
        os.symlink( options.bamFileIndex, bai_link_filename)
    else:
        create_bam_idx(link_filename, target_dir = tmp_dir)

    ## run the bamtools split cmd
    cmd = "samtools view -H %s | grep \"\@SQ\" | sed 's/^.*SN://g' | cut -f 1 | grep -v GL | grep -v NC_ |xargs -I {} -n 1 -P 32 sh -c \"samtools view -b %s {} > %s/{}.bam; samtools index %s/{}.bam\"" % (link_filename, link_filename, options.extra_files_path, options.extra_files_path)
    print cmd

    stderr_name = tempfile.NamedTemporaryFile( prefix = "cmd_stderr" ).name
    proc = subprocess.Popen( args=cmd, shell=True, stderr=open( stderr_name, 'wb' ) )
    exit_code = proc.wait()
    if exit_code:
        for line in open( stderr_name ):
            print >> sys.stderr, line
            os.unlink( stderr_name ) #clean up
            raise Exception( "Error running command " )
        os.unlink( stderr_name ) #clean up

    # delete the BAM link file
    os.unlink(link_filename)
    cleanup_before_exit(tmp_dir)

    # create an output file with the name of all files created and stored in the library, and name of library
    #Create HTML file
    f = open(options.outputFile,'w')
    f.write(galhtmlprefix)
    #f.write("<p>Files generated from BAM file: %s</p>\n" % (os.path.basename(options.bamFile))
    for dirfile in os.listdir("%s" % (options.extra_files_path)):
        f.write('<a href="%s">%s</a><br>\n' % (dirfile, dirfile))
    f.write(galhtmlpostfix)
    f.close()

    print "\n\nEnd time:"
    print time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))

if __name__=="__main__": __main__()

