#!/usr/bin/env python

"""
Classes encapsulating decypher tool.
James E Johnson - University of Minnesota
"""
import pkg_resources;
import logging, os, string, sys, tempfile, glob, shutil, types, urllib
import shlex, subprocess
from optparse import OptionParser, OptionGroup
from stat import *


log = logging.getLogger( __name__ )

assert sys.version_info[:2] >= ( 2, 4 )

def stop_err( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit()

def __main__():
    #Parse Command Line
    s = 'preppereads.py:  argv = %s\n' % (sys.argv)
    # print >> sys.stderr, s # so will appear as blurb for file
    argcnt = len(sys.argv)
    working_dir = sys.argv[1]
    input1 = sys.argv[2]
    input2 = sys.argv[3]
    outpe = sys.argv[4]
    #outsingletons = sys.argv[5]

    # create symlinks in a temp directory for the input files
    tmpdir = tempfile.mkdtemp(dir=working_dir)
    link1 = "%s/input1.fq" % tmpdir
    link2 = "%s/input2.fq" % tmpdir
    os.symlink(input1, link1)
    os.symlink(input2, link2)

#    cmdline = '/nfs/software/galaxy/tools/sr_assembly//prepare_pe_reads_for_velvet.sh %s %s %s %s %s %s > /dev/null' % (input1, '1.fastq', input2, '2.fastq', 'single.fq', 'pe.fq')
    cmdline = '/nfs/software/galaxy/tools/sr_assembly//prepare_pe_reads_for_velvet.sh %s %s %s %s %s %s > /dev/null' % (link1, '1.fastq', link2, '2.fastq', 'single.fq', 'pe.fq')
    print cmdline
    #print >> sys.stderr, cmdline # so will appear as blurb for file
    try:
        proc = subprocess.Popen( args=cmdline, shell=True, stderr=subprocess.PIPE )
        returncode = proc.wait()
        # get stderr, allowing for case where it's very large
        stderr = ''
        buffsize = 1048576
        try:
            while True:
                stderr += proc.stderr.read( buffsize )
                if not stderr or len( stderr ) % buffsize != 0:
                    break
        except OverflowError:
            pass
        if returncode != 0:
            raise Exception, stderr
    except Exception, e:
        stop_err( 'Error running preppereads.sh ' + str( e ) )

    out = open(outpe,'w')
    #outpe_path = os.path.join(working_dir,'')
    for line in open('pe.fq'):
        out.write( "%s" % (line) )
    out.close()

#    out = open(outsingletons,'w')
#    #singletons_path = os.path.join(working_dir,'')
#    for line in open('single.fq'):
#        out.write( "%s" % (line) )
#    out.close()

if __name__ == "__main__": __main__()
