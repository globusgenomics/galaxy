#!/usr/bin/env python


import optparse, os, sys, subprocess, tempfile, shutil
#from galaxy import eggs
#import pkg_resources; pkg_resources.require( "bx-python" )
#from bx.cookbook import doc_optparse
#from galaxy import util

def stop_err( msg ):
    sys.stderr.write( '%s\n' % msg )
    sys.exit()

def __main__():
    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '', '--input', dest='input', help='The input BAM dataset' )
    parser.add_option( '', '--output', dest='output', help='The output BAM dataset' )
    parser.add_option( '', '--order', dest='order', help='The order to be sorted in' )
    ( options, args ) = parser.parse_args()

    # output version # of tool
    try:
        tmp = tempfile.NamedTemporaryFile().name
        tmp_stdout = open( tmp, 'wb' )
        proc = subprocess.Popen( args='samtools 2>&1', shell=True, stdout=tmp_stdout )
        tmp_stdout.close()
        returncode = proc.wait()
        stdout = None
        for line in open( tmp_stdout.name, 'rb' ):
            if line.lower().find( 'version' ) >= 0:
                stdout = line.strip()
                break
        if stdout:
            sys.stdout.write( 'Samtools %s\n' % stdout )
        else:
            raise Exception
    except:
        sys.stdout.write( 'Could not determine Samtools version\n' )

    try:
        tmp_dir = tempfile.mkdtemp()
        input_bam = "%s/input.bam" % tmp_dir
        os.symlink(options.input, input_bam )
        order_param = ""
        if options.order == "lexicographically":
            order_param = "-n"
        command = 'samtools sort -m 8G -@ 28 -O bam -T tempsort -o %s %s %s ' % ( options.output, order_param, input_bam )
        print command
        tmp_dir = tempfile.mkdtemp()
        tmp = tempfile.NamedTemporaryFile( dir=tmp_dir ).name
        tmp_stderr = open( tmp, 'wb' )
        proc = subprocess.Popen( args=command, shell=True, cwd=tmp_dir, stderr=tmp_stderr.fileno() )
        returncode = proc.wait()
        tmp_stderr.close()

#    subprocess.call(command,shell=True)
        #generate_bam = options.output + ".bam"

        #shutil.move(generate_bam,options.output)
    except:
        sys.stdout.write( 'Could not run Samtools sort\n' )
    



if __name__=="__main__": __main__()

