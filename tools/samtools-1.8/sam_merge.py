#!/usr/bin/env python

"""
Merges any number of BAM files
usage: %prog [options]
    output1
    input1
    input2
    [input3[,input4[,input5[,...]]]]
"""

import optparse, os, subprocess, sys, tempfile, shutil

def stop_err( msg ):
    sys.stderr.write( '%s\n' % msg )
    sys.exit()

def __main__():
    parser = optparse.OptionParser()
    parser.add_option( '', '--rmdup', dest='rmdup', action='store_true', help='Remove duplicates' )
    (options, args) = parser.parse_args()
    #print args

    infiles = []
    infiles.append(args[1])
    outfile = args[0]
    if len( args ) < 3:
        stop_err( 'There are not enough files to merge' )
    infiles.extend(args[2:])
    print "INFILES: %s" % infiles
    #filenames = args[2:]
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

    tmp_merged_bam = tempfile.NamedTemporaryFile(prefix="merged_bam", suffix=".bam").name
    tmp_rmdup_bam = tempfile.NamedTemporaryFile(prefix="rmdup_bam", suffix=".bam").name
    cmd = 'samtools merge -@ 24 -f %s %s ' % ( tmp_merged_bam, ' '.join( infiles ) ) # Galaxy must be touching the output files before tool executed.  W/o -f, samtools merge complains about not being able to overwrite the existing output file.
    if options.rmdup:
        # build the header by taking the RG field from each of the BAM's header
        rgfield = []
        headerlines = []
        pgfield = []
        count = 0
        #infiles = filenames
        #infiles.append(infile)
        for i in infiles:
            pipe = subprocess.Popen('samtools view -H %s' % i, shell=True, stdout=subprocess.PIPE )
            for line in pipe.stdout:
                if line.startswith('@RG'):
                    rgfield.append(line)
                elif line.startswith('@PG') and count == 0:
                    pgfield.append(line)
                elif count == 0:
                    headerlines.append(line)
            count += 1
        tmp_header_sam = tempfile.NamedTemporaryFile(prefix="header", suffix=".sam").name
        f = open(tmp_header_sam, 'w')
        f.write('%s%s%s' % ( "".join(headerlines), "".join(rgfield), "".join(pgfield) ))
        f.close()
        cmd = 'samtools merge -@ 24 -h %s  -f - %s | samtools rmdup - %s' % (tmp_header_sam, ' '.join( infiles ), tmp_rmdup_bam)
    tmp = tempfile.NamedTemporaryFile().name
    try:
        tmp_stderr = open( tmp, 'wb' )
        print "command is: %s" % cmd
        proc = subprocess.Popen( args=cmd, shell=True, stderr=tmp_stderr.fileno() )
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
        if os.path.exists( tmp ):
            os.unlink( tmp )
    except Exception, e:
        if os.path.exists( tmp ):
            os.unlink( tmp )
        stop_err( 'Error running SAMtools merge tool\n' + str( e ) )
    # copy the output to the outfile location
    if options.rmdup:
        shutil.copy(tmp_rmdup_bam, outfile)
    else:
        shutil.copy(tmp_merged_bam, outfile)

    if os.path.getsize( outfile ) > 0:
        sys.stdout.write( '%s files merged.' % ( len( sys.argv ) - 2 ) )
    else:
        stop_err( 'The output file is empty, there may be an error with one of your input files.' )

if __name__ == "__main__" : __main__()
