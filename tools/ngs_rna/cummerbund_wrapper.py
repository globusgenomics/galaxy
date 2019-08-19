#!/usr/bin/env python

### Runs "r_script" and generates a HTML report
### Inspired on cuffdiff_wrapper.py and gatk_wrapper.py
### Carlos Borroto <carlos.borroto@gmail.com>

import optparse, os, shutil, subprocess, sys, tempfile

def stop_err( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit()

def html_report_from_directory( html_out, dir ):
    html_out.write( '<html>\n<head>\n<title>Galaxy - cummeRbund Output</title>\n</head>\n<body>\n<p/>\n<ul>\n' )
    for fname in sorted( os.listdir( dir ) ):
        html_out.write(  '<li><a href="%s">%s</a></li>\n' % ( fname, fname ) )
        html_out.write(  '<img src="%s" alt="%s" style="width:304px;height:228px;">\n' % (fname, fname) )
    html_out.write( '</ul>\n</body>\n</html>\n' )

def __main__():
    #Parse Command Line
    parser = optparse.OptionParser()

    # wrapper options
    parser.add_option('', '--r-script', dest='r_script', help='R script')
    parser.add_option('', '--html-report-from-directory', dest='html_report_from_directory', type="string", nargs=2, help='"Target HTML File" "Directory"')
        
    (options, args) = parser.parse_args()
    
    (html_filename, html_dir) = options.html_report_from_directory
    
    # Make html report directory for output.
    os.mkdir( html_dir )
    
    # Make a tmp dir
    tmp_dir = tempfile.mkdtemp( prefix='tmp-cummeRbund-' )
    
    # Build command.
    cmd = ( "Rscript --vanilla %s" % options.r_script )
    
    # Debugging.
    print cmd

#liubo added, for test, look at the generated R script
#    shutil.copy2(options.r_script, '/nfs/software/galaxy/r_script')

    
    # Run command.
    try:
        tmp_name = tempfile.NamedTemporaryFile( dir=tmp_dir ).name
        tmp_stderr = open( tmp_name, 'wb' )
        proc = subprocess.Popen( args=cmd, shell=True, cwd=html_dir, stderr=tmp_stderr.fileno() )
        returncode = proc.wait()
        tmp_stderr.close()
        
        # Get stderr, allowing for case where it's very large.
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
        
        # Error checking.
        if returncode != 0:
            raise Exception, stderr
    except Exception, e:
        stop_err( 'Error running R script. ' + str( e ) )
    
    # write the html report
    html_report_from_directory( open( html_filename, 'wb' ), html_dir )
    
    # Clean up temp dirs
    if os.path.exists( tmp_dir ):
        shutil.rmtree( tmp_dir )

if __name__=="__main__": __main__()
