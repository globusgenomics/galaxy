#/usr/bin/python

import optparse, sys, subprocess, tempfile
import shutil, os

CHUNK_SIZE = 2**20 #1mbi

def cleanup_before_exit( tmp_dir ):
    if tmp_dir and os.path.exists( tmp_dir ):
        shutil.rmtree( tmp_dir )

def __main__():
    parser = optparse.OptionParser()
    parser.add_option( '-p', '--pass_through', dest='pass_through_options', action='append', type="string", help='These options are passed through directly to GATK, without any modification.' )
    parser.add_option('', '--of', dest='output_format', action='store', type="string", help='output format')
    parser.add_option('-i', '--input-files', dest='input_files', action='store', type="string", help='input files')
    parser.add_option( '-o', '--output-file', dest='output_file', action='store', type="string", default=None, help='bed output' )

    (options, args) = parser.parse_args()
    pass_through = ""
    if options.pass_through_options:
        pass_through = ' '.join( options.pass_through_options )

    cmd = "fseq %s -of %s %s" % (pass_through, options.output_format, options.input_files)
    print cmd
    tmp_dir = tempfile.mkdtemp( prefix='tmp-tool-' )

    stdout = tempfile.NamedTemporaryFile( prefix="fseq-stdout-", dir=tmp_dir )
    stderr = tempfile.NamedTemporaryFile( prefix="fseq-stderr-", dir=tmp_dir )
    proc = subprocess.Popen( args=cmd, stdout=stdout, stderr=stderr, shell=True, cwd=tmp_dir )
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
    stdout.close()
#    print options.output_format
    print tmp_dir
    if options.output_format == 'bed':
        os.system("cat %s/*.bed > %s/output.bed" % (tmp_dir, tmp_dir)) 
#        print "IAM HERE: %s/output.bed" % tmp_dir, options.output_file
        shutil.copy("%s/output.bed" % tmp_dir, options.output_file)
    else: 
        os.system("cat %s/*.wig > %s/output.wig" % (tmp_dir, tmp_dir))
#        print "%s/output.wig" % tmp_dir, options.output_file
        shutil.copy("%s/output.wig" % tmp_dir, options.output_file)

    cleanup_before_exit( tmp_dir )

if __name__ == "__main__": __main__()
