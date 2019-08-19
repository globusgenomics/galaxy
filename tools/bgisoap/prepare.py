"""
prepare.py
A wrapper script for prepare
Copyright Peter Li peter@gigasciencejournal.com
"""

import sys, optparse, os, tempfile, shutil, subprocess

def stop_err(msg):
    sys.stderr.write('%s\n' % msg)
    sys.exit()


def cleanup_before_exit(tmp_dir):
    if tmp_dir and os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)


def __main__():
    #Parse command line
    parser = optparse.OptionParser()

    #Input
    parser.add_option('', '--input_contigs', dest='input_contigs', help='Input contings FASTA')
    parser.add_option('', '--kmer_length', dest='kmer_length', help='Kmer length')
    parser.add_option('', '--output_prefix', dest='output_prefix', help='Prefix of output')

    #Outputs
    parser.add_option("-o", "--outfile", dest='outfile', help="Sorted tabular data")
    opts, args = parser.parse_args()

    #Need to write inputs to a temporary directory
    dirpath = tempfile.mkdtemp()

    #Set up command line call
    cmd = "prepare -c %s -K %s -g %s" % (opts.input_contigs, opts.kmer_length, opts.output_prefix)
    print cmd

    #Run
    try:
        tmp_out_file = tempfile.NamedTemporaryFile(dir=dirpath).name
        tmp_stdout = open(tmp_out_file, 'wb')
        tmp_err_file = tempfile.NamedTemporaryFile(dir=dirpath).name
        tmp_stderr = open(tmp_err_file, 'wb')

        #Perform SOAPsnp call
        print "Doing prepare call..."
        soap_proc = subprocess.Popen(args=cmd, shell=True, cwd=dirpath, stdout=tmp_stdout, stderr=tmp_stderr)
        returncode = soap_proc.wait()
        tmp_stderr.close()
        # get stderr, allowing for case where it's very large
        tmp_stderr = open(tmp_err_file, 'rb')
        stderr = ''
        buffsize = 1048576
        try:
            while True:
                stderr += tmp_stderr.read(buffsize)
                if not stderr or len(stderr) % buffsize != 0:
                    break
        except OverflowError:
            pass
        tmp_stdout.close()
        tmp_stderr.close()
        if returncode != 0:
            raise Exception, stderr

    # TODO: look for errors in program output.
    except Exception, e:
        raise Exception, 'Problem performing prepare process ' + str(e)

    #Clean up temp files
    cleanup_before_exit(dirpath)
    stop_err('Error in running prepare from (%s), %s' % (opts.infile, str(e)))

    #Check results in output file
    if os.path.getsize(opts.outfile) > 0:
        sys.stdout.write('Status complete')
    else:
        stop_err("The output is empty")

if __name__ == "__main__": __main__()
