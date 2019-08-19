"""
msort.py
A wrapper script for msort
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
    parser.add_option('', '--infile', dest='infile', help='Input tabular data')
    parser.add_option('-k', '--sort_specified_fields', dest='sort_specified_fields', help='Specify fields for sorting')

    #Custom params
    parser.add_option("", "--default_full_settings_type", dest="default_full_settings_type")
    parser.add_option("-l", "--line_delimiters", dest="line_delimiters", help="Specify line delimiters")
    parser.add_option('-L', '--num_lines_per_block', dest='num_lines_per_block', help='Specify multiple lines as blocks')
    parser.add_option('-t', '--field_delimiters', dest='field_delimiters', help='Specify delimiters for fields')
    parser.add_option('-r', '--reverse_sort', dest='reverse_sort', help='Sort input data in reverse')
    parser.add_option('-f', '--ignore_char_case', dest='ignore_char_case', help='Ignore character cases')
    parser.add_option('-n', '--treat_fields_as_num', dest='treat_fields_as_num', help='Treat fields as numbers')
    parser.add_option('-m', '--treat_fields_as_num_or_string', dest='treat_fields_as_num_or_string', help='Treat fields as numbers or strings')

    #Outputs
    parser.add_option("-o", "--outfile", dest='outfile', help="Sorted tabular data")
    opts, args = parser.parse_args()

    #Need to write inputs to a temporary directory
    dirpath = tempfile.mkdtemp()
    print opts.default_full_settings_type
    print opts.line_delimiters
    #Set up command line call
    #Need to fix to account for user configuration at run time
    if opts.default_full_settings_type == "default":
        cmd = "msort -k %s %s > %s" % (opts.sort_specified_fields, opts.infile, opts.outfile)
    elif opts.default_full_settings_type == "full":
        cmd = "msort %s -l %s -L %s -t %s -k %s > %s" % (opts.infile, opts.line_delimiters, opts.num_lines_per_block, opts.field_delimiters, opts.sort_specified_fields, opts.outfile)
    print cmd

    #Run
    try:
        tmp_out_file = tempfile.NamedTemporaryFile(dir=dirpath).name
        tmp_stdout = open(tmp_out_file, 'wb')
        tmp_err_file = tempfile.NamedTemporaryFile(dir=dirpath).name
        tmp_stderr = open(tmp_err_file, 'wb')

        #Perform SOAPsnp call
        print "Doing msort call..."
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
        raise Exception, 'Problem performing msort process ' + str(e)

    #Clean up temp files
    cleanup_before_exit(dirpath)
    #12/12/2012, commented out by Huayan
    #This will cause an error
    #stop_err('Error in running msort from (%s), %s' % (opts.infile, str(e)))

    #Check results in output file
    if os.path.getsize(opts.outfile) > 0:
        sys.stdout.write('Status complete')
    else:
        stop_err("The output is empty")

if __name__ == "__main__": __main__()
