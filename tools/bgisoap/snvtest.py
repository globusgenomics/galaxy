"""
snvtest.py
A wrapper script for SNVtest
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
    #Generic input params
    parser.add_option("", "--input1", dest="input1", help="First input pileup file")
    parser.add_option("", "--input2", dest="input2", help="Second input pileup file")

    parser.add_option("", "--default_full_settings_type", dest="default_full_settings_type")
    #Custom params
    parser.add_option("-d", "--max_depth", dest="max_depth", help="Maximum depth requirement, sites with higher depth would be filtered")
    parser.add_option("-p", "--pvalue", dest="pvalue", help="p-value")
    parser.add_option("-q", "--cutoff", dest="cutoff", help="Quality cutoff for bases")

    #Outputs
    parser.add_option("-o", "--output", dest='output', help="Results from SNVtest process")
    opts, args = parser.parse_args()

    #Create temp directory
    tmp_dir = tempfile.mkdtemp()

    #Set up command line call
    if opts.default_full_settings_type == "default":
        cmd = "SNVtest %s %s -o %s" % (opts.input1, opts.input2, opts.output)
    elif opts.default_full_settings_type == "full":
        cmd = "SNVtest %s %s -o %s -d %s -p %s -q %s" % (opts.input1, opts.input2, opts.output, opts.max_depth, opts.pvalue, opts.cutoff)

    print cmd

    #Run
    try:
        tmp_out = tempfile.NamedTemporaryFile(dir=tmp_dir).name
        tmp_stdout = open(tmp_out, 'wb')
        tmp_err = tempfile.NamedTemporaryFile().name
        tmp_stderr = open(tmp_err, 'wb')
        proc = subprocess.Popen( args=cmd, shell=True, cwd=".", stdout=tmp_stdout, stderr=tmp_stderr )
        returncode = proc.wait()
        tmp_stderr.close()
        # get stderr, allowing for case where it's very large
        tmp_stderr = open(tmp_err, 'rb')
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
        #Clean up temp files
        cleanup_before_exit(tmp_dir)
        stop_err('Error in running SNVtest from (%s), %s' % (opts.output, str(e)))

    #Clean up temp files
    cleanup_before_exit(tmp_dir)
    #Check results in output file
    if os.path.getsize(opts.output) > 0:
        sys.stdout.write('Status complete')
    else:
        stop_err("The output is empty")

if __name__ == "__main__": __main__()
