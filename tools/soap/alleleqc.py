"""
alleleqc.py
A wrapper script for alleleQC
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
    parser.add_option("", "--pileup_input", dest="pileup_input", help="A pileup file")

    parser.add_option("", "--default_full_settings_type", dest="default_full_settings_type")
    #Custom params
    parser.add_option("-m", "--min_depth", dest="min_depth")
    parser.add_option("-x", "--max_depth", dest="max_depth")
    parser.add_option("-l", "--test_sites", dest="test_sites")

    #Outputs
    parser.add_option("-o", "--alleleqc_output", dest='alleleqc_output', help="alleleqc output file from alleleQC process")
    opts, args = parser.parse_args()

    #Create temp directory
    tmp_dir = tempfile.mkdtemp()

    tmp_output = tempfile.NamedTemporaryFile(dir=tmp_dir).name+".gz"

    #Set up command line call
    if opts.default_full_settings_type == "default":
        cmd = "alleleQC %s %s" % (opts.pileup_input, tmp_output)
    elif opts.default_full_settings_type == "full":
        cmd = "alleleQC %s %s -m %s -x %s -l %s" % (opts.pileup_input, opts.tmp_output, opts.min_depth, opts.max_depth, opts.test_sites)

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
        stop_err('Error in running alleleQC from (%s), %s' % (tmp_output, str(e)))

    #pass the tmp_output to opts.alleleqc_output
    gz_output = open(opts.alleleqc_output, 'wb')
    with open(tmp_output, mode='rb') as file: # b is important -> binary
        fileContent = file.read()
        gz_output.write(fileContent)
    gz_output.close()
    file.close()
    

    #Clean up temp files
    cleanup_before_exit(tmp_dir)
    #Check results in output file
    if os.path.getsize(opts.alleleqc_output) > 0:
        sys.stdout.write('Status complete')
    else:
        stop_err("The output is empty")

if __name__ == "__main__": __main__()
