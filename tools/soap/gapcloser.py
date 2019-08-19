"""
gapcloser.py
A wrapper script for GapCloser
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
    parser.add_option('-a', '--scaffold', dest='scaffold', help='Scaffold file')
    parser.add_option('-b', '--configinfo', dest='configinfo', help='Configuration file used by SOAPdenovo')

    #Custom params
    parser.add_option('', "--default_full_settings_type", dest="default_full_settings_type")
    parser.add_option('-l', '--max_read_len', dest='max_read_len', help='Maximum read length')
    parser.add_option('-p', '--overlap_param', dest='overlap_param', help='Overlap parameter')
    parser.add_option('-t', '--thread_num', dest='thread_num', help='Number of threads')

    #Outputs
    parser.add_option("-o", "--outfile", dest='outfile', help="Scaffold with gaps filled")
    parser.add_option('', "--fill", dest='fill', help="Information about gaps in scaffold")
    opts, args = parser.parse_args()

    #Need to write inputs to a temporary directory
    dirpath = tempfile.mkdtemp()
    scaffold_data = open(opts.scaffold, 'r')
    scaffold_file = open(dirpath + "/input.scaffold", "w")
    for line in scaffold_data:
        scaffold_file.write(line)
    scaffold_data.close()
    scaffold_file.close()

    configinfo_data = open(opts.configinfo, 'r')
    configinfo_file = open(dirpath + "/input.configinfo", "w")
    for line in configinfo_data:
        configinfo_file.write(line)
    configinfo_data.close()
    configinfo_file.close()

    #Set up command line call
    #Need to fix to account for user configuration at run time
    if opts.default_full_settings_type == "default":
        cmd = "GapCloser -a %s -b %s -o %s" % (dirpath + "/input.scaffold", dirpath + "/input.configinfo", dirpath + "/output.gapcloser")
    elif opts.default_full_settings_type == "full":
        cmd = "GapCloser -a %s -b %s -o %s -l %s -p %s -t %s" % (dirpath + "/input.scaffold", dirpath + "/input.configinfo", "output.gapcloser", opts.max_read_len, opts.overlap_param, opts.thread_num)

    print cmd

    #Run
    try:
        tmp_out_file = tempfile.NamedTemporaryFile(dir=dirpath).name
        tmp_stdout = open(tmp_out_file, 'wb')
        tmp_err_file = tempfile.NamedTemporaryFile().name
        tmp_stderr = open(tmp_err_file, 'wb')

        #Perform SOAPsnp call
        print "Doing GapCloser call..."
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
        raise Exception, 'Problem performing GapCloser process ' + str(e)

    #Read scaffold file into its output
    outfile_out = open(opts.outfile, 'w')
    file = open(dirpath + "/output.gapcloser")
    for line in file:
        outfile_out.write(line)
    outfile_out.close()
    file.close()

    #Read fill file into its output
    fill_out = open(opts.fill, 'w')
    file = open(dirpath + "/output.gapcloser.fill")
    for line in file:
        fill_out.write(line)
    fill_out.close()
    file.close()

    #Clean up temp files
    cleanup_before_exit(dirpath)
    stop_err('Error in running GapCloser from (%s), %s' % (opts.scaffold, str(e)))

    #Check results in output file
    if os.path.getsize(opts.scaffold) > 0:
        sys.stdout.write('Status complete')
    else:
        stop_err("The output is empty")

if __name__ == "__main__": __main__()