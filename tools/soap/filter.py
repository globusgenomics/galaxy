"""
filter.py
A wrapper script for Filter tool from SOAPdenovo2
Peter Li peter@gigasciencejournal.com
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
    parser.add_option("", "--read1", dest="read1", help="Forward set of reads")
    parser.add_option("", "--read2", dest="read2", help="Reverse set of reads")

#    parser.add_option("", "--analysis_settings_type", dest="analysis_settings_type")
    parser.add_option("", "--default_full_settings_type", dest="default_full_settings_type")
    #Custom params
    parser.add_option("-a", "--read1_5prime_trim_length", dest="read1_5prime_trim_length", help="Length of forward reads to be trimmed from 5' end")
    parser.add_option("-b", "--read1_3prime_trim_length", dest="read1_3prime_trim_length", help="Length of forward reads to be trimmed from 3' end")
    parser.add_option("-c", "--read2_5prime_trim_length", dest="read2_5prime_trim_length", help="Length of reverse reads to be trimmed from 5' end")
    parser.add_option("-d", "--read2_3prime_trim_length", dest="read2_3prime_trim_length", help="Length of reverse reads to be trimmed from 3' end")
    parser.add_option("-f", "--trim_flag", dest="trim_flag", help="Trim flag")
    parser.add_option("-q", "--quality_shift", dest="quality_shift", help="Quality shift value")
    parser.add_option("-m", "--reads_pair_number", dest="reads_pair_number", help="Number of read pairs in buffer")
    parser.add_option("-t", "--num_threads", dest="num_threads", help="Number of threads to use")
    parser.add_option("-B", "--filter_low_quality_bases", dest="filter_low_quality_bases", help="Filter low quality bases")
    parser.add_option("-l", "--library_insert_size", dest="library_insert_size", help="Library insert size")
    parser.add_option("-w", "--filter_percent_N_bases", dest="filter_percent_N_bases", help="Percentage cutoff of N bases to filter")
    parser.add_option("-y", "--filter_reads_with_adapter_seq", dest="filter_reads_with_adapter_seq", help="Filter reads containing adapter sequences")
    parser.add_option("-F", "--read1_adapter_seq", dest="read1_adapter_seq", help="Forward read adapter sequence")
    parser.add_option("-R", "--read2_adapter_seq", dest="read2_adapter_seq", help="Reverse read adapter sequence")
    parser.add_option("-z", "--filter_small_reads", dest="filter_small_reads", help="Filter small-sized reads")
    parser.add_option("-p", "--filter_PCR_duplications", dest="filter_PCR_duplications", help="Filter PCR duplications")
    parser.add_option("-g", "--compress_output_read_file", dest="compress_output_read_file", help="Compress output read files")

    #Outputs
    parser.add_option("", "--stat", dest='stat', help="Provides statistics on the reads being filtered")
    parser.add_option("", "--read1_clean", dest='read1_clean', help="Filtered forward set of reads")
    parser.add_option("", "--read2_clean", dest='read2_clean', help="Filtered reverse set of reads")
    opts, args = parser.parse_args()

    #Create temp directory
    tmp_dir = tempfile.mkdtemp()

    #Set up command line call
    if opts.default_full_settings_type == "default":
        cmd = "SOAPfilter_v2.0 %s %s %s %s %s" % (opts.read1, opts.read2, opts.stat, opts.read1_clean, opts.read2_clean)
    elif opts.default_full_settings_type == "full":
        #Create cmd string
        cmd = "SOAPfilter_v2.0 -t %s -m %s" % (opts.num_threads, opts.reads_pair_number)

        if opts.filter_reads_with_adapter_seq == "yes":
            cmd =  cmd + " -y -F %s -R %s" % (opts.read1_adapter_seq, opts.read2_adapter_seq)

        if opts.filter_small_reads == "yes":
            cmd =  cmd + " -z"

        if opts.filter_PCR_duplications == "yes":
            cmd =  cmd + " -p"

        if opts.compress_output_read_file == "yes":
            cmd =  cmd + " -g"

        cmd =  cmd + " -f %s" % (opts.trim_flag)

        if opts.read1_5prime_trim_length != 0:
            cmd =  cmd + " -a %s" % opts.read1_5prime_trim_length

        if opts.read1_3prime_trim_length != 0:
            cmd =  cmd + " -b %s" % (opts.read1_3prime_trim_length)

        if opts.read2_5prime_trim_length > 0:
            cmd =  cmd + " -c %s" % (opts.read2_5prime_trim_length)

        if opts.read2_3prime_trim_length > 0:
            cmd =  cmd + " -d %s" % (opts.read2_3prime_trim_length)

        cmd = cmd + " -q %s -B %s -l %s -w %s %s %s %s %s %s" % (opts.quality_shift, opts.filter_low_quality_bases, opts.library_insert_size, opts.filter_percent_N_bases, opts.read1, opts.read2, opts.stat, opts.read1_clean, opts.read2_clean)

    print cmd

    #Run
    try:
        tmp_out_file = tempfile.NamedTemporaryFile(dir=tmp_dir).name
        tmp_stdout = open(tmp_out_file, 'wb')
        tmp_err_file = tempfile.NamedTemporaryFile().name
        tmp_stderr = open(tmp_err_file, 'wb')

        #Perform Filter call
        print "Doing Filter call..."
        proc = subprocess.Popen(args=cmd, shell=True, cwd=tmp_dir, stdout=tmp_stdout, stderr=tmp_stderr)
        returncode = proc.wait()
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
        #Clean up temp files
        cleanup_before_exit(tmp_dir)
        stop_err('Error in running Filter from %s' % (str(e)))

    #Clean up temp files
    cleanup_before_exit(tmp_dir)
    #Check results in output file
    if os.path.getsize(opts.stat) > 0:
        sys.stdout.write('Status complete')
    else:
        stop_err("The output is empty")

if __name__ == "__main__": __main__()
