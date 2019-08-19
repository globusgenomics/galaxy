"""
soap2.py
A wrapper script for SOAP2
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
    parser.add_option('', '--fileSource', dest='fileSource', help='Source of reference dataset index')
    parser.add_option('', '--ref', dest='ref', help='The reference genome to use or index' )
    parser.add_option('', '--dbkey', dest='dbkey', help='Dbkey for reference genome' )
    parser.add_option('', '--do_not_build_index', dest='do_not_build_index', action='store_true', help="Don't build index" )

    parser.add_option("", "--analysis_settings_type", dest="analysis_settings_type")
    parser.add_option("", "--default_full_settings_type", dest="default_full_settings_type")
    #Single-end params
    parser.add_option("-a", "--forward_set", dest="forward_set", help="Forward set of reads")
    #Paired-end params
    parser.add_option("-b", "--reverse_set", dest="reverse_set", help="Reverse set of reads for paired end analysis")
    parser.add_option("-m", "--min_insert_size", dest="min_insert_size", help="Minimal insert size allowed")
    parser.add_option("-x", "--max_insert_size", dest="max_insert_size", help="Maximum insert size allowed")
    #Custom params
    parser.add_option("-n", "--filter", dest="filter", help="Filter low-quality reads containing >n Ns")
    parser.add_option("-t", "--read_id", dest="read_id", help="Read ID in output file")
    parser.add_option("-r", "--report_repeats", dest="report_repeats", help="Report repeat hits")
    parser.add_option("-R", "--long_insert_align", dest="long_insert_align", help="alignment for long insert size(>=  2k  bps)  PE  data, [none] FR alignment")
    parser.add_option("-l", "--high_error_rate", dest="high_error_rate", help="For  long  reads  with  high  error rate at 3'-end, those can't align whole length, then  first  align  5'  INT  bp subsequence as a seed, [256] use whole length of the read")
    parser.add_option("-v", "--allow_all_mismatches", dest="allow_all_mismatches", help="Totally allowed mismatches in one read")
    parser.add_option("-M", "--match_mode", dest="match_mode", help="Match mode for each read or the seed part of read, which shouldn't contain more than 2 mismaches")
    parser.add_option("-p", "--num_threads", dest="num_threads", help="Multithreads, n threads")

    #Outputs
    parser.add_option("-o", "--alignment_out", dest='alignment_out', help="An alignment of paired reads on a reference sequence")
    parser.add_option("-2", "--unpaired_alignment_out", dest='unpaired_alignment_out', help="Unpaired alignment hits")
    opts, args = parser.parse_args()

    #Create temp directory
    tmp_dir = tempfile.mkdtemp()
    print tmp_dir

    #Get reference index directory
    ref_index_filename = opts.ref

    #Set up command line call
    if opts.analysis_settings_type == "single" and opts.default_full_settings_type == "default":
        soap_cmd = "soap2 -a %s -D " %(opts.forward_set)+ ref_index_filename + " -o %s" % ( opts.alignment_out)
    elif  opts.analysis_settings_type == "paired" and opts.default_full_settings_type == "default":
        soap_cmd = "soap2 -a %s -b %s -D " %  (opts.forward_set, opts.reverse_set) + ref_index_filename + " -o %s -2 %s -m %s -x %s" % (opts.alignment_out, opts.unpaired_alignment_out, opts.min_insert_size, opts.max_insert_size)
    elif opts.analysis_settings_type == "single" and opts.default_full_settings_type == "full":
        soap_cmd = "soap2 -a %s -D " %(opts.forward_set) + ref_index_filename + " -o %s -n %s -t %s -r %s -v %s -M %s -p %s" % (opts.alignment_out, opts.filter, opts.read_id, opts.report_repeats, opts.allow_all_mismatches, opts.match_mode, opts.num_threads)
    elif opts.analysis_settings_type == "paired" and opts.default_full_settings_type == "full":
        soap_cmd = "soap2  -a %s -b %s -D " %(opts.forward_set, opts.reverse_set) + ref_index_filename + " -o %s -2 %s -m %s -x %s -n %s -t %s -r %s -v %s -M %s -p %s" % (opts.alignment_out, opts.unpaired_alignment_out, opts.min_insert_size, opts.max_insert_size, opts.filter, opts.read_id, opts.report_repeats, opts.allow_all_mismatches, opts.match_mode, opts.num_threads)

    print soap_cmd

    #Need to check format of reference sequence
    #Has to be in fasta format

    #Run
    try:
        tmp_out_file = tempfile.NamedTemporaryFile(dir=tmp_dir).name
        tmp_stdout = open(tmp_out_file, 'wb')
        tmp_err_file = tempfile.NamedTemporaryFile().name
        tmp_stderr = open(tmp_err_file, 'wb')

        #Perform SOAP2 call
        print "Doing SOAP2 call..."
        soap_proc = subprocess.Popen(args=soap_cmd, shell=True, cwd=".", stdout=tmp_stdout, stderr=tmp_stderr)
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
        #Clean up temp files
        cleanup_before_exit(tmp_dir)
        stop_err('Error in running soap1 from (%s), %s' % (opts.alignment_out, str(e)))

    #Clean up temp files
    cleanup_before_exit(tmp_dir)
    #Check results in output file
    if os.path.getsize(opts.alignment_out) > 0:
        sys.stdout.write('Status complete')
    else:
        stop_err("The output is empty")

if __name__ == "__main__": __main__()
