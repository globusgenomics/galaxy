"""
soap1.py
A wrapper script for SOAP1
Copyright Peter Li peter@gigasciencejournal.com
"""

import sys, optparse, os, tempfile, shutil, subprocess, ConfigParser

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
    parser.add_option("-d", "--ref_seq", dest="ref_seq", help="A reference genome in FASTA format")
    parser.add_option("", "--analysis_settings_type", dest="analysis_settings_type")
    parser.add_option("", "--default_full_settings_type", dest="default_full_settings_type")
    #Single-end params
    parser.add_option("-a", "--forward_set", dest="forward_set", help="Forward set of reads")
    #Paired-end params
    parser.add_option("-b", "--reverse_set", dest="reverse_set", help="Reverse set of reads for paired end analysis")
    parser.add_option("-m", "--min_insert_size", dest="min_insert_size", help="Minimal insert size allowed")
    parser.add_option("-x", "--max_insert_size", dest="max_insert_size", help="Maximum insert size allowed")
    #Custom params
    parser.add_option("-s", "--seed_size", dest="seed_size", help="Seed size")
    parser.add_option("-v", "--max_mismatches", dest="max_mismatches", help="Maximum number of mismatches allowed on a read")
    parser.add_option("-g", "--max_gap_size", dest="max_gap_size", help="Maximum gap size allowed on a read")
    parser.add_option("-w", "--max_best_hits", dest="max_best_hits", help="Maximum number of equal best hits to count")
    parser.add_option("-e", "--gap_exist", dest="gap_exist", help="Will not allow gap exist inside n-bp edge of a read")
    parser.add_option("-z", "--initial_quality", dest="initial_quality", help="Illumina uses '@', Sanger Institute uses '!'")
    parser.add_option("-c", "--trim", dest="trim", help="Trim 3' end due to low quality")
    parser.add_option("-f", "--filter", dest="filter", help="Filter low-quality reads containing >n Ns")
    parser.add_option("-r", "--report_repeats", dest="report_repeats", help="Report repeat hits")
    parser.add_option("-t", "--read_id", dest="read_id", help="Read ID in output file")
    parser.add_option("-n", "--ref_chain_align", dest="ref_chain_align", help="Do alignment on which reference chain")
    parser.add_option("-p", "--num_processors", dest="num_processors", help="Number of processors to use")
    #Outputs
    parser.add_option("-o", "--alignment_out", dest='alignment_out', help="An alignment of paired reads on a reference sequence")
    parser.add_option("-2", "--unpaired_alignment_out", dest='unpaired_alignment_out', help="Unpaired alignment hits")
    opts, args = parser.parse_args()

    #Create temp directory
    tmp_dir = tempfile.mkdtemp()

    #Set up command line call
    if opts.analysis_settings_type == "single" and opts.default_full_settings_type == "default":
        cmd = "soap1 -d %s -a %s -o %s" % (opts.ref_seq, opts.forward_set, opts.alignment_out)
    elif  opts.analysis_settings_type == "paired" and opts.default_full_settings_type == "default":
        cmd = "soap1 -d %s -a %s -b %s -o %s -2 %s -m %s -x %s" % (opts.ref_seq, opts.forward_set, opts.reverse_set, opts.alignment_out, opts.unpaired_alignment_out, opts.min_insert_size, opts.max_insert_size)
    elif opts.analysis_settings_type == "single" and opts.default_full_settings_type == "full":
        cmd = "soap1 -d %s -a %s -o %s -s %s -v %s -g %s -w %s -e %s -z %s -c %s -f %s -r %s -t %s -n %s -p %s" % (opts.ref_seq, opts.forward_set, opts.alignment_out, opts.seed_size, opts.max_mismatches, opts.max_gap_size, opts.max_best_hits, opts.gap_exist, opts.initial_quality, opts.trim, opts.filter, opts.report_repeats, opts.read_id, opts.ref_chain_align, opts.num_processors)
    elif opts.analysis_settings_type == "paired" and opts.default_full_settings_type == "full":
        cmd = "soap1 -d %s -a %s -b %s -o %s -2 %s -m %s -x %s -s %s -v %s -g %s -w %s -e %s -z %s -c %s -f %s -r %s -t %s -n %s -p %s" % (opts.ref_seq, opts.forward_set, opts.reverse_set, opts.alignment_out, opts.unpaired_alignment_out, opts.min_insert_size, opts.max_insert_size, opts.seed_size, opts.max_mismatches, opts.max_gap_size, opts.max_best_hits, opts.gap_exist, opts.initial_quality, opts.trim, opts.filter, opts.report_repeats, opts.read_id, opts.ref_chain_align, opts.num_processors)

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
        stop_err('Error in running soap1 from (%s), %s' % (opts.alignment_out, str(e)))

    #Clean up temp files
    cleanup_before_exit(tmp_dir)
    #Check results in output file
    if os.path.getsize(opts.alignment_out) > 0:
        sys.stdout.write('Status complete')
    else:
        stop_err("The output is empty")

if __name__ == "__main__": __main__()
