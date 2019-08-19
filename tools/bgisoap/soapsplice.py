"""
A wrapper script for SOAPsplice
Copyright Peter Li peter@gigasciencejournal.com

Runs SOAPsplice on paired-end data. Need to ask Ruibang if it works on single-end data
Works with SOAPsplice version 1.9.

usage: soapsplice.py [options]

"""

import optparse, os, shutil, subprocess, sys, tempfile

def get_version():
    sys.stdout.write('SOAPsplice version 1.9')


def stop_err( msg ):
    sys.stderr.write('%s\n' % msg)
    sys.exit()


def __main__():
    #Parse command line input parameters
    parser = optparse.OptionParser()
    parser.add_option('', '--fileSource', dest='fileSource', help='Source of reference dataset index')
    parser.add_option('', '--ref', dest='ref', help='The reference genome to use or index')
    parser.add_option('', '--dbkey', dest='dbkey', help='Dbkey for reference genome')
    parser.add_option('', '--do_not_build_index', dest='do_not_build_index', action='store_true',
        help="Don't build index")
    parser.add_option('-1', '--input1', dest='fastq', help='The forward fastq file to use for the mapping')
    parser.add_option('-2', '--input2', dest='rfastq',
        help='The reverse fastq file to use for mapping if paired-end data')
    parser.add_option('-o', '--output_dir_prefix', dest='output_dir_prefix',
        help='Prefix of output files, can also include a file directory')
    parser.add_option('-I', '--insert_length', dest='insert_length', help='Insert length of paired-end reads')
    #Remaining custom parameters
    parser.add_option("", "--analysis_settings_type", dest="analysis_settings_type")
    parser.add_option("", "--default_full_settings_type", dest="default_full_settings_type")
    parser.add_option('-p', '--num_threads', dest='num_threads', help='Number of threads to use')
    parser.add_option('-S', '--forward_reverse_both', dest='forward_reverse_both', help='Forward, reverse or both')
    parser.add_option('-m', '--max_mismatch', dest='max_mismatch', help='Maximum mismatch for one-segment alignment')
    parser.add_option('-g', '--max_indel', dest='max_indel', help='Maximum indel for one-segment alignment')
    parser.add_option('-i', '--ignore_tail_length', dest='ignore_tail_length',
        help='Length of tail that can be ignored in one-segment alignment')
    parser.add_option('-t', '--longest_gap_length', dest='longest_gap_length',
        help='Longest gap between two segments in two-segment alignment')
    parser.add_option('-a', '--shortest_segment_length', dest='shortest_segment_length',
        help='Shortest length of a segment in two-segment alignment')
    parser.add_option('-c', '--output_read_and_quality', dest='output_read_and_quality',
        help='Output read and quality information?')
    parser.add_option('-f', '--output_format', dest='output_format', help='Format of output files')
    parser.add_option('-s', '--set_mapq', dest='set_mapq', help='Set mapping quality value')
    parser.add_option('-q', '--input_quality_type', dest='input_quality_type', help='Input quality type in FASTQ file')
    parser.add_option('-L', '--max_distance_bet_paired_ends', dest='max_distance_bet_paired_ends',
        help='Maximum distance between paired-end reads')
    parser.add_option('-l', '--min_distance_bet_paired_ends', dest='min_distance_bet_paired_ends',
        help='Minimum distance between paired-end reads')
    parser.add_option('-j', '--output_junction_info', dest='output_junction_info', help='Output junction information?')
    #Individual outputs
    parser.add_option('', '--forward_2segs', dest='forward_2segs', help='Two-segment alignment result for forward reads')
    parser.add_option('', '--forward_out', dest='forward_out', help='One-segment alignment result for forward reads')
    parser.add_option('', '--reverse_2segs', dest='reverse_2segs', help='Two-segment alignment result for reverse reads')
    parser.add_option('', '--reverse_out', dest='reverse_out', help='One-segment alignment result for reverse reads')
    parser.add_option('', '--junction', dest='junction', help='Junction file')

    options, args = parser.parse_args()

    #Create temporary directory to hold results
    tmp_dir = tempfile.mkdtemp()
    print 'Temporary dir: ' + tmp_dir

    #Name of reference index
    ref_file_name = options.ref

    #Create soapsplice command call
    default_cmd = 'soapsplice -d %s -1 %s -2 %s -I %s -o %s' % (
        ref_file_name, options.fastq, options.rfastq, options.insert_length, tmp_dir + '/result')
    print default_cmd

    full_cmd = 'soapsplice -d %s -1 %s -2 %s -I %s -o %s -p %s -S %s -m %s -g %s -i %s -t %s -a %s -c %s -f %s -s %s -q %s -L %s -l %s -j %s' % (
        ref_file_name, options.fastq, options.rfastq, options.insert_length, tmp_dir, options.num_threads, options.forward_reverse_both, options.max_mismatch, options.max_indel, options.ignore_tail_length, options.longest_gap_length, options.shortest_segment_length, options.output_read_and_quality, options.output_format, options.set_mapq, options.input_quality_type, options.max_distance_bet_paired_ends, options.min_distance_bet_paired_ends, options.output_junction_info)
    print full_cmd

    #Perform soapsplice analysis
    buffsize = 1048576
    try:
        #Create file in temporary directory
        tmp = tempfile.NamedTemporaryFile(dir=tmp_dir).name
        #Open a stream to file
        tmp_stderr = open(tmp, 'wb')
        #Call soapsplice
        proc = subprocess.Popen(args=default_cmd, shell=True, cwd=tmp_dir, stderr=tmp_stderr.fileno())
        returncode = proc.wait()
        #Close stream
        tmp_stderr.close()
        #Get stderr, allowing for case where it's very large
        tmp_stderr = open(tmp, 'rb')
        stderr = ''
        try:
            while True:
                stderr += tmp_stderr.read(buffsize)
                if not stderr or len(stderr) % buffsize != 0:
                    break
        except OverflowError:
            pass
        tmp_stderr.close()
        if returncode != 0:
            raise Exception, stderr
    except Exception, e:
        raise Exception, 'Error aligning sequence. ' + str(e)

    #Check that there are results in the output file
    #    if os.path.getsize(tmp_input1_2segs) > 0:
    #        sys.stdout.write('SOAPsplice run on paired-end data')
    #    else:
    #        raise Exception, 'The output file is empty. You may simply have no matches, or there may be an error with your input file or settings.'

    #Read soapsplice results into outputs
    out_forward_2segs = open(options.forward_2segs, 'wb')
    file = open(tmp_dir + '/result_1.2Segs')
    for line in file:
        print line
        out_forward_2segs.write(line)
    out_forward_2segs.close()

    out_forward_out = open(options.forward_out, 'wb')
    file = open(tmp_dir + '/result_1.out')
    for line in file:
        print line
        out_forward_out.write(line)
    out_forward_out.close()

    out_reverse_2segs = open(options.reverse_2segs, 'wb')
    file = open(tmp_dir + '/result_2.2Segs')
    for line in file:
        print line
        out_reverse_2segs.write(line)
    out_reverse_2segs.close()

    out_reverse_out = open(options.reverse_out, 'wb')
    file = open(tmp_dir + '/result_2.out')
    for line in file:
        print line
        out_reverse_out.write(line)
    out_reverse_out.close()

    junction = open(options.junction, 'wb')
    file = open(tmp_dir + '/result.junc')
    for line in file:
        print line
        junction.write(line)
    junction.close()

    # clean up temp dir
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

if __name__ == "__main__": __main__()
