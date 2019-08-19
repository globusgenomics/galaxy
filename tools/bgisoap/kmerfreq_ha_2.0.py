"""
kmerfreq-ha-2.0.py
A wrapper script for kmerfreq
Peter Li - GigaScience, BGI-HK
"""

import optparse, os, shutil, subprocess, sys, tempfile

def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit()


def cleanup_before_exit(tmp_dir):
    if tmp_dir and os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)


def main():
    #Parse command line
    parser = optparse.OptionParser()
    #Make list of params
    parser.add_option("", "--format_of_data", action="append", type="string", dest="format_of_data", help="Format of data")
    #Data inputs
    parser.add_option("", "--paired_fastq_input1", action="append", type="string", dest="paired_fastq_input1_list")
    parser.add_option("", "--paired_fastq_input2", action="append", type="string", dest="paired_fastq_input2_list")
    parser.add_option("", "--paired_fasta_input1", action="append", type="string", dest="paired_fasta_input1_list")
    parser.add_option("", "--paired_fasta_input2", action="append", type="string", dest="paired_fasta_input2_list")

    parser.add_option("", "--default_full_settings_type", dest="default_full_settings_type")

    #Custom params
    parser.add_option("-k", "--kmer_size", dest="kmer_size")
    parser.add_option("-r", "--read_length", dest="read_length")
    parser.add_option("-a", "--ignore_first", dest="ignore_first")
    parser.add_option("-d", "--ignore_last", dest="ignore_last")
    parser.add_option("-g", "--use_num_bases", dest="use_num_bases")
    #Don't need to expose -l input read list file param
    #Don't need to expose output prefix
    parser.add_option("-i", "--hash_size", dest="hash_size")
    parser.add_option("-t", "--thread_num", dest="thread_num")
    parser.add_option("-L", "--max_read_length", dest="max_read_length")
    parser.add_option("-f", "--bloom_filter", dest="bloom_filter")
    parser.add_option("-s", "--bloom_array_size", dest="bloom_array_size")
    parser.add_option("-b", "--num_processing_steps", dest="num_processing_steps")

    #Outputs
    parser.add_option("", "--stat", dest='stat', help="Statistical information")
    parser.add_option("", "--gz_freq", dest='gz_freq', help="Zipped frequency information")
    parser.add_option("", "--filelist", dest='filelist', help="List of files processed by KmerFreq")
    opts, args = parser.parse_args()

    #Create temp directory for performing analysis
    tmp_dir = tempfile.mkdtemp()
    #Create temp file for configuration
    tmp_file = tempfile.NamedTemporaryFile(dir=tmp_dir).name
    try:
        fout = open(tmp_file,'wb')
        read_index = 0
#        print opts.format_of_data
        for index in range(len(opts.paired_fastq_input1_list)):
            if opts.format_of_data[index]  == "fastq":
                #Loop thru files
                path = opts.paired_fastq_input1_list[index]
                fout.write(path)
                fout.write("\n")
                path = opts.paired_fastq_input2_list[index]
                fout.write(path)
                fout.write("\n")
            elif opts.format_of_data  == "fasta":
                #Loop thru files
                path = opts.paired_fasta_input1_list[index]
                fout.write(path)
                fout.write("\n")
                path = opts.paired_fasta_input2_list[index]
                fout.write(path)
                fout.write("\n")
        fout.close()
    except Exception, e:
        stop_err("Output config file cannot be opened for writing. " + str(e))

    #Files for std out and std error
    tmp_out_file = tempfile.NamedTemporaryFile(dir=tmp_dir).name
    tmp_stdout = open(tmp_out_file, 'wb')
    tmp_err_file = tempfile.NamedTemporaryFile(dir=tmp_dir).name
    tmp_stderr = open(tmp_err_file, 'wb')

    #Set up command line call - need to remove hard coded path
    if opts.default_full_settings_type == "default":
        cmd = "KmerFreq_HA_v2.0 -l %s >%s 2>%s" % (tmp_file, tmp_out_file, tmp_err_file)
    elif opts.default_full_settings_type == "full":
        cmd = "KmerFreq_HA_v2.0 -l %s -k %s" % (tmp_file, opts.kmer_size)

        if opts.read_length != "":
            cmd =  cmd + " -r %s" % opts.read_length

        if opts.use_num_bases != "":
            cmd =  cmd + " -g %s" % opts.use_num_bases

        cmd = cmd + " -a %s -d %s -i %s -t %s -L %s -f %s -s %s -b %s >%s 2>%s" % (opts.ignore_first, opts.ignore_last, opts.hash_size, opts.thread_num, opts.max_read_length, opts.bloom_filter, opts.bloom_array_size, opts.num_processing_steps, tmp_out_file, tmp_err_file)

    print "Command executed: ", cmd

    #Perform kmer analysis
    buffsize = 1048576
    #Temp directory to perform processing
    dirpath = tempfile.mkdtemp()
    print "Working directory: ", dirpath

    try:
        proc = subprocess.Popen(args=cmd, shell=True, cwd=dirpath, stdout=tmp_stdout, stderr=tmp_stderr)
        returncode = proc.wait()
        #Get stdout, allowing for case where it's very large
        tmp_stdout = open(tmp_out_file, 'rb')
        stdout = ''
        try:
            while True:
                stdout += tmp_stdout.read(buffsize)
                if not stdout or len(stdout) % buffsize != 0:
                    break
        except OverflowError:
            pass

        print stdout

        #Get stderr, allowing for case where it's very large
        tmp_stderr = open(tmp_err_file, 'rb')
        stderr = ''
        try:
            while True:
                stderr += tmp_stderr.read(buffsize)
                if not stderr or len(stderr) % buffsize != 0:
                    break
        except OverflowError:
            pass

        #Close streams
        tmp_stdout.close()
        tmp_stderr.close()
        if returncode != 0:
            raise Exception, stderr

    except Exception, e:
        raise Exception, 'Problem performing KmerFreq process ' + str(e)

    #Read KmerFreq results into outputs
    stat_out = open(opts.stat, 'wb')
    file = open(dirpath + '/output.freq.stat')
    for line in file:
        stat_out.write(line)
    stat_out.close()
    file.close()

    freq_out = open(opts.gz_freq, 'wb')
    with open(dirpath + "/output.freq.gz", mode='rb') as file: # b is important -> binary
        fileContent = file.read()
        freq_out.write(fileContent)
    freq_out.close()
    file.close()

    filelist_out = open(opts.filelist, 'w')
    read_list_handle = open(tmp_file, 'r')
    paths = read_list_handle.read()
    filelist_out.write(paths)
    filelist_out.close()
    read_list_handle.close()

    #Clean up temp files
    cleanup_before_exit(dirpath)
    #Check results in output file
    if os.path.getsize(opts.stat) > 0:
        sys.stdout.write('Status complete')
    else:
        stop_err("The output is empty")

if __name__ == "__main__": main()
