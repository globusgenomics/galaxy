"""
soapdenovo2.py
A wrapper script for SOAPdenovo2
Peter Li - GigaScience and BGI-HK
Huayan Gao - CUHK
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
    parser.add_option("", "--max_read_length", dest="max_read_length", help="Maximum read length")
    parser.add_option("", "--file_source", dest="file_source", help="Select a configuration file from history or create a new one?")
    parser.add_option("", "--configuration", dest="configuration", help="Select a configuration file from history or create a new one?")

    #Make list of params
    parser.add_option("", "--avg_ins", action="append", type="string", dest="avg_insert_list", help="Average insert size")
    parser.add_option("", "--reverse_seq", action="append", type="string", dest="reverse_seq_list", help="Reverse sequence?")
    parser.add_option("", "--asm_flags", action="append", type="string", dest="asm_flags_list", help="Which operations should the reads be used for?")
    parser.add_option("", "--rd_len_cutoff", action="append", type="string", dest="rd_len_cutoff_list",
        help="Number of base pairs to use from reads")
    parser.add_option("", "--rank", action="append", type="string", dest="rank_list", help="Which order are the reads used while scaffolding")
    parser.add_option("", "--pair_num_cutoff", action="append", type="string", dest="pair_num_cutoff_list",
        help="Pair number cutoff for a reliable connection")
    parser.add_option("", "--map_len", action="append", type="string", dest="map_len_list",
        help="Length of contig that has to be aligned for a reliable read location")

    #Data inputs
    parser.add_option("", "--type_of_data", action="append", type="string", dest="type_of_data_list")
    parser.add_option("", "--format_of_data", action="append", type="string", dest="format_of_data_list")
    parser.add_option("", "--single_fastq_input1", action="append", type="string", dest="single_fastq_input1_list")
    parser.add_option("", "--single_fasta_input1", action="append", type="string", dest="single_fasta_input1_list")
    parser.add_option("", "--single_bam_input1", action="append", type="string", dest="single_bam_input1_list")

    parser.add_option("", "--paired_fastq_input1", action="append", type="string", dest="paired_fastq_input1_list")
    parser.add_option("", "--paired_fastq_input2", action="append", type="string", dest="paired_fastq_input2_list")
    parser.add_option("", "--paired_fasta_input1", action="append", type="string", dest="paired_fasta_input1_list")
    parser.add_option("", "--paired_fasta_input2", action="append", type="string", dest="paired_fasta_input2_list")
    parser.add_option("", "--paired_bam_input1", action="append", type="string", dest="paired_bam_input1_list")
    parser.add_option("", "--paired_bam_input2", action="append", type="string", dest="paired_bam_input2_list")

    parser.add_option("", "--analysis_settings_type", dest="analysis_settings_type")
    parser.add_option("", "--default_full_settings_type", dest="default_full_settings_type")

    #Custom params
    parser.add_option("-K", "--kmer_size", dest="kmer_size", help="kmer size")
    parser.add_option("-p", "--ncpu", dest="ncpu", help="Number of cpu for use")
    parser.add_option("-a", "--init_memory_assumption", dest="init_memory_assumption", help="Memory assumption initialized to avoid further reallocation, unit G")
    parser.add_option("-d", "--kmer_freq_cutoff", dest="kmer_freq_cutoff", help="kmers with frequency no larger than KmerFreqCutoff will be deleted")
    parser.add_option("-R", "--resolve_repeats", dest="resolve_repeats", help="Resolve repeats by reads")
    parser.add_option("-D", "--edge_cov_cutoff", dest="edge_cov_cutoff", help="Edges with coverage no larger than EdgeCovCutoff will be deleted")
    parser.add_option("-M", "--merge_level", dest="merge_level", help="The strength of merging similar sequences during contiging")
    parser.add_option("-m", "--max_k", dest="max_k", help="max k when using multi kmer")
    parser.add_option("-e", "--weight", dest="weight", help="weight to filter arc when linearize two edges")
    #parser.add_option("-r", "--keep_avail_read", dest="keep_avail_read", help="keep available read")
    parser.add_option("-E", "--merge_clean_bubble", dest="merge_clean_bubble", help="Merge clean bubble before iterate")
    #parser.add_option("-f", "--output_gap_related_reads", dest="output_gap_related_reads", help="Output gap related reads in map step for using SRkgf to fill gap")
    parser.add_option("-k", "--kmer_r2c", dest="kmer_r2c", help="kmer size used for mapping read to contig")
    parser.add_option("-F", "--fill_gaps", dest="fill_gaps", help="Fill gaps in scaffold")
    parser.add_option("-u", "--unmask_contigs", dest="unmask_contigs", help="Unmask contigs with high/low coverage before scaffolding")
    parser.add_option("-w", "--keep_contigs_connected", dest="keep_contigs_connected",
        help="Keep contigs weakly connected to other contigs in scaffold")
    parser.add_option("-G", "--gap_len_diff", dest="gap_len_diff", help="Allowed length difference between estimated and filled gap")
    parser.add_option("-L", "--min_contig_len", dest="min_contig_len", help="Shortest contig for scaffolding")
    parser.add_option("-c", "--min_contig_cvg", dest="min_contig_cvg", help="minimum contig coverage (c*avgCvg) contigs shorter than 100bp with coverage smaller than c*avgCvg will be masked before scaffolding unless -u is set")
    parser.add_option("-C", "--max_contig_cvg", dest="max_contig_cvg", help="maximum contig coverage (C*avgCvg), contigs with coverage larger than C*avgCvg or contigs shorter than 100bp with coverage larger than 0.8*C*avgCvg will be masked before scaffolding unless -u is set")
    parser.add_option("-b", "--insert_size_upper_bound", dest="insert_size_upper_bound", help="will be used as upper bound of insert size for large insert size ( > 1000) when handling pair-end connections between contigs if b is set to larger than 1")
    parser.add_option("-B", "--bubble_coverage", dest="bubble_coverage", help="remove contig with lower cvoerage in bubble structure if both contigs' coverage are smaller than bubbleCoverage*avgCvg")
    parser.add_option("-N", "--genome_size", dest="genome_size", help="Genome size for statistics")
    parser.add_option("-V", "--ass_visual", dest="ass_visual", help="Output visualization information of assembly")

    #Outputs
    parser.add_option("", "--contig", dest='contig', help="Contig sequence file")
    parser.add_option("", "--scafseq", dest='scafseq', help="Scaffold sequence file")
    opts, args = parser.parse_args()

    #Create temp directory for performing analysis
    tmp_dir = tempfile.mkdtemp()

    if opts.file_source == "history":
        config_file = opts.configuration
    else:
        #Create temp file for configuration
        config_file = tempfile.NamedTemporaryFile(dir=tmp_dir).name
        try:
            fout = open(config_file,'wb')
            fout.write("max_rd_len=%s\n" % opts.max_read_length)
            #Calculate how sets of data there are - use avg_ins as a measure of this
            #Loop thru this number of times
            #Also use separate index to keep count of reads
            single_read_index = 0
            paired_read_index = 0
            for index in range(len(opts.avg_insert_list)):
                print "single_read_index ", single_read_index
                print "paired_read_index ", paired_read_index
                fout.write("[LIB]\n")
                fout.write("avg_ins=%s\n" % (opts.avg_insert_list)[index])
                fout.write("reverse_seq=%s\n" % opts.reverse_seq_list[index])
                fout.write("asm_flags=%s\n" % opts.asm_flags_list[index])
                fout.write("rd_len_cutoff=%s\n" % opts.rd_len_cutoff_list[index])
                fout.write("rank=%s\n" % opts.rank_list[index])
                fout.write("pair_num_cutoff=%s\n" % opts.pair_num_cutoff_list[index])
                fout.write("map_len=%s\n" % opts.map_len_list[index])
                #Add data file configuration - needs careful looping due to single and paired reads
                print opts.type_of_data_list[index]
                print opts.format_of_data_list[index]
                if opts.type_of_data_list[index] == "single":  #then only one read
                    if (opts.format_of_data_list)[index] == "fastq":
                        fout.write("q=%s\n" % (opts.single_fastq_input1_list)[single_read_index])
                    elif opts.format_of_data == "fasta":
                        fout.write("f=%s\n" % opts.single_fasta_input1_list[single_read_index])
                    else:
                        fout.write("b=%s\n" % opts.single_bam_input1_list[single_read_index])
                    single_read_index = single_read_index + 1
                elif opts.type_of_data_list[index] == "paired":
                    if opts.format_of_data_list[index] == "fastq":
                        fout.write("q1=%s\n" % (opts.paired_fastq_input1_list)[paired_read_index])
                        fout.write("q2=%s\n" % (opts.paired_fastq_input2_list)[paired_read_index])
                    elif opts.format_of_data_list[index] == "fasta":
                        fout.write("f1=%s\n" % opts.paired_fasta_input1_list[paired_read_index])
                        fout.write("f2=%s\n" % opts.paired_fasta_input2_list[paired_read_index])
                    else:
                        fout.write("b1=%s\n" % opts.paired_fasta_input1_list[paired_read_index])
                        fout.write("b2=%s\n" % opts.paired_fasta_input2_list[paired_read_index])
                    paired_read_index = paired_read_index + 1
            fout.close()
        except Exception, e:
            stop_err("config file cannot be opened for writing" + str(e))

    #Set up command line call - assumes path to executable has been defined in users environment
    #Need to check kmer size to decide if to use SOAPdenovo-63mer or SOAPdenovo-127mer
    if opts.kmer_size <= 63 and opts.default_full_settings_type == "default":
        cmd = "SOAPdenovo-63mer all -s %s -o %s" % (config_file, tmp_dir + "/result")
    elif opts.kmer_size <= 63 and opts.default_full_settings_type == "full":
        cmd = "SOAPdenovo-63mer all -s %s -o %s -K %s -p %s -a %s -d %s -R %s -D %s -M %s -m %s -e %s -E %s -k %s -F %s -u %s -w %s -G %s -L %s -c %s -C %s -b %s -B %s -N %s -V %s" % (config_file, tmp_dir + "/result", opts.kmer_size, opts.ncpu, opts.init_memory_assumption, opts.kmer_freq_cutoff, opts.resolve_repeats, opts.edge_cov_cutoff, opts.merge_level, opts.max_k, opts.weight, opts.merge_clean_bubble, opts.kmer_r2c, opts.fill_gaps, opts.unmask_contigs, opts.keep_contigs_connected, opts.gap_len_diff, opts.min_contig_len, opts.min_contig_cvg, opts.max_contig_cvg, opts.insert_size_upper_bound, opts.bubble_coverage, opts.genome_size, opts.ass_visual)
    elif opts.kmer_size >= 64 and opts.default_full_settings_type == "default":
        cmd = "SOAPdenovo-127mer all -s %s -o %s" % (config_file, tmp_dir + "/result")
    else:
        cmd = "SOAPdenovo-127mer all -s %s -o %s -K %s -p %s -a %s -d %s -R %s -D %s -M %s -m %s -e %s -E %s -k %s -F %s -u %s -w %s -G %s -L %s -c %s -C %s -b %s -B %s -N %s -V %s" % (config_file, tmp_dir + "/result", opts.kmer_size, opts.ncpu, opts.init_memory_assumption, opts.kmer_freq_cutoff, opts.resolve_repeats, opts.edge_cov_cutoff, opts.merge_level, opts.max_k, opts.weight, opts.merge_clean_bubble, opts.kmer_r2c, opts.fill_gaps, opts.unmask_contigs, opts.keep_contigs_connected, opts.gap_len_diff, opts.min_contig_len, opts.min_contig_cvg, opts.max_contig_cvg, opts.insert_size_upper_bound, opts.bubble_coverage, opts.genome_size, opts.ass_visual)

    print cmd

    #Perform SOAPdenovo2 analysis
    buffsize = 1048576
    try:
        #Create file in temporary directory
        tmp = tempfile.NamedTemporaryFile(dir=tmp_dir).name
        #Open a stream to file
        tmp_stderr = open(tmp, 'wb')
        #Call SOAPdenovo
        proc = subprocess.Popen(args=cmd, shell=True, cwd=tmp_dir, stderr=tmp_stderr.fileno())
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
        raise Exception, 'Problem performing SOAPdenovo process ' + str(e)

    #Read SOAPdenovo-trans results into outputs
    contig_out = open(opts.contig, 'wb')
    file = open(tmp_dir + '/result.contig')
    for line in file:
        print line
        contig_out.write(line)
    contig_out.close()

    scafseq_out = open(opts.scafseq, 'wb')
    file = open(tmp_dir + '/result.scafSeq')
    for line in file:
        print line
        scafseq_out.write(line)
    scafseq_out.close()

    #Clean up temp files
    cleanup_before_exit(tmp_dir)

    #Check results in output file
    if os.path.getsize(opts.contig) > 0:
        sys.stdout.write('SOAPdenovo process finished')
    else:
        stop_err("Problem with SOAPdenovo process")

if __name__ == "__main__": main()
