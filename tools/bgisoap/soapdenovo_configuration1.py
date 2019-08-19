"""
soapdenovo_configuration.py
create configuration file for SOAPdenovo modules
Copyright   Huayan Gao - CUHK and CBIIT
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
        help="Length of contig to be aligned for a reliable read location")

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
    parser.add_option("-a", "--init_mem_assumption", dest="init_mem_assumption", help="Initialize memory assumption to avoid further reallocation, unit G")
    parser.add_option("-d", "--kmer_freq_cutoff", dest="kmer_freq_cutoff", help="kmers with frequency no larger than KmerFreqCutoff will be deleted")
    parser.add_option("-R", "--output_extra_info", dest="output_extra_info", help="Output extra information for resolving repeats in contig step")

    #Outputs
    parser.add_option("", "--soap_config", dest='soap_config')

    opts, args = parser.parse_args()

    #Need a temporary directory to perform processing
    dirpath = tempfile.mkdtemp()

    #Create temp file to store soapdenovo2 running configuration
    config_file = tempfile.NamedTemporaryFile(dir=dirpath, prefix="soap_",suffix=".config").name

    try:
        fout = open(config_file,'w')
        fout.write("max_rd_len=%s\n" % opts.max_read_length)
        #Calculate how many sets of data there are - use avg_ins as a measure of this
        #Separate indices required to keep count of reads
        single_read_index = 0
        paired_read_index = 0
        for index in range(len(opts.avg_insert_list)):
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
        stop_err("File cannot be opened for writing soap.config " + str(e))

    config_out = open(opts.soap_config, 'wb')
    file = open(config_file)
    for line in file:
        config_out.write(line)
    config_out.close()
    file.close()

    #Clean up temp files
    #cleanup_before_exit(tmp_dir)
    #Check results in output file
    if os.path.getsize(opts.soap_config) > 0:
        sys.stdout.write('Status complete')
    else:
        stop_err("The output is empty")

if __name__ == "__main__": main()
