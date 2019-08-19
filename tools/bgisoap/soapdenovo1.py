"""
soapdenovo.py
A wrapper script for SOAPdenovo-trans
Copyright   Peter Li - GigaScience and BGI-HK
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
    parser.add_option("", "--rank", action="append", type="string", dest="rank_list", help="Which order are the reads used while scaffolding")
    #Data inputs
    parser.add_option("", "--type_of_data", action="append", type="string", dest="type_of_data_list")
    parser.add_option("", "--format_of_data", action="append", type="string", dest="format_of_data_list")
    parser.add_option("", "--single_fastq_input1", action="append", type="string", dest="single_fastq_input1_list")
    parser.add_option("", "--single_fasta_input1", action="append", type="string", dest="single_fasta_input1_list")

    parser.add_option("", "--paired_fastq_input1", action="append", type="string", dest="paired_fastq_input1_list")
    parser.add_option("", "--paired_fastq_input2", action="append", type="string", dest="paired_fastq_input2_list")
    parser.add_option("", "--paired_fasta_input1", action="append", type="string", dest="paired_fasta_input1_list")
    parser.add_option("", "--paired_fasta_input2", action="append", type="string", dest="paired_fasta_input2_list")
    #TODO: Remaining SOAPdenovo parameters need to be implemented
    #Outputs
    parser.add_option("", "--contig", dest='contig', help="Contig sequence file")
    parser.add_option("", "--scafseq", dest='scafseq', help="Scaffold sequence file")
    opts, args = parser.parse_args()

    #Create temp directory for performing analysis
    tmp_dir = tempfile.mkdtemp()
    #Create temp file for configuration
    
    if opts.file_source == "history":
        config_file = opts.configuration
    else:
        #Create temp file for configuration
        config_file = tempfile.NamedTemporaryFile(dir=tmp_dir).name
        try:
            fout = open(config_file,'w')
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
		fout.write("rank=%s\n" % opts.rank_list[index])
		#Add data file configuration - needs careful looping due to single and paired reads
		print len(opts.type_of_data_list)
		print opts.type_of_data_list[index]
		print opts.format_of_data_list[index]
		if opts.type_of_data_list[index] == "single":  #then only one read
		    if (opts.format_of_data_list)[index] == "fastq":
		        fout.write("q=%s\n" % (opts.single_fastq_input1_list)[single_read_index])
		    elif opts.format_of_data == "fasta":
		        fout.write("f=%s\n" % opts.single_fasta_input1_list[single_read_index])
		    single_read_index = single_read_index + 1
		elif opts.type_of_data_list[index] == "paired":
		    if opts.format_of_data_list[index] == "fastq":
		        fout.write("q1=%s\n" % (opts.paired_fastq_input1_list)[paired_read_index])
		        fout.write("q2=%s\n" % (opts.paired_fastq_input2_list)[paired_read_index])
		    elif opts.format_of_data_list[index] == "fasta":
		        fout.write("f1=%s\n" % opts.paired_fasta_input1_list[paired_read_index])
		        fout.write("f2=%s\n" % opts.paired_fasta_input2_list[paired_read_index])
		    paired_read_index = paired_read_index + 1
	    fout.close()
	except Exception, e:
	    stop_err("Output configuration file cannot be opened for writing." + str(e))


    #Set up command line call - need to remove hard coded path
#    cmd = "/usr/local/bgisoap/soapdenovo/current/bin/SOAPdenovo-63mer all -s %s -o %s" % (tmp_file, tmp_dir + "/result")
    cmd = "SOAPdenovo-63mer all -s %s -o %s" % (config_file, tmp_dir + "/result")
    print cmd

    #Perform SOAPdenovo-trans analysis
    buffsize = 1048576
    try:
        #Create file in temporary directory
        tmp = tempfile.NamedTemporaryFile(dir=tmp_dir).name
        #Open a stream to file
        tmp_stderr = open(tmp, 'wb')
        #Call SOAPdenovo-trans
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
        raise Exception, 'Problem performing de novo assembly ' + str(e)

    #Read SOAPdenovo-trans results into outputs
    
    contig_out = open(opts.contig, 'wb')
    file = open(tmp_dir + '/result.contig')
    for line in file:
        #print line
        contig_out.write(line)
    contig_out.close()

    scafseq_out = open(opts.scafseq, 'wb')
    file = open(tmp_dir + '/result.scafSeq')
    for line in file:
        #print line
        scafseq_out.write(line)
    scafseq_out.close()

    #Clean up temp files
    cleanup_before_exit(tmp_dir)
    #Check results in output file
    if os.path.getsize(opts.contig) > 0:
        sys.stdout.write('Status complete')
    else:
        stop_err("The output is empty")

if __name__ == "__main__": main()
