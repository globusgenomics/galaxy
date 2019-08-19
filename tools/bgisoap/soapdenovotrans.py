"""
soapdenovotrans.py
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
    parser.add_option("", "--avg_ins", dest="avg_ins", help="Average insert size")
    parser.add_option("", "--reverse_seq", dest="reverse_seq", help="Reverse sequence?")
    parser.add_option("", "--asm_flags", dest="asm_flags", help="Which operations should the reads be used for?")
    parser.add_option("", "--rank", dest="rank", help="Which order are the reads used while scaffolding")
    #Data inputs
    parser.add_option("", "--type_of_data", dest="type_of_data")
    parser.add_option("", "--format_of_data", dest="format_of_data")
    parser.add_option("", "--single_fastq_input1", dest="single_fastq_input1")
    parser.add_option("", "--single_fastq_input2", dest="single_fastq_input2")
    parser.add_option("", "--single_fasta_input1", dest="single_fasta_input1")
    parser.add_option("", "--single_fasta_input2", dest="single_fasta_input2")
    #TODO: Remaining SOAPdenovo-trans parameters need to be implemented
    #Outputs
    parser.add_option("", "--contig", dest='contig', help="Contig sequence file")
    parser.add_option("", "--scafseq", dest='scafseq', help="Scaffold sequence file")
    opts, args = parser.parse_args()

    #Create temp directory for performing analysis
    tmp_dir = tempfile.mkdtemp()
    #Create temp file for configuration
    tmp_file = tempfile.NamedTemporaryFile(dir=tmp_dir).name
    try:
        fout = open(tmp_file,'wb')
        fout.write("max_rd_len=%s\n" % opts.max_read_length)
        fout.write("[LIB]\n")
        fout.write("avg_ins=%s\n" % opts.avg_ins)
        fout.write("reverse_seq=%s\n" % opts.reverse_seq)
        fout.write("asm_flags=%s\n" % opts.asm_flags)
        fout.write("rank=%s\n" % opts.rank)
        #Add data file configuration
        if opts.type_of_data == "single":
            if opts.format_of_data == "fastq":
                fout.write("q=%s\n" % opts.single_fastq_input1)
            elif opts.format_of_data == "fasta":
                fout.write("f=%s\n" % opts.single_fasta_input1)
        elif opts.type_of_data == "paired":
            if opts.format_of_data == "fastq":
                fout.write("q1=%s\n" % opts.single_fastq_input1)
                fout.write("q2=%s\n" % opts.single_fastq_input2)
            elif opts.format_of_data == "fasta":
                fout.write("f1=%s\n" % opts.single_fasta_input1)
                fout.write("f2=%s\n" % opts.single_fasta_input2)
        fout.close()
    except Exception, e:
        stop_err("Output config file cannot be opened for writing." + str(e))

    #Set up command line call - need to remove hard coded path
    cmd = "SOAPdenovo-Trans-31kmer all -s %s -o %s" % (tmp_file, tmp_dir + "/result")
#    cmd = "SOAPdenovo-Trans-31kmer all -s %s -o %s" % (tmp_file, tmp_dir + "/result")

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
        raise Exception, 'Problem performing transcriptome sequencing. ' + str(e)

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
        sys.stdout.write('Status complete')
    else:
        stop_err("The output is empty")

if __name__ == "__main__": main()
