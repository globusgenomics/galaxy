"""
soapdenovo2_map.py
A wrapper script for SOAPdenovo2 map module
Copyright   Peter Li - GigaScience and BGI-HK
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
    #Inputs
    parser.add_option('', '--contig', dest='contig')
    parser.add_option('', '--contigindex', dest='contigindex')
    parser.add_option('', '--soap_config', dest='soap_config')

    parser.add_option("", "--analysis_settings_type", dest="analysis_settings_type")
    parser.add_option("", "--default_full_settings_type", dest="default_full_settings_type")

    parser.add_option("-p", "--ncpu", dest="ncpu", help="Number of cpu for use")
    parser.add_option("-f", "--output_gap_related_reads", dest="output_gap_related_reads")
    parser.add_option("-k", "--kmer_r2c", dest="kmer_r2c")

    #Outputs
    parser.add_option("", "--pegrads", dest='pegrads')
    parser.add_option("", "--read_on_contig", dest='read_on_contig')
    parser.add_option("", "--read_in_gap", dest='read_in_gap')

    opts, args = parser.parse_args()

    #Need to write inputs to files in a temporary directory for use by executable
    dirpath = tempfile.mkdtemp()

    contig_data = open(opts.contig, 'r')
    contig_file = open(dirpath + "/out.contig", "w")
    for line in contig_data:
        contig_file.write(line)
    contig_data.close()
    contig_file.close()

    contigindex_data = open(opts.contigindex, 'r')
    contigindex_file = open(dirpath + "/out.ContigIndex", "w")
    for line in contigindex_data:
        contigindex_file.write(line)
    contigindex_data.close()
    contigindex_file.close()

    #Set up command line call
    #TODO - remove hard coded path; shouldn't need the full path to the executable as SOAPdenovo-63mer is declared on the path
    #TODO - Need to declare path to Galaxy
    #Code for adding directory path to other file required as output
    if opts.default_full_settings_type == "default":
        cmd = "SOAPdenovo-63mer_v2.0 map -s %s -g %s" % (opts.soap_config, dirpath + "/out")
    elif opts.default_full_settings_type == "full":
        cmd = "SOAPdenovo-63mer_v2.0 map -s %s -g %s -f %s -p % -k %s" % (opts.soap_config, dirpath + "/out", opts.output_gap_related_reads, opts.ncpu, opts.kmer_r2c)

    #Check
    print cmd

    #Perform SOAPdenovo2_map analysis
    buffsize = 1048576
    #Create temp directory for standard error and out
    tmp_dir = tempfile.mkdtemp()
    try:
        tmp_out_file = tempfile.NamedTemporaryFile(dir=tmp_dir).name
        tmp_stdout = open(tmp_out_file, 'wb')
        tmp_err_file = tempfile.NamedTemporaryFile(dir=tmp_dir).name
        tmp_stderr = open(tmp_err_file, 'wb')

        #Call SOAPdenovo2
        proc = subprocess.Popen(args=cmd, shell=True, cwd=dirpath, stdout=tmp_stdout, stderr=tmp_stderr.fileno())
        returncode = proc.wait()
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
        raise Exception, 'Problem performing map process ' + str(e)

    pegrads_out = open(opts.pegrads, 'w')
    file = open(dirpath + "/out.peGrads")
    for line in file:
        pegrads_out.write(line)
    pegrads_out.close()
    file.close()

    read_on_contig_out = open(opts.read_on_contig, 'wb')
    with open(dirpath + "/out.readOnContig.gz", mode='rb') as file: # b is important -> binary
        fileContent = file.read()
        read_on_contig_out.write(fileContent)
    read_on_contig_out.close()
    file.close()

    read_in_gap_out = open(opts.read_in_gap, 'wb')
    with open(dirpath + "/out.readInGap.gz", mode='rb') as file: # b is important -> binary
        fileContent = file.read()
        read_in_gap_out.write(fileContent)
    read_in_gap_out.close()
    file.close()

    #Clean up temp files
    cleanup_before_exit(tmp_dir)
    #Check results in output file
    if os.path.getsize(opts.pegrads) > 0:
        sys.stdout.write('Status complete')
    else:
        stop_err("The output is empty")

if __name__ == "__main__": main()
