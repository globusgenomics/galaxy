"""
snvdetect.py
A wrapper script for snvdetect
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

    ## Inputs
    parser.add_option("", "--chr_length", dest="chr_length", help="chromosome length info file")
    parser.add_option("", "--input1_alleleqc", dest="input1_alleleqc", help="First input alleleQC file")
    parser.add_option("", "--input2_alleleqc", dest="input2_alleleqc", help="Second input alleleQC file")
    parser.add_option("", "--snvtest", dest="snvtest", help="SNVtest file")

    parser.add_option("", "--default_full_settings_type", dest="default_full_settings_type")
    #Custom params
    parser.add_option("-m", "--min_depth", dest="min_depth", help="Minimum depth requirement")
    parser.add_option("-x", "--max_depth", dest="max_depth", help="Maximum depth requirement")
    parser.add_option("-c", "--max_copy_num", dest="max_copy_num", help="Maximum copy number")
    parser.add_option("-r", "--pvalue_rank_sum", dest="pvalue_rank_sum", help="p-value cutoff for quality rank sum test")
    parser.add_option("-u", "--pvalue_mismatch_enrichment_repetitive_hits", dest="pvalue_mismatch_enrichment_repetitive_hits", help="p-value cutoff for mismatch enrichment test on repetitive hits")
    parser.add_option("-e", "--pvalue_mismatch_enrichment_read_ends", dest="pvalue_mismatch_enrichment_read_ends", help="p-value cutoff for mismatch enrichment test on 5' or 3' end of reads")
    parser.add_option("-l", "--relax_depth_requirement", dest="relax_depth_requirement", help="Relax the depth requirement. If this parameter was not given, sites which unique aligned reads less than minimum depth requirement would be filtered.")
    parser.add_option("-i", "--min_allele_freq", dest="min_allele_freq", help="Minimum allele frequency requirement")
    parser.add_option("-f", "--allele_freq_requirement", dest="allele_freq_requirement", help="Allele frequency requirement for confident genotype call")
    parser.add_option("-t", "--max_allele_freq", dest="max_allele_freq", help="Maximum allowed mutant allele frequency in normal for somatic mutation detection")
    parser.add_option("-s", "--pvalue_somatic_mut_detection", dest="pvalue_somatic_mut_detection", help="p-value cutoff for somatic mutation detection")

    #Outputs
    parser.add_option("-o", "--output", dest='output', help="Results from SNVdetect process")
    opts, args = parser.parse_args()

    #Create temp directory
    tmp_dir = tempfile.mkdtemp()

    #tmp_output = tempfile.NamedTemporaryFile(dir=tmp_dir).name+".gz"
    
    #Set up command line call
    if opts.default_full_settings_type == "default":
        cmd = "SNVdetect %s %s %s %s %s" % (opts.chr_length, opts.input1_alleleqc, opts.input2_alleleqc, opts.snvtest, opts.output)
    elif opts.default_full_settings_type == "full":
        cmd = "SNVdetect %s %s %s %s %s -m %s -x %s -c %s -r %s -u %s -e %s -l %s -i %s -f %s -t %s -s %s" % (opts.chr_length, opts.input1_alleleqc, opts.input2_alleleqc, opts.snvtest, opts.output, opts.min_depth, opts.max_depth, opts.max_copy_num, opts.pvalue_rank_sum, opts.pvalue_mismatch_enrichment_repetitive_hits, opts.pvalue_mismatch_enrichment_read_ends, opts.relax_depth_requirement, opts.min_allele_freq, opts.allele_freq_requirement, opts.max_allele_freq, opts.pvalue_somatic_mut_detection)

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
        stop_err('Error in running SNVtest from (%s), %s' % (opts.output, str(e)))

    #pass the tmp_output to opts.output
#    gz_output = open(opts.output, 'wb')
#    with open(tmp_output, mode='rb') as file: # b is important -> binary
#        fileContent = file.read()
#        gz_output.write(fileContent)
#    gz_output.close()
#    file.close()
    
    #Clean up temp files
    cleanup_before_exit(tmp_dir)
    #Check results in output file
    if os.path.getsize(opts.output) > 0:
        sys.stdout.write('Status complete')
    else:
        stop_err("The output is empty")

if __name__ == "__main__": __main__()
