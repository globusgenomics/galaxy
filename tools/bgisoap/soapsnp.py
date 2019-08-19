"""
soapsnp.py
A wrapper script for SOAPsnp
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
    parser.add_option('-i', '--soap_alignment', dest='soap_alignment', help='A sorted SOAP alignment result')
    parser.add_option('-d', '--ref_seq', dest='ref_seq', help='A reference DNA sequence in FASTA format')
    parser.add_option('', '--genome_type', dest='genome_type',
        help='Whether performing diploid or monoploid re-sequencing')

    #Custom params
    parser.add_option("", "--default_full_settings_type", dest="default_full_settings_type")
    parser.add_option('', '--quality_calibration_matrix_setting', dest='quality_calibration_matrix_setting',
        help='Check if user wants to include a quality calibration matrix')
    parser.add_option('', '--output_calibration_matrix_setting', dest='output_calibration_matrix_setting',
        help='Check if user wants to output a quality calibration matrix')

    parser.add_option("", "--include_snp_info", dest="include_snp_info",
        help="Check if user wants to include pre-formatted known SNP information in analysis")
    parser.add_option("-s", "--snp_info", dest="snp_info", help="Pre-formatted known SNP information")

    parser.add_option('', '--call_consensus_setting', dest='call_consensus_setting',
        help='Check if user wants to call consensus on specific regions')
    parser.add_option('-T', '--specific_regions', dest='specific_regions',
        help='File containing specific chromosome regions')

    parser.add_option("-z", "--quality_score_char", dest="quality_score_char",
        help="ASCII character that stands for quality score==0")
    parser.add_option("-g", "--global_error_dependency_coefficient", dest="global_error_dependency_coefficient",
        help="Global error dependency coefficient, 0.0(complete dependent)~1.0(complete independent)")
    parser.add_option("-p", "--pcr_error_dependency_coefficient", dest="pcr_error_dependency_coefficient",
        help="PCR error dependency coefficient")
    parser.add_option("-r", "--novel_althom_prior_probability", dest="novel_althom_prior_probability",
        help="Novel altHOM prior probability")
    parser.add_option("-e", "--novel_het_prior_probability", dest="novel_het_prior_probability",
        help="Novel HET prior probability")
    parser.add_option("-t", "--ratio", dest="ratio",
        help="Set transition/transversion ratio to 2:1 in prior probability")

    parser.add_option("-2", "--refine_snp_calling", dest="refine_snp_calling",
        help="Specifying this option will refine SNP calling using known SNP information")
    parser.add_option("-a", "--validated_het_prior", dest="validated_het_prior",
        help="Validated HET prior, if no allele frequency known")
    parser.add_option("-b", "--validated_althom_prior", dest="validated_althom_prior",
        help="Validated altHOM prior, if no allele frequency known")
    parser.add_option("-j", "--unvalidated_het_prior", dest="unvalidated_het_prior",
        help="Unvalidated HET prior, if no allele frequency known")
    parser.add_option("-k", "--unvalidated_althom_rate", dest="unvalidated_althom_rate",
        help="Unvalidated altHOM rate, if no allele frequency known")
    parser.add_option("-u", "--enable_rank_sum", dest="enable_rank_sum",
        help="Enable rank sum test that check whether the two allele of a possible HET call have same sequencing quality to give HET further penalty for better accuracy")
    parser.add_option("-n", "--enable_binom_calc", dest="enable_binom_calc",
        help="Enable binomial probability calculation (that check whether the two allele are observed equally)to give HET further penalty for better accuracy")
    parser.add_option("-m", "--enable_monoploid_calling", dest="enable_monoploid_calling",
        help="Enable monoploid calling mode, this will ensure all consensus as HOM and you probably should SPECIFY higher altHOM rate")
    parser.add_option("-q", "--output_potential_snps", dest="output_potential_snps",
        help="Only output potential SNPs. Useful in Text output mode")
    parser.add_option("-L", "--max_length_short_read", dest="max_length_short_read",
        help="Maximum length of read. Please note that once length of some reads exceeds the parameter will probably collapse the program.")
    parser.add_option("-Q", "--max_fastq_score", dest="max_fastq_score", help="Maximum FASTQ quality score")
    parser.add_option("-F", "--output_format", dest="output_format", help="Output format")

    #Outputs
    parser.add_option("-o", "--consensus_out", dest='consensus_out',
        help="An alignment of paired reads on a reference sequence")
    #If this is set by the user
    parser.add_option("-M", "--quality_calibration_matrix_out", dest='quality_calibration_matrix_out',
        help="Quality calibration matrix")
    opts, args = parser.parse_args()

    #Create temp directory
    tmp_dir = tempfile.mkdtemp()
    print tmp_dir

    #Get reference index directory
    genome_type = opts.genome_type
    print genome_type

    #Set up command line call
    #Need to fix to account for user configuration at run time
    #12/12/2012, Huayan
    # use the default 45 for opts.max_length_short_read since it is not set up
    cmd="test only. empty output."
    print opts.default_full_settings_type
    print opts.output_calibration_matrix_setting
    print cmd
    if opts.default_full_settings_type == "default":
        cmd = "soapsnp -i %s -d %s -o %s -r %s -e %s -t %s -u %s -L 45 -M %s -m %s" % (
        opts.soap_alignment, opts.ref_seq, opts.consensus_out, opts.novel_althom_prior_probability,
        opts.novel_het_prior_probability, opts.ratio, opts.enable_rank_sum,
        opts.quality_calibration_matrix_out, opts.enable_monoploid_calling)
    elif opts.default_full_settings_type == "full" and opts.output_calibration_matrix_setting == "yes":
        cmd = "soapsnp -i %s -d %s -o %s -z %s -g %s -p %s -r %s -e %s -t %s -s %s -2 %s -a %s -b %s -j %s -k %s -u %s -n %s -m %s -q %s -M %s -L %s -Q %s -F %s" % (
        opts.soap_alignment, opts.ref_seq, opts.consensus_out, opts.quality_score_char, opts.global_error_dependency_coefficient,
        opts.pcr_error_dependency_coefficient, opts.novel_althom_prior_probability, opts.novel_het_prior_probability,
        opts.ratio, opts.snp_info, opts.refine_snp_calling, opts.validated_het_prior, opts.validated_althom_prior,
        opts.unvalidated_het_prior, opts.unvalidated_althom_rate, opts.enable_rank_sum, opts.enable_binom_calc,
        opts.enable_monoploid_calling, opts.output_potential_snps, opts.quality_calibration_matrix_out,
        opts.max_length_short_read, opts.max_fastq_score, opts.output_format)
    else:
        cmd = "soapsnp -i %s -d %s -o %s -z %s -g %s -p %s -r %s -e %s -t %s -s %s -2 %s -a %s -b %s -j %s -k %s -u %s -n %s -m %s -q %s -L %s -Q %s -F %s" % (
        opts.soap_alignment, opts.ref_seq, opts.consensus_out, opts.quality_score_char, opts.global_error_dependency_coefficient,
        opts.pcr_error_dependency_coefficient, opts.novel_althom_prior_probability, opts.novel_het_prior_probability,
        opts.ratio, opts.snp_info, opts.refine_snp_calling, opts.validated_het_prior, opts.validated_althom_prior,
        opts.unvalidated_het_prior, opts.unvalidated_althom_rate, opts.enable_rank_sum, opts.enable_binom_calc,
        opts.enable_monoploid_calling, opts.output_potential_snps, 
        opts.max_length_short_read, opts.max_fastq_score, opts.output_format)

    print cmd

    #Need to check format of reference sequence
    #Has to be in fasta format

    #Run
    try:
        tmp_out_file = tempfile.NamedTemporaryFile(dir=tmp_dir).name
        tmp_stdout = open(tmp_out_file, 'wb')
        tmp_err_file = tempfile.NamedTemporaryFile().name
        tmp_stderr = open(tmp_err_file, 'wb')

        #Perform SOAPsnp call
        print "Doing SOAPsnp call..."
        soap_proc = subprocess.Popen(args=cmd, shell=True, cwd=".", stdout=tmp_stdout, stderr=tmp_stderr)
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
        stop_err('Error in running SOAPsnp from (%s), %s' % (opts.consensus_out, str(e)))

    #Clean up temp files
    cleanup_before_exit(tmp_dir)
    #Check results in output file
    if os.path.getsize(opts.consensus_out) > 0:
        sys.stdout.write('Status complete')
    else:
        stop_err("The output is empty")

if __name__ == "__main__": __main__()
