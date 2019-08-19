# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""
Runs the ITMI QC pipeline in optimized mode in one Node
See below for options
"""

import workflows, time, optparse, os, shutil, subprocess, sys, tempfile, multiprocessing, gzip
CHUNK_SIZE = 2**20 #1mb

def create_wf_meta(options, output_dir, tmp_dir):

    # reads
    sample_name = options.sample_name

    #lcm = 1
    lcm = 4

    # get ncores
    ncpu = workflows.get_ncores()/lcm
    nthreads = ncpu/2

    # get memory
    jvm_heap = workflows.get_linuxRAM()/lcm

    # get the reference datasets
    picard_reference_db = options.picard_ref
    JAVA_JAR_PATH = "/mnt/galaxyTools/tools/picard/1.134"
    fastqc_path = "/mnt/galaxyTools/tools/FastQC/0.11.3/"
    input_bam = options.input_bam
    input_vcf = options.input_vcf
    input_tabix = options.input_tabix

    ## Setup workflow
    command_meta = {}

    #
    step = 0
    command_meta[step] = []
    inputF = input_bam
    #outputF = "%s/%s-.picard.alignmentSummaryMetrics" % (output_dir, sample_name)
    outputF = options.output_picard_alignmentSummary
    outputMetricsDir = options.output_picard_alignmentSummary_directory
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/picard/picard_wrapper.py -i %s -d %s -t %s --assumesorted true -b false --adaptors \"\" --maxinsert 100000 -n \"Picard Alignment Summary Metrics\" --datatype bam -j \"%s/picard.jar CollectAlignmentSummaryMetrics\" --tmpdir %s --ref %s; cp %s/CollectAlignmentSummaryMetrics.metrics.txt %s" % (inputF, outputMetricsDir, outputF, JAVA_JAR_PATH, tmp_dir, picard_reference_db, outputMetricsDir, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 1
    command_meta[step] = []
    inputF = input_vcf
    inputF2 = input_tabix
    outputF = "%s/%s.vcftools_tstv.txt" % (output_dir, sample_name)
    outputFrest = options.output_vcftools_tstv_stats
    input_files = [inputF, inputF2]
    output_files = [outputF]
    cmd = "ln -s %s input.vcf.gz; ln -s %s input.vcf.tbi; vcftools --gzvcf input.vcf.gz --TsTv-summary > %s 2> /dev/null; cp %s/out.TsTv.summary %s" % (inputF, inputF2, outputF, output_dir, outputFrest)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 2
    command_meta[step] = []
    inputF = input_vcf
    inputF2 = input_tabix
    outputF = "%s/%s.vcflib_hethom_ratio.txt" % (output_dir, sample_name)
    outputFrest = options.output_vcflib_hethom_ratio
    input_files = [inputF, inputF2]
    output_files = [outputF]
    cmd = "cp %s input2.vcf.gz; cp %s input2.vcf.tbi; gunzip input2.vcf.gz; vcfhethomratio input2.vcf > %s; cp %s %s" % (inputF, inputF2, outputF, outputF, outputFrest)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 3
    command_meta[step] = []
    inputF = input_bam
    input_files = [inputF]
    outputHtml = "%s/%s.fastqc.html" % (output_dir, sample_name)
    outputHtmlRest = options.output_fastqc_html
    outputDir = "%s/%s.fastqc" % (output_dir, sample_name)
    outputText = "%s/%s.fastqc.txt" % (output_dir, sample_name)
    outputTextRest = options.output_fastqc_text
    output_files = [outputHtml, outputText]
    cmd = "python /opt/galaxy/tools/rgenetics/rgFastQC.py -i %s -d %s -o %s -t %s -f bam -j BAM_Input_File -e %s/fastqc; cp %s %s; cp %s %s" % (inputF, outputDir, outputHtml, outputText, fastqc_path, outputHtml, outputHtmlRest, outputText, outputTextRest)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    return command_meta

def __main__():
    descr = "snpir_optimized_workflow.py: Runs the Meissnert SNPiR RNA-variant discovery pipeline in optimized mode in one Node"
    parser = optparse.OptionParser(description=descr)
    parser.add_option( '--sample-name', dest="sample_name", help='sample name' )
    parser.add_option( '--picard-ref', dest="picard_ref", help='The picard reference genome to use or index' )
    parser.add_option( '', '--input-bam', dest="input_bam", help='The path to input bam' )
    parser.add_option( '', '--input-vcf', dest="input_vcf", help='The path to input vcf' )
    parser.add_option( '', '--input-tabix', dest="input_tabix", help='The path to input tabix' )
    parser.add_option( '--output-vcflib-hethomratio', dest="output_vcflib_hethom_ratio", help='The file to save' )
    parser.add_option( '--output-vcftools-tstvstats', dest="output_vcftools_tstv_stats", help='The file to save' )
    parser.add_option( '--output-picard-alignmentSummary', dest="output_picard_alignmentSummary", help='The file to save' )
    parser.add_option( '--output-picard-alignmentSummary-directory', dest="output_picard_alignmentSummary_directory", help='The dir output' )
    parser.add_option( '--output-fastqc-html', dest="output_fastqc_html", help="File to save" )
    parser.add_option( '--output-fastqc-text', dest="output_fastqc_text", help="File to save" )
    parser.add_option( '--output-log', dest="output_log", help='The file output' )

    parser.add_option("-u", "--username", dest="username", help="username")
    parser.add_option("-c", "--cert", dest="cert_file", help="client cert file", metavar="CERT_FILE")
    parser.add_option("-k", "--key", dest="key_file", help="client key file", metavar="KEY_FILE")
    parser.add_option("-b", "--base-url", dest="base_url", help="alternate base URL", metavar="URL")
    parser.add_option("-o", "--out-transfer-log", dest="output_transfer_log", help="write log output to PATH", metavar="PATH")
    parser.add_option("--source-ep", dest="source_ep", help="Endpoint to transfer from")
    parser.add_option("--destination-ep", dest="destination_ep", help="Endpoint to transfer to")
    parser.add_option("--source-path-bam", dest="source_path_bam", help="BAM: Source endpoint filepath to transfer")
    parser.add_option("--destination-path-bam", dest="destination_path_bam", help="BAM: Destination endpoint filepath to transfer")
    parser.add_option("--source-path-vcf", dest="source_path_vcf", help="VCF: Source endpoint filepath to transfer")
    parser.add_option("--destination-path-vcf", dest="destination_path_vcf", help="VCF: Destination endpoint filepath to transfer")
    parser.add_option("--source-path-tbi", dest="source_path_tbi", help="TBI: Source endpoint filepath to transfer")
    parser.add_option("--destination-path-tbi", dest="destination_path_tbi", help="TBI: Destination endpoint filepath to transfer")

    parser.add_option("-d", "--deadline", dest="deadline", help="Deadline for transfer in minutes.")
    parser.add_option('-a', '--goauth-token', dest="goauth_token", help="Use the Globus Access Token as the authentication method for this transfer")
    parser.add_option("--type", dest="path_type", help="directory or file")
    (options, args) = parser.parse_args()

    final_outputs = [options.output_picard_alignmentSummary, options.output_picard_alignmentSummary_directory, options.output_vcflib_hethom_ratio, options.output_vcftools_tstv_stats, options.output_fastqc_html, options.output_fastqc_text, options.output_log ]


    if len(args) > 0:
        parser.error('Wrong number of arguments')

#    # make temp directory 
#    output_dir = tempfile.mkdtemp(dir="/scratch/galaxy/tmp", prefix="optimized")
#    tmp_dir = tempfile.mkdtemp(dir=output_dir, prefix="optimized-tmp-")
    # make temp directory 
    output_dir = tempfile.mkdtemp(prefix="optimized-")
    tmp_dir = tempfile.mkdtemp(prefix="optimized-tmp-")

    # create the meta dictionary of the commands to run
    wf_meta = create_wf_meta(options, output_dir, tmp_dir)

    #print wf_meta
    os.chdir(output_dir)
    for step_id, step_list in wf_meta.iteritems():
        print "START STEP %s: %s" % (step_id, workflows.time_stamp())
        step_input_files = []
        if len(step_list) == 1:
            for job in step_list:
                print "run step %s:\n%s" % (step_id, job['cl'])
                workflows.run_cmd ( job['cl'], tmp_dir, "running some job")
                for ifile in job['input_files']:
                    if ifile not in step_input_files:
                        step_input_files.append(ifile)
        else:   # run jobs in parallel
            ps = {}
            for job in step_list:
                print "run step %s:\n%s" % (step_id, job['cl'])
                p = workflows.run_cmd_parallel ( job['cl'], tmp_dir, "running some job")
                ps[p.pid] = p
            print "Waiting for %d processes..." % len(ps)
            while ps:
                pid, status = os.wait()
                if pid in ps:
                    del ps[pid]
                    print "Waiting for %d processes..." % len(ps)
        print "END STEP %s: %s" % (step_id, workflows.time_stamp())


        # job completed, check to see if there are inputs that should be deleted
        if step_id > 0:
            keep_files = workflows.job_step_dependency(wf_meta, step_id)
            #workflows.delete_file(step_input_files, keep_files, step_id, final_outputs)


    #clean up
    #shutil.rmtree(output_dir)

if __name__ == "__main__":
    __main__()
