# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""
Runs the whole genome jointGenotyping pipeline 
See below for options
"""

import glob, workflows, time, optparse, os, shutil, subprocess, sys, tempfile, multiprocessing, gzip
CHUNK_SIZE = 2**20 #1mb

def _sort_input_fastq(input_dir):
    forward = []
    reverse = []
    for infile in sorted(os.listdir(input_dir)):
        # assume files are fastq
        if infile.endswith(".fastq.gz"):
            if "_R1_" in infile or "_1_" in infile or "_1.fastq.gz" in infile:
                forward.append("%s/%s" % (input_dir, infile))
            elif "_R2_" in infile or "_2_" in infile or "_2.fastq.gz" in infile:
                reverse.append("%s/%s" % (input_dir, infile))
    return (forward, reverse)

def create_wf_meta(options, output_dir, tmp_dir, input_dir):

    # get sample name
    sample_name = "sampleName"

    # get ncores
    ncpu = workflows.get_ncores()
    nthreads = ncpu/2

    # get memory
    jvm_heap = workflows.get_linuxRAM()

    # get the reference datasets
    reference_db = options.gatk_ref
    dbsnp="/mnt/galaxyIndices/genomes/Hsapiens/hg38/annotation/dbsnp_144.hg38.vcf"
    #dbsnp="/mnt/galaxyIndices/genomes/Hsapiens/hg38_ucsc/annotation/Homo_sapiens_assembly38.dbsnp138.vcf"
    mills="/mnt/galaxyIndices/genomes/Hsapiens/hg38/annotation/Mills_and_1000G_gold_standard.indels.hg38.vcf"
    #mills="/mnt/galaxyIndices/genomes/Hsapiens/hg38_ucsc/annotation/Mills_and_1000G_gold_standard.indels.hg38.vcf"
    GATK3_PATH = "/mnt/galaxyTools/tools/gatk3/GenomeAnalysisTK-3.5"

    ## Setup workflow
    command_meta = {}

    #
    step = 0
    command_meta[step] = []
    input_files = sorted(glob.glob("%s/*/*.[Vv][Cc][Ff]" % input_dir))
    outputF = "%s/%s-joint.vcf" % (output_dir, sample_name)
    outputLog = "%s/%s-joint.log" % (output_dir, sample_name)
    output_files = [outputF, outputLog]
    input_cmd = ""
    number = 0
    for inputF in input_files:
        input_cmd += " -d \"--variant\" \"%s\" \"vcf\" \"gatk_input_%s\"" % (inputF, str(number))
        number += 1
    cmd = "python /opt/galaxy/tools/gatk3/gatk3_wrapper.py --max_jvm_heap %sm --stdout %s %s -p \"java -jar %s/GenomeAnalysisTK.jar -T GenotypeGVCFs --num_threads %s --out %s -R %s --disable_auto_index_creation_and_locking_when_reading_rods\"" % (jvm_heap, outputLog, input_cmd, GATK3_PATH, ncpu, outputF, reference_db)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 1
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-recal.txt" % (output_dir, sample_name)
    outputTranches = "%s/%s-tranches.txt" % (output_dir, sample_name)
    outputRscript = "%s/%s-rscript.txt" % (output_dir, sample_name)
    outputLog = "%s/%s-variantRecalibrator.log" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF, outputTranches, outputRscript, outputLog]
    cmd = "python /opt/galaxy/tools/gatk3/gatk3_wrapper.py --max_jvm_heap %sm --stdout %s -d \"--input\" \"%s\" \"vcf\" \"gatk_input_0\" -p \"java -jar %s/GenomeAnalysisTK.jar -T VariantRecalibrator -nct %s --recal_file %s --tranches_file %s --rscript_file %s -R %s --disable_auto_index_creation_and_locking_when_reading_rods --mode SNP -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 %s \"" % (jvm_heap, outputLog, inputF, GATK3_PATH, ncpu, outputF, outputTranches, outputRscript, reference_db, dbsnp )
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 2
    command_meta[step] = []
    inputV = command_meta[step-2][0]["output_file"]
    inputF = command_meta[step-1][0]["output_files"][0]
    inputT = command_meta[step-1][0]["output_files"][1]
    outputF = "%s/%s-avr.vcf" % (output_dir, sample_name)
    outputFRest = options.output_vcf
    outputLog = "%s/%s-avr.log" % (output_dir, sample_name)
    output_files = [outputF, outputLog]
    input_files = [inputF, inputT, inputV]
    cmd = "python /opt/galaxy/tools/gatk3/gatk3_wrapper.py --max_jvm_heap %sm --stdout %s -d \"--input\" \"%s\" \"vcf\" \"gatk_input_0\" -p \"java -jar %s/GenomeAnalysisTK.jar -T ApplyRecalibration -nct %s --recal_file %s --tranches_file %s -R %s --disable_auto_index_creation_and_locking_when_reading_rods -mode SNP --ts_filter_level 99.0\"; cp %s %s" % (jvm_heap, outputLog, inputV, GATK3_PATH, ncpu, inputF, inputT, reference_db, outputF, outputFRest)

    return command_meta


def __main__():
    descr = "snpir_optimized_workflow.py: Runs the Meissnert SNPiR RNA-variant discovery pipeline in optimized mode in one Node"
    parser = optparse.OptionParser(description=descr)
    parser.add_option("-u", "--username", dest="username", help="username")
    parser.add_option("-c", "--cert", dest="cert_file", help="client cert file", metavar="CERT_FILE")
    parser.add_option("-k", "--key", dest="key_file", help="client key file", metavar="KEY_FILE")
    parser.add_option("-b", "--base-url", dest="base_url", help="alternate base URL", metavar="URL")
    parser.add_option("-o", "--out-transfer-log", dest="output_transfer_log", help="write log output to PATH", metavar="PATH")
    parser.add_option("--source-ep", dest="source_ep", help="Endpoint to transfer from")
    parser.add_option("--source-path", dest="source_path", help="Source endpoint filepath to transfer")
    parser.add_option("--extra-source-path", dest="extra_source_path", help="Source endpoint filepath to transfer for BAM Index file")
    parser.add_option("--destination-ep", dest="destination_ep", help="Endpoint to transfer to")
    parser.add_option("--destination-path", dest="destination_path", help="Destination endpoint filepath to transfer")
    parser.add_option("-d", "--deadline", dest="deadline", help="Deadline for transfer in minutes.")
    parser.add_option("-g", "--galaxy-dataset-id", dest="galaxy_dataset_id", help="Galaxy Dataset Id For This Transfer")
    parser.add_option('-a', '--goauth-token', dest="goauth_token", help="Use the Globus Access Token as the authentication method for this transfer")
    parser.add_option("--type", dest="path_type", help="directory or file")
    parser.add_option( '--input-path', dest="input_dir", help='The path to the input GVCF files', default=None)
    parser.add_option( '--gatk-ref', dest="gatk_ref", help='The GATK reference genome to use or index' )
    parser.add_option( '', '--output-vcf', dest="output_vcf", help='The file to save the output (interval format)' )
    parser.add_option( '-l', '--output-log', dest="output_log", help='The log output (Txt format)' )
    (options, args) = parser.parse_args()

    final_outputs = [options.output_vcf, options.output_log]
    if len(args) > 0:
        parser.error('Wrong number of arguments')

#    # make temp directory 
    #output_dir = tempfile.mkdtemp(dir="/scratch/galaxy/tmp", prefix="optimized-")
    #tmp_dir = tempfile.mkdtemp(dir="/scratch/galaxy/tmp", prefix="optimized-tmp-")
    # make temp directory 
    output_dir = tempfile.mkdtemp(prefix="optimized-")
    tmp_dir = tempfile.mkdtemp(prefix="optimized-tmp-")

    if options.input_dir is not None:
        worker_input_dir = options.input_dir
    else:
        worker_input_dir = "%s/input" % output_dir
        workflows.transfer_data(options, output_dir, worker_input_dir)

    # create the meta dictionary of the commands to run
    wf_meta = create_wf_meta(options, output_dir, tmp_dir, worker_input_dir)

    #print wf_meta
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
            workflows.delete_file(step_input_files, keep_files, step_id, final_outputs)


    #clean up
    shutil.rmtree(output_dir)
    shutil.rmtree(tmp_dir)

if __name__ == "__main__":
    __main__()
