# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""
Runs the whole genome pipeline up to haplotype caller 
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
    JAVA_JAR_PATH = "/mnt/galaxyTools/tools/picard/1.134"

    bwa_index = options.bwa_ref
    ## Setup workflow
    command_meta = {}

    read_group = "--rgid=\"%s\" --rglb=\"%s\" --rgpl=\"%s\" --rgsm=\"%s\"" % (options.rgid, options.rglb, options.rgpl, options.rgsm)
    if options.rgcn:
        read_group += " --rgcn=\"%s\"" % options.rgcn
    if options.rgds:
        read_group += " --rgds=\"%s\"" % options.rgds
    if options.rgdt:
        read_group += " --rgdt=\"%s\"" % options.rgdt
    if options.rgfo:
        read_group += " --rgfo=\"%s\"" % options.rgfo
    if options.rgks:
        read_group += " --rgks=\"%s\"" % options.rgks
    if options.rgpg:
        read_group += " --rgpg=\"%s\"" % options.rgpg
    if options.rgpi:
        read_group += " --rgpi=\"%s\"" % options.rgpi
    if options.rgpu:
        read_group += " --rgpu=\"%s\"" % options.rgpu

    input_files = sorted(glob.glob("%s/*.fastq.gz" % input_dir))
    forward, reverse = _sort_input_fastq(input_dir)

    #
    step = 0
    command_meta[step] = []
    input_files = forward
    outputF1 = "%s/%s-R1.fastq.gz" % (output_dir, sample_name)
    fastq = "%s/%s-R1.fastq.gz" % (output_dir, sample_name)
    output_files = [outputF1]
    cmd = "python /opt/galaxy/tools/filters/catWrapper.py %s %s" % (outputF1, " ".join(input_files))
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF1, "output_files": output_files})

    #step = 0
    #command_meta[step] = []
    input_files = reverse
    outputF1 = "%s/%s-R2.fastq.gz" % (output_dir, sample_name)
    rfastq = "%s/%s-R2.fastq.gz" % (output_dir, sample_name)
    output_files = [outputF1]
    cmd = "python /opt/galaxy/tools/filters/catWrapper.py %s %s" % (outputF1, " ".join(input_files))
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF1, "output_files": output_files})

    #
    step = 1
    command_meta[step] = []
    ncpu_step = int(ncpu) - 2
    input1 = fastq
    input2 = rfastq
    input_files = [input1, input2]
    outputF = "%s/%s-bwa-out-sort.bam" % (output_dir, sample_name)
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/sr_mapping/bwa_wrapper.py --threads=%s --fileSource=indexed --ref=%s --do_not_build_index --input1=%s --input2=%s --output=%s --output-format=BAM --genAlignType=paired --params=pre_set %s; rm -rf %s %s" % (ncpu_step, bwa_index, input1, input2, outputF, read_group, input1, input2)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files": output_files})

    #
    step = 2
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-bwa-out.sort-fixmate.bam" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "samtools fixmate %s %s; rm -rf %s" % (inputF, outputF, inputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 3
    command_meta[step] = []
    ncpu_step = int(ncpu) - 2
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-bwa-out.sort-fixmate-filtered.bam" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "samtools view -@ %s -bF 4 %s > %s; rm -rf %s" % (str(ncpu_step), inputF, outputF, inputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 4
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-bwa-out.sort-fixmate-filtered-sorted.bam" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/sambamba/sambamba_sort.py --memory %sM --input=%s --order=coordinate --output=%s; rm -rf %s" % (jvm_heap, inputF, outputF, inputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 5
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-bwa-out.sort.dedup.bam" % (output_dir, sample_name)
    outputMetrics = "%s/%s-SM1.dups" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF, outputMetrics]
    cmd = "python /opt/galaxy/tools/picard/picard_wrapper.py --maxjheap %sm -i %s -n Dupes_Marked --tmpdir %s -o %s --remdups false --assumesorted true --readregex \"[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*\" --optdupdist 100 -j \"%s/picard.jar MarkDuplicates\" -d %s -t %s -e bam; rm -rf %s" % (jvm_heap, inputF, tmp_dir, outputF, JAVA_JAR_PATH, tmp_dir, outputMetrics, inputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 6
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-bwa-out.sort.rd.bam.bai" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "samtools index %s" % (inputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 7
    command_meta[step] = []
    inputF = command_meta[step-2][0]["output_file"]
    inputIndex = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-output.intervals" % (output_dir, sample_name)
    outputLog = "%s/%s-output.intervals.log" % (output_dir, sample_name)
    input_files = [inputF, inputIndex]
    output_files = [outputF, outputLog]
    cmd = "python /opt/galaxy/tools/gatk3/gatk3_wrapper.py --max_jvm_heap_fraction 1 --stdout %s -d \"-I\" \"%s\" \"bam\" \"gatk_input\" -d \"\" \"%s\" \"bam_index\" \"gatk_input\" -p \"java -jar %s/GenomeAnalysisTK.jar -T RealignerTargetCreator -o %s --num_threads %s -R %s\"" % (outputLog, inputF, inputIndex, GATK3_PATH, outputF, ncpu, reference_db)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 8
    command_meta[step] = []
    inputBam = command_meta[step-3][0]["output_file"]
    inputF = command_meta[step-1][0]["output_file"]
    inputIndex = command_meta[step-2][0]["output_file"]
    outputF = "%s/%s-bwa-out.sort.rd.realigned.bam" % (output_dir, sample_name)
    outputLog = "%s/%s-bwa-out.sort.rd.realigned.log" % (output_dir, sample_name)
    input_files = [inputF, inputIndex, inputBam]
    output_files = [outputF, outputLog]
    cmd = "python /opt/galaxy/tools/gatk3/gatk3_wrapper.py --max_jvm_heap_fraction 1 --stdout %s -d \"-I\" \"%s\" \"bam\" \"gatk_input\" -d \"\" \"%s\" \"bam_index\" \"gatk_input\" -p \"java -jar %s/GenomeAnalysisTK.jar -T IndelRealigner -o %s -R %s\" -d \"-targetIntervals\" \"%s\" \"gatk_interval\" \"gatk_target_intervals\"; rm -rf %s %s %s" % (outputLog, inputBam, inputIndex, GATK3_PATH, outputF, reference_db, inputF, inputBam, inputIndex, inputF )
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 9
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-bwa-out.sort.rd.realigned.bam.bai" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "samtools index %s" % (inputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 10
    command_meta[step] = []
    inputBam = command_meta[step-2][0]["output_file"]
    inputIndex = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-recal_data.table" % (output_dir, sample_name)
    outputLog = "%s/%s-bwa-out.sort.rd.realigned.log" % (output_dir, sample_name)
    input_files = [inputBam, inputIndex]
    output_files = [outputF, outputLog]
    cmd = "python /opt/galaxy/tools/gatk3/gatk3_wrapper.py --max_jvm_heap_fraction 1 --stdout %s -d \"-I\" \"%s\" \"bam\" \"gatk_input\" -d \"\"  \"%s\" \"bam_index\" \"gatk_input\" -p \"java -jar %s/GenomeAnalysisTK.jar -T BaseRecalibrator -nct %s -R %s --out %s\" -d \"--knownSites:dbsnp,%%(file_type)s\" \"%s\" \"vcf\" \"input_dbsnp_0\" -d \"--knownSites:indels,%%(file_type)s\" \"%s\" \"vcf\" \"input_indels_0\"" % (outputLog, inputBam, inputIndex, GATK3_PATH, ncpu, reference_db, outputF, dbsnp, mills )
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 11
    command_meta[step] = []
    inputBam = command_meta[step-3][0]["output_file"]
    inputTable = command_meta[step-1][0]["output_file"]
    inputIndex = command_meta[step-2][0]["output_file"]
    outputF = "%s/%s-bwa-out.sort.rd.realigned.recal.bam" % (output_dir, sample_name)
    outputFRest = options.output_bam
    outputLog = "%s/%s-bwa-out.sort.rd.realigned.recal.log" % (output_dir, sample_name)
    input_files = [inputTable, inputBam, inputIndex]
    output_files = [outputF, outputLog, outputFRest]
    cmd = "python /opt/galaxy/tools/gatk3/gatk3_wrapper.py --max_jvm_heap_fraction 1 --stdout %s -d \"-I\" \"%s\" \"bam\" \"gatk_input\" -d \"\" \"%s\" \"bam_index\" \"gatk_input\" -p \"java -jar %s/GenomeAnalysisTK.jar -T PrintReads -o %s -nct %s -R %s --BQSR %s\"; rm -rf %s %s %s; cp %s %s" % (outputLog, inputBam, inputIndex, GATK3_PATH, outputF, ncpu, reference_db, inputTable, inputBam, inputTable, inputIndex, outputF, outputFRest)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 12
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-bwa-out.sort.rd.realigned.recal.bam.bai" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "samtools index %s" % (inputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 13
    command_meta[step] = []
    inputBam = command_meta[step-2][0]["output_file"]
    inputIndex = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-raw_variants.g.vcf" % (output_dir, sample_name)
    outputLog = "%s/%s-raw_variants.g.vcf.log" % (output_dir, sample_name)
    outputFRest = options.output_vcf
    input_files = [inputBam, inputIndex]
    output_files = [outputF, outputLog, outputFRest, inputIndex]
    cmd = "python /opt/galaxy/tools/gatk3/gatk3_wrapper.py --max_jvm_heap %sm --stdout %s -d \"-I\" \"%s\" \"bam\" \"gatk_input_0\"  -d \"\" \"%s\" \"bam_index\" \"gatk_input_0\" -p \"java -jar %s/GenomeAnalysisTK.jar -T HaplotypeCaller -nct %s --out %s -R %s --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000\"; cp %s %s" % (jvm_heap, outputLog, inputBam, inputIndex, GATK3_PATH, ncpu, outputF, reference_db, outputF, outputFRest)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

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

    parser.add_option( '--bwa-ref', dest="bwa_ref", help='The BWA reference genome to use or index' )
    parser.add_option( '--gatk-ref', dest="gatk_ref", help='The GATK reference genome to use or index' )
    parser.add_option( '-f', '--fastq', help='The (forward) fastq file to use for the mapping' )
    parser.add_option( '-F', '--rfastq', help='The reverse fastq file to use for mapping if paired-end data' )
    parser.add_option( '', '--output-bam', dest="output_bam", help='The file to save the output (VCF format)' )
    parser.add_option( '', '--output-vcf', dest="output_vcf", help='The file to save the output (interval format)' )
    parser.add_option( '-l', '--output-log', dest="output_log", help='The log output (Txt format)' )
    parser.add_option( '--rgid', help='Read group identifier' )
    parser.add_option( '--rgsm', help='Sample' )
    parser.add_option( '--rgpl', choices=[ 'CAPILLARY', 'LS454', 'ILLUMINA', 'SOLID', 'HELICOS', 'IONTORRENT', 'PACBIO' ], help='Platform/technology used to produce the reads' )
    parser.add_option( '--rglb', help='Library name' )
    parser.add_option( '--rgpu', help='Platform unit (e.g. flowcell-barcode.lane for Illumina or slide for SOLiD)' )
    parser.add_option( '--rgcn', help='Sequencing center that produced the read' )
    parser.add_option( '--rgds', help='Description' )
    parser.add_option( '--rgdt', help='Date that run was produced (ISO8601 format date or date/time, like YYYY-MM-DD)' )
    parser.add_option( '--rgfo', help='Flow order' )
    parser.add_option( '--rgks', help='The array of nucleotide bases that correspond to the key sequence of each read' )
    parser.add_option( '--rgpg', help='Programs used for processing the read group' )
    parser.add_option( '--rgpi', help='Predicted median insert size' )
    (options, args) = parser.parse_args()

    final_outputs = [options.output_bam, options.output_vcf, options.output_log]
    if len(args) > 0:
        parser.error('Wrong number of arguments')

#    # make temp directory 
    #output_dir = tempfile.mkdtemp(dir="/scratch/galaxy/tmp", prefix="optimized-")
    #tmp_dir = tempfile.mkdtemp(dir="/scratch/galaxy/tmp", prefix="optimized-tmp-")
    # make temp directory 
    output_dir = tempfile.mkdtemp(prefix="optimized-")
    tmp_dir = tempfile.mkdtemp(prefix="optimized-tmp-")

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
