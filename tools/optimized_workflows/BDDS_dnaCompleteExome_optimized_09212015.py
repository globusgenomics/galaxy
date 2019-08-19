# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""
Runs the Meissnert SNPiR RNA-variant discovery pipeline in optimized mode in one Node
See below for options
"""

import time, optparse, os, shutil, subprocess, sys, tempfile, multiprocessing, gzip
CHUNK_SIZE = 2**20 #1mb

def time_stamp():
    return time.strftime('%d/%m/%Y %H:%M:%S\n', time.localtime(time.time()))


def file_type(filename):
    magic_dict = {
        "\x1f\x8b\x08": "gz",
        "\x42\x5a\x68": "bz2",
        "\x50\x4b\x03\x04": "zip"
        }

    max_len = max(len(x) for x in magic_dict)

    with open(filename) as f:
        file_start = f.read(max_len)
    for magic, filetype in magic_dict.items():
        if file_start.startswith(magic):
            return filetype
    return "txt"

def get_ncores():
    return multiprocessing.cpu_count()

def get_linuxRAM():
    # return the available memory on the instance in MB
    totalMemory = os.popen("free -m").readlines()[1].split()[1]
    return int(totalMemory)

def get_readLength(input_fastq):
    if file_type(input_fastq) == "gz":
        fh_zip = gzip.open(input_fastq)
        fh_zip.readline()
        return len(fh_zip.readline().rstrip()) 
    else:
        return len(os.popen("head -n 2 %s" % input_fastq).readlines()[1].rstrip())

def get_sampleName(input_fastq):
    if input_fastq.endswith("fastq.gz"):
        name = os.path.basename(input_fastq)
        return (".").join(name.split(".")[:-2])
    else:
        name = os.path.basename(input_fastq)
        return (".").join(name.split(".")[:-1])

def run_cmd_parallel( cmd, wd_tmpdir, descriptor):
    print cmd
    stderr = tempfile.NamedTemporaryFile( prefix = "cmd_stderr" )
    p = subprocess.Popen(args=cmd, stderr=stderr, shell=True)
    return p

def run_cmd ( cmd , wd_tmpdir, descriptor):
    print cmd
    stderr = tempfile.NamedTemporaryFile( prefix = "cmd_stderr" )
    proc = subprocess.Popen( args=cmd, stderr=stderr, shell=True )

    exit_code = proc.wait()

    if exit_code:
        stderr_target = sys.stderr
    else:
        stderr_target = sys.stdout
    stderr.flush()
    stderr.seek(0)
    while True:
        chunk = stderr.read( CHUNK_SIZE )
        if chunk:
            stderr_target.write( chunk )
        else:
            break
    stderr.close()

def create_wf_meta(options, output_dir, tmp_dir):

    # reads
    bam = options.bam_input

    # get ncores
    ncpu = get_ncores()
    nthreads = ncpu/2

    # get memory
    jvm_heap = get_linuxRAM()

    # get the reference datasets
    reference_db = options.gatk_ref
    if "Hsapiens/human_g1k_v37" in reference_db:
        dbsnp="/mnt/galaxyIndices/genomes/Hsapiens/human_g1k_v37/annotation/dbsnp_138.b37.vcf"
        mills="/mnt/galaxyIndices/genomes/Hsapiens/human_g1k_v37/annotation/Mills_and_1000G_gold_standard.indels.b37.vcf"
    elif "Hsapiens/hg19" in reference_db:
        dbsnp="/mnt/galaxyIndices/genomes/Hsapiens/hg19/annotation/dbsnp_138.hg19.vcf"
        mills="/mnt/galaxyIndices/genomes/Hsapiens/hg19/annotation/Mills_and_1000G_gold_standard.indels.hg19.vcf"

    GATK3_PATH = "/mnt/galaxyTools/tools/gatk3/GenomeAnalysisTK-3.4-46"
    #JAVA_JAR_PATH = os.environ['JAVA_JAR_PATH']
    JAVA_JAR_PATH = "/mnt/galaxyTools/tools/picard/1.134"

    ## Setup workflow
    command_meta = {}
    sample_name = "exome_seq"

    #
    step = 0
    command_meta[step] = []
    inputF = bam
    outputF = "%s/%s-bwa-out.sort.bam" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/sambamba/sambamba_sort.py --memory %sM --input=%s --order=coordinate --output=%s" % (jvm_heap, inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 1
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-bwa-out.sort.rd.bam" % (output_dir, sample_name)
    outputMetrics = "%s/%s-SM1.dups" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF, outputMetrics]
    cmd = "python /opt/galaxy/tools/picard/picard_wrapper.py --maxjheap %sm -i %s -n Dupes_Marked --tmpdir %s -o %s --remdups true --assumesorted true --readregex \"[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*\" --optdupdist 100 -j \"%s/picard.jar MarkDuplicates\" -d %s -t %s -e bam" % (jvm_heap, inputF, tmp_dir, outputF, JAVA_JAR_PATH, tmp_dir, outputMetrics)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 2
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-bwa-out.sort.rd.reorder.bam" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/picard/picard_wrapper.py --input=%s --ref=%s --allow-inc-dict-concord=false --allow-contig-len-discord=false --output-format=bam --output=%s --tmpdir %s -j \"%s/picard.jar ReorderSam\"" % (inputF, reference_db, outputF, tmp_dir, JAVA_JAR_PATH)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 3
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-bwa-out.sort.rd.reorder.bam.bai" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "samtools index %s" % (inputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 4
    command_meta[step] = []
    inputF = command_meta[step-2][0]["output_file"]
    inputIndex = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-output.intervals" % (output_dir, sample_name)
    outputLog = "%s/%s-output.intervals.log" % (output_dir, sample_name)
    input_files = [inputF, inputIndex]
    output_files = [outputF, outputLog]
    cmd = "python /opt/galaxy/tools/gatk3/gatk3_wrapper.py --max_jvm_heap_fraction 1 --stdout %s -d \"-I\" \"%s\" \"bam\" \"gatk_input\" -d \"\" \"%s\" \"bam_index\" \"gatk_input\" -p \"java -jar %s/GenomeAnalysisTK.jar -T RealignerTargetCreator -ip %s -o %s --num_threads %s -R %s\"  -d \"--intervals\" \"%s\" \"bed\" \"input_intervals_0\" -d \"-known:dbsnp,%%(file_type)s\" \"%s\" \"vcf\" \"input_dbsnp_0\" -d \"-known:indels,%%(file_type)s\" \"%s\" \"vcf\" \"input_indels_0\"" % (outputLog, inputF, inputIndex, GATK3_PATH, options.ip, outputF, ncpu, reference_db, options.bed, dbsnp, mills)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 5
    command_meta[step] = []
    inputBam = command_meta[step-3][0]["output_file"]
    inputF = command_meta[step-1][0]["output_file"]
    inputIndex = command_meta[step-2][0]["output_file"]
    outputF = "%s/%s-bwa-out.sort.rd.reorder.realigned.bam" % (output_dir, sample_name)
    outputLog = "%s/%s-bwa-out.sort.rd.reorder.realigned.log" % (output_dir, sample_name)
    input_files = [inputF, inputIndex, inputBam]
    output_files = [outputF, outputLog]
    cmd = "python /opt/galaxy/tools/gatk3/gatk3_wrapper.py --max_jvm_heap_fraction 1 --stdout %s -d \"-I\" \"%s\" \"bam\" \"gatk_input\" -d \"\" \"%s\" \"bam_index\" \"gatk_input\" -p \"java -jar %s/GenomeAnalysisTK.jar -T IndelRealigner -ip %s -o %s -R %s -LOD 5.0\" -d \"--intervals\" \"%s\" \"bed\" \"input_intervals_0\" -d \"-known:dbsnp,%%(file_type)s\" \"%s\" \"vcf\" \"input_dbsnp_0\" -d \"-known:indels,%%(file_type)s\" \"%s\" \"vcf\" \"input_indels_0\"  -d \"-targetIntervals\" \"%s\" \"gatk_interval\" \"gatk_target_intervals\"" % (outputLog, inputBam, inputIndex, GATK3_PATH, options.ip, outputF, reference_db, options.bed, dbsnp, mills, inputF )
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 6
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-bwa-out.sort.rd.reorder.realigned.bam.bai" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "samtools index %s" % (inputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 7
    command_meta[step] = []
    inputBam = command_meta[step-2][0]["output_file"]
    inputIndex = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-recal_data.table" % (output_dir, sample_name)
    outputLog = "%s/%s-bwa-out.sort.rd.reorder.realigned.log" % (output_dir, sample_name)
    input_files = [inputBam, inputIndex]
    output_files = [outputF, outputLog]
    cmd = "python /opt/galaxy/tools/gatk3/gatk3_wrapper.py --max_jvm_heap_fraction 1 --stdout %s -d \"-I\" \"%s\" \"bam\" \"gatk_input\" -d \"\"  \"%s\" \"bam_index\" \"gatk_input\" -p \"java -jar %s/GenomeAnalysisTK.jar -T BaseRecalibrator -ip %s -nct %s -R %s --out %s\" -d \"--intervals\" \"%s\" \"bed\" \"input_intervals_0\" -d \"--knownSites:dbsnp,%%(file_type)s\" \"%s\" \"vcf\" \"input_dbsnp_0\" -d \"--knownSites:indels,%%(file_type)s\" \"%s\" \"vcf\" \"input_indels_0\"" % (outputLog, inputBam, inputIndex, GATK3_PATH, options.ip, ncpu, reference_db, outputF, options.bed, dbsnp, mills )
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 8
    command_meta[step] = []
    inputBam = command_meta[step-3][0]["output_file"]
    inputTable = command_meta[step-1][0]["output_file"]
    inputIndex = command_meta[step-2][0]["output_file"]
    #outputF = "%s/%s-bwa-out.sort.rd.realigned.recal.bam" % (output_dir, sample_name)
    outputF = options.output_bam
    outputLog = "%s/%s-bwa-out.sort.rd.reorder.realigned.recal.log" % (output_dir, sample_name)
    input_files = [inputTable, inputBam, inputIndex]
    output_files = [outputF, outputLog]
    cmd = "python /opt/galaxy/tools/gatk3/gatk3_wrapper.py --max_jvm_heap_fraction 1 --stdout %s -d \"-I\" \"%s\" \"bam\" \"gatk_input\" -d \"\" \"%s\" \"bam_index\" \"gatk_input\" -d \"--intervals\" \"%s\" \"bed\" \"input_intervals_0\" -p \"java -jar %s/GenomeAnalysisTK.jar -T PrintReads -ip %s -o %s -nct %s -R %s --BQSR %s\"" % (outputLog, inputBam, inputIndex, options.bed, GATK3_PATH, options.ip, outputF, ncpu, reference_db, inputTable)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 9
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-bwa-out.sort.rd.reorder.realigned.recal.bam.bai" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "samtools index %s" % (inputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 10
    command_meta[step] = []
    inputBam = command_meta[step-2][0]["output_file"]
    inputIndex = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-raw_variants.vcf" % (output_dir, sample_name)
    outputLog = "%s/%s-raw_variants.log" % (output_dir, sample_name)
    input_files = [inputBam, inputIndex]
    output_files = [outputF, outputLog, inputIndex]
    cmd = "python /opt/galaxy/tools/gatk3/gatk3_wrapper.py --max_jvm_heap %sm --stdout %s -d \"-I\" \"%s\" \"bam\" \"gatk_input_0\"  -d \"\" \"%s\" \"bam_index\" \"gatk_input_0\" -d \"--dbsnp:dbsnp,%%(file_type)s\" \"%s\" \"vcf\" \"input_dbsnp_0\" -d \"--intervals\" \"%s\" \"bed\" \"input_intervals_0\" -p \"java -jar %s/GenomeAnalysisTK.jar -T HaplotypeCaller -ip %s -nct %s --out %s -R %s --output_mode EMIT_VARIANTS_ONLY\"" % (jvm_heap, outputLog, inputBam, inputIndex, dbsnp, options.bed, GATK3_PATH, options.ip, ncpu, outputF, reference_db)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    step = 11
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    #outputF = "%s/%s-raw_variants.filter.select.vcf" % (output_dir, sample_name)
    outputF = options.output_vcf
    outputLog = "%s/%s-raw_variants.filter.select.log" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF, outputLog]
    cmd = "python /opt/galaxy/tools/gatk3/gatk3_wrapper.py --max_jvm_heap %sm --stdout %s -d \"--variant:variant,%%(file_type)s\" \"%s\" \"vcf\" \"input_variant\" -d \"--intervals\" \"%s\" \"bed\" \"input_intervals_0\" -p \"java -jar %s/GenomeAnalysisTK.jar -T SelectVariants -ip %s --num_threads %s -o %s -R %s\"" % (jvm_heap, outputLog, inputF, options.bed, GATK3_PATH, options.ip, ncpu, outputF, reference_db)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    step = 12
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-raw_variants.filter.select.table.txt" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/gatk3/gatk3_wrapper.py --max_jvm_heap_fraction %s -d \"--variant:variant,%%(file_type)s\" \"%s\" \"vcf\" \"gatk_input\" -p \"java -jar %s/GenomeAnalysisTK.jar -T VariantsToTable -o %s -R %s -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F AC -F AN -F DB -F DP -F MQ -GF GT -GF PL -GF GQ --allowMissingData\"" % (jvm_heap, inputF, GATK3_PATH, outputF, reference_db)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    step = 13
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = options.output_csv
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/globus/vcfTable_to_bddsCsv.py -i %s -o %s" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})


    return command_meta

def job_step_dependency(meta, current_step_id):
    step_meta = meta[current_step_id]
    input_files = []
    for job in step_meta:
        for ifile in job['input_files']:
            if ifile not in input_files:
                input_files.append(ifile)

    #print "INPUT_FILES: %s" % input_files

    keep_file = []
    for input_file in input_files:
        flag = None
        for step_id, step_meta in meta.iteritems():
            if step_id > current_step_id:
                for job in step_meta:
                    if input_file in job['input_files']:
                        keep_file.append(input_file)
                        flag = 1
                        break
                if flag:
                    #print "Can't delete file uses in step %s: %s" % (step_id, input_file)
                    break
    return keep_file

def delete_file(input_files, keep_files, step_id, final_outputs):
    for ifile in input_files:
        if ifile not in keep_files and ifile not in final_outputs:
            print "DELETING INPUT: %s" % ifile
            #os.remove(ifile)

def __main__():
    descr = "snpir_optimized_workflow.py: Runs the Meissnert SNPiR RNA-variant discovery pipeline in optimized mode in one Node"
    parser = optparse.OptionParser(description=descr)
    parser.add_option( '--bwa-ref', dest="bwa_ref", help='The BWA reference genome to use or index' )
    parser.add_option( '--gatk-ref', dest="gatk_ref", help='The GATK reference genome to use or index' )
    parser.add_option( '-f', '--fastq', help='The (forward) fastq file to use for the mapping' )
    parser.add_option( '-F', '--rfastq', help='The reverse fastq file to use for mapping if paired-end data' )
    parser.add_option( '', '--bam', dest="bam_input", help='The bam input' )
    parser.add_option( '', '--output-vcf', dest="output_vcf", help="The file to save the output (BAM format)" )
    parser.add_option( '', '--output-csv', dest="output_csv", help="The file to save the output (CSV format)" )
    parser.add_option( '-u', '--output-bam', dest="output_bam", help='The file to save the output (VCF format)' )
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
    parser.add_option( '-b', '--bed', dest="bed", help='The bed file' )
    parser.add_option( '', '--ip', dest="ip", help='The padding for GATK command' )
    parser.add_option( '--username', dest="username", help='Globus username for file transfer' )
    parser.add_option( '--goauth-token', dest="goauth_token", help='Globus token' )
    parser.add_option( '--source-ep', dest="source_ep", help='Globus source endpoint' )
    parser.add_option( '--destination-ep', dest="destination_ep", help='Globus destination endpoint' )
    parser.add_option( '--destination-dir', dest="destination_dir", help='Globus destination dir' )
    parser.add_option( '--deadline', dest="deadline", help='Globus transfer deadline' )

    (options, args) = parser.parse_args()

    final_outputs = [options.output_bam, options.output_log, options.output_vcf]
    if len(args) > 0:
        parser.error('Wrong number of arguments')

#    # make temp directory 
#    output_dir = tempfile.mkdtemp(dir="/scratch/galaxy/tmp", prefix="optimized-dnaExome")
#    tmp_dir = tempfile.mkdtemp(dir=output_dir, prefix="optimized-tmp-")
    # make temp directory 
    output_dir = tempfile.mkdtemp(prefix="optimized-")
    tmp_dir = tempfile.mkdtemp(prefix="optimized-tmp-")

    # create the meta dictionary of the commands to run
    wf_meta = create_wf_meta(options, output_dir, tmp_dir)

    #print wf_meta
    for step_id, step_list in wf_meta.iteritems():
        print "START STEP %s: %s" % (step_id, time_stamp())
        step_input_files = []
        if len(step_list) == 1:
            for job in step_list:
                print "run step %s:\n%s" % (step_id, job['cl'])
                run_cmd ( job['cl'], tmp_dir, "running some job")
                for ifile in job['input_files']:
                    if ifile not in step_input_files:
                        step_input_files.append(ifile)
        else:   # run jobs in parallel
            ps = {}
            for job in step_list:
                print "run step %s:\n%s" % (step_id, job['cl'])
                p = run_cmd_parallel ( job['cl'], tmp_dir, "running some job")
                ps[p.pid] = p
            print "Waiting for %d processes..." % len(ps)
            while ps:
                pid, status = os.wait()
                if pid in ps:
                    del ps[pid]
                    print "Waiting for %d processes..." % len(ps)
        print "END STEP %s: %s" % (step_id, time_stamp())


        # job completed, check to see if there are inputs that should be deleted
        if step_id > 0:
            keep_files = job_step_dependency(wf_meta, step_id)
            delete_file(step_input_files, keep_files, step_id, final_outputs)


    #clean up
    #shutil.rmtree(output_dir)

if __name__ == "__main__":
    __main__()
