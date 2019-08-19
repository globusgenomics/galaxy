# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""
Runs the HLI Exome-SEQ QC pipeline in optimized mode in one Node
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
    fastq = None
    if os.path.exists(options.inputDir):
        fastq = os.listdir(options.inputDir)[0]

    # get sample name
    sample_name = None
    if not options.sample_name:
        sample_name = get_sampleName(fastq)
    else:
        sample_name = options.sample_name

    lcm = 4

    # get ncores
    ncpu = get_ncores()/lcm
    nthreads = ncpu/2

    # get memory
    jvm_heap = get_linuxRAM()/lcm

    # get the reference datasets
    picard_reference_db = options.picard_ref
    bwa_reference_db = options.bwa_ref
    paths = bwa_reference_db.split("bwa")
    anno_path = "%s/annotation" % paths[0]
    if options.vendor == "Illumina":
        baitfile = "%s/NexteraRapidCapture_Exome_Probes_v1.2_picard.txt" % anno_path
        targetfile = "%s/nexterarapidcapture_exome_targetedregions_v1.2_picard.txt" % anno_path
    elif options.vendor == "IDT":
        baitfile = "%s/xgen-exome-research-panel-probes_picard.bed" % anno_path
        targetfile = "%s/xgen-exome-research-panel-targets_picard.bed" % anno_path
    elif options.vendor == "Agilent":
        baitfile = "%s/S07604514_Regions_picard.bed" % anno_path
        targetfile = "%s/S07604514_Covered_picard.bed" % anno_path
    elif options.vendor == "Roche":
        baitfile = "%s/MedExome_capture_targets_picard.bed" % anno_path
        targetfile = "%s/MedExome_empirical_targets_picard.bed" % anno_path

    #baitfile = options.bait_file
    #targetfile = options.target_file
    refseqfile = options.refseq_file
    genetargetfile = options.gene_targetfile
    #JAVA_JAR_PATH = os.environ['JAVA_JAR_PATH']
    JAVA_JAR_PATH = "/mnt/galaxyTools/tools/picard/1.134"
    qualimap_jar = "/mnt/galaxyTools/tools/qualimap/qualimap_v2.1.1"
    GATK_JAR_PATH = "/mnt/galaxyTools/tools/gatk3/GenomeAnalysisTK-3.3-0"
    read_group = "--rgid=\"%s\" --rglb=\"%s\" --rgpl=\"%s\" --rgsm=\"%s\"" % (options.rgid, options.rglb, options.rgpl, options.rgsm)

    ## Setup workflow
    command_meta = {}

    #
    step = 0
    command_meta[step] = []
    inputF = options.inputDir
    input_files = [inputF]
    outputF1 = "%s/%s-R1.fastq.gz" % (output_dir, sample_name)
    outputF2 = "%s/%s-R2.fastq.gz" % (output_dir, sample_name)
    output_files = [outputF1, outputF2]
    cmd = "python /opt/galaxy/tools/filters/fastqcatWrapper.py %s %s %s" % (outputF1, outputF2, inputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF1, "output_files": output_files})

    #
    step = 1
    ncpu_step = int(ncpu)
    command_meta[step] = []
    input1 = command_meta[0][0]["output_files"][0]
    input2 = command_meta[0][0]["output_files"][1]
    input_files = [input1, input2]
    outputF = "%s/%s-bwa-out.bam" % (output_dir, sample_name)
    outputFrest = options.output_bwa_bam
    #outputLog = "%s/%s-tophat.log" % (output_dir, sample_name)
    output_files = [outputF, options.output_bwa_bam]
    cmd = "python /opt/galaxy/tools/sr_mapping/bwa_mem.py --threads=%s --fileSource=indexed --ref=%s --fastq=%s --rfastq=%s --output=%s --genAlignType=paired --params=full --minSeedLength 19 --bandWidth 100 --offDiagonal 100 --internalSeeds 1.5 --seedsOccurrence 10000 --seqMatch 1 --mismatch 4 --gapOpen 6 --gapExtension 1 --clipping 5 --unpairedReadpair 17 --minScore 30 %s -M --bam" % (ncpu_step, bwa_reference_db, input1, input2, outputF, read_group)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files": output_files})

    #
    step = 2
    command_meta[step] = []
    inputF = command_meta[1][0]["output_file"]
    outputF = "%s/%s-bwa-out.rd.bam" % (output_dir, sample_name)
    outputFrest = options.output_bwa_bam
    outputMetrics = options.output_picard_markDuplicates_html
    outputMetricsDir = options.output_picard_markDuplicates_directory
    #outputMetrics = "%s/%s-SM1.dups" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF, outputMetrics]
    cmd = "python /opt/galaxy/tools/picard/picard_wrapper.py --maxjheap %sm -i %s -n Dupes_Marked --tmpdir %s -o %s --assumesorted true --readregex \"[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*\" --remdups false --optdupdist 100 -j \"%s/picard.jar MarkDuplicates\" -d %s -t %s -e bam; cp %s %s" % (jvm_heap, inputF, tmp_dir, outputF, JAVA_JAR_PATH, outputMetricsDir, outputMetrics, outputF, outputFrest)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 3
    command_meta[step] = []
    inputF = command_meta[2][0]["output_file"]
    #outputF = "%s/%s-tophat-out.rd.picard.alignmentSummaryMetrics" % (output_dir, sample_name)
    outputF = options.output_picard_alignmentSummary
    outputMetricsDir = options.output_picard_alignmentSummary_directory
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/picard/picard_wrapper.py -i %s -d %s -t %s --assumesorted true -b false --adaptors \"\" --maxinsert 100000 -n \"Picard Alignment Summary Metrics\" --datatype bam -j \"%s/picard.jar CollectAlignmentSummaryMetrics\"  --tmpdir %s --ref %s" % (inputF, outputMetricsDir, outputF, JAVA_JAR_PATH, tmp_dir, picard_reference_db)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    #step = 4
    #command_meta[step] = []
    inputF = command_meta[2][0]["output_file"]
    #outputF = "%s/%s-tophat-out.rd.picard.insertSizeMetrics" % (output_dir, sample_name)
    outputF = options.output_picard_insertSize
    outputMetricsDir = options.output_picard_insertSize_directory
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/picard/picard_wrapper.py -i %s -n \"Insertion size metrics\" --tmpdir %s --deviations 10.0 --histwidth 0 --minpct 0.05 --malevel ALL_READS -j \"%s/picard.jar CollectInsertSizeMetrics\" -d %s -t %s" % (inputF, tmp_dir, JAVA_JAR_PATH, outputMetricsDir, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 4
    command_meta[step] = []
    inputF = command_meta[2][0]["output_file"]
    #outputF = "%s/%s-tophat-out.rd.picard.gcBiasMetrics" % (output_dir, sample_name)
    outputF = options.output_picard_shortReadGCBias
    outputMetricsDir = options.output_picard_shortReadGCBias_directory
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/picard/picard_wrapper.py -i %s -d %s -t %s --windowsize 100 --mingenomefrac 1e-05 -n \"Short Read GC Bias Metrics\" --tmpdir %s -j \"%s/picard.jar CollectGcBiasMetrics\"  --ref %s" % (inputF, outputMetricsDir, outputF, tmp_dir, JAVA_JAR_PATH, picard_reference_db)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 5
    command_meta[step] = []
    inputF = command_meta[2][0]["output_file"]
    #outputF = "%s/%s-bwa-out.bam.bai" % (output_dir, sample_name)
    outputF = "%s/%s-bwa-out.rd.bam.bai" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "samtools index %s" % (inputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 6
    command_meta[step] = []
    inputF = command_meta[2][0]["output_file"]
    inputIndex = command_meta[5][0]["output_file"]
    outputF = "%s/%s-depthOfCoverage.txt" % (output_dir, sample_name)
    outputGeneSummary = options.output_gatk_depthOfCoverage_geneSummary
    outputIntervalSummary = options.output_gatk_depthOfCoverage_intervalSummary
    input_files = [inputF]
    output_files = [outputGeneSummary]
    cmd = "python /opt/galaxy/tools/gatk3/gatk3_wrapper.py --max_jvm_heap_fraction 1 -d \"-I\" \"%s\" \"bam\" \"gatk_input_0\" -d \"\" \"%s\" \"bam_index\" \"gatk_input_0\" -p \"java -jar %s/GenomeAnalysisTK.jar -T DepthOfCoverage -R %s --calculateCoverageOverGenes %s --partitionType sample --out %s --outputFormat table\" -d \"--intervals\" \"%s\" \"bed\" \"input_intervals_0\"; mv %s.sample_gene_summary %s; mv %s.sample_interval_summary %s" % (inputF, inputIndex, GATK_JAR_PATH, picard_reference_db, refseqfile, outputF, genetargetfile, outputF, outputGeneSummary, outputF, outputIntervalSummary)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 7
    command_meta[step] = []
    inputF = command_meta[2][0]["output_file"]
    inputIndex = command_meta[5][0]["output_file"]
    outputF = options.output_gatk_diagnoseTargets_vcf
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/gatk3/gatk3_wrapper.py --max_jvm_heap_fraction 1 -d \"-I\" \"%s\" \"bam\" \"gatk_input_0\" -d \"\" \"%s\" \"bam_index\" \"gatk_input_0\" -p \"java -jar %s/GenomeAnalysisTK.jar -T DiagnoseTargets -R %s --out %s\" -d \"--intervals\" \"%s\" \"bed\" \"input_intervals_0\"" % (inputF, inputIndex, GATK_JAR_PATH, picard_reference_db, outputF, genetargetfile)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 8
    command_meta[step] = []
    inputF = command_meta[2][0]["output_file"]
    outputF = options.output_picard_HSMetrics
    outputMetricsDir = options.output_picard_HSMetrics_directory
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/picard/picard_wrapper.py -i %s -d %s -t %s --datatype bam --baitbed %s --targetbed %s -n \"Picard HS Metrics\" --tmpdir %s -j \"%s/picard.jar CalculateHsMetrics\"" % (inputF, outputMetricsDir, outputF, baitfile, targetfile, tmp_dir, JAVA_JAR_PATH)
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
    parser.add_option( '--sample-name', dest="sample_name", help='sample name' )
    parser.add_option( '--bwa-ref', dest="bwa_ref", help='The tophat2 reference genome to use or index' )
    parser.add_option( '--picard-ref', dest="picard_ref", help='The picard reference genome to use or index' )
    parser.add_option( '-i', '--inputDir', dest="inputDir", help='The path to fastq files to use for the mapping' )
    parser.add_option( '-b', '--bait-file', dest="bait_file", help='The bait file' )
    parser.add_option( '', '--target-file', dest="target_file", help='The target file' )
    parser.add_option( '', '--refseq-file', dest="refseq_file", help='The refseq file' )
    parser.add_option( '', '--gene-targets', dest="gene_targetfile", help="The gene target file")
    parser.add_option( '', '--vendor', dest="vendor", help="the vendor used for exome capture" )
    parser.add_option( '--output-bwa-bam', dest="output_bwa_bam", help="the file to save")
    parser.add_option( '--output-picard-markDuplicates-html', dest="output_picard_markDuplicates_html", help='The file to save' )
    parser.add_option( '--output-picard-markDuplicates-directory', dest="output_picard_markDuplicates_directory", help='The dir output' )
    parser.add_option( '--output-picard-shortReadGCBias', dest="output_picard_shortReadGCBias", help='The file to save' )
    parser.add_option( '--output-picard-shortReadGCBias-directory', dest="output_picard_shortReadGCBias_directory", help='The dir output' )
    parser.add_option( '--output-picard-insertSize', dest="output_picard_insertSize", help='The file to save' )
    parser.add_option( '--output-picard-insertSize-directory', dest="output_picard_insertSize_directory", help='The dir output' )
    parser.add_option( '--output-picard-alignmentSummary', dest="output_picard_alignmentSummary", help='The file to save' )
    parser.add_option( '--output-picard-alignmentSummary-directory', dest="output_picard_alignmentSummary_directory", help='The dir output' )
    parser.add_option( '--output-picard-HSMetrics', dest="output_picard_HSMetrics", help='The file to save' )
    parser.add_option( '--output-picard-HSMetrics-directory', dest="output_picard_HSMetrics_directory", help='The file output' )
    parser.add_option( '--output-gatk-depthOfCoverage-geneSummary', dest="output_gatk_depthOfCoverage_geneSummary", help='The file output' )
    parser.add_option( '--output-gatk-depthOfCoverage-intervalSummary', dest="output_gatk_depthOfCoverage_intervalSummary", help='The file output' )
    parser.add_option( '--output-gatk-diagnoseTargets-vcf', dest="output_gatk_diagnoseTargets_vcf", help='The file to save' )
    parser.add_option( '--output-log', dest="output_log", help='The file output' )
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
    parser.add_option( '--username', dest="username", help='Globus username for file transfer' )
    parser.add_option( '--goauth-token', dest="goauth_token", help='Globus token' )
    parser.add_option( '--source-ep', dest="source_ep", help='Globus source endpoint' )
    parser.add_option( '--destination-ep', dest="destination_ep", help='Globus destination endpoint' )
    parser.add_option( '--destination-dir', dest="destination_dir", help='Globus destination dir' )
    parser.add_option( '--deadline', dest="deadline", help='Globus transfer deadline' )
    (options, args) = parser.parse_args()

    final_outputs = [ options.output_bwa_bam, options.output_picard_markDuplicates_html, options.output_picard_markDuplicates_directory, options.output_picard_shortReadGCBias, options.output_picard_shortReadGCBias_directory, options.output_picard_insertSize, options.output_picard_insertSize_directory, options.output_picard_alignmentSummary, options.output_picard_alignmentSummary_directory, options.output_picard_HSMetrics, options.output_picard_HSMetrics_directory, options.output_gatk_depthOfCoverage_geneSummary, options.output_gatk_depthOfCoverage_intervalSummary, options.output_gatk_diagnoseTargets_vcf ]

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
            #delete_file(step_input_files, keep_files, step_id, final_outputs)


    #clean up
    #shutil.rmtree(output_dir)

if __name__ == "__main__":
    __main__()
