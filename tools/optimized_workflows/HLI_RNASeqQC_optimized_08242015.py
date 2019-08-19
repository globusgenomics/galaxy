# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""
Runs the HLI RNA-SEQ QC pipeline in optimized mode in one Node
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

    # get ncores
    ncpu = get_ncores()
    nthreads = ncpu/2

    # get memory
    jvm_heap = get_linuxRAM()

    # get the reference datasets
    picard_reference_db = options.picard_ref
    tophat2_reference_db = options.tophat2_ref
    refGTF= options.gtf_file
    refFlat = "/mnt/galaxyIndices/genomes/Hsapiens/hg38/annotation/Homo_sapiens.GRCh38.77.refFlat"
    #JAVA_JAR_PATH = os.environ['JAVA_JAR_PATH']
    JAVA_JAR_PATH = "/mnt/galaxyTools/tools/picard/1.134"
    qualimap_jar = "/mnt/galaxyTools/tools/qualimap/qualimap_v2.1.1"

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

    input_files = [refGTF]
    outputF = "%s/refFlat.txt" % (output_dir)
    tmp_out = "%s/refFlat.tmp.txt" % (output_dir)
    output_files = [outputF]
    cmd = "gtfToGenePred -genePredExt -geneNameAsName2 %s %s; paste <(cut -f 12 %s) <(cut -f 1-10 %s) > %s; rm -rf %s" % (refGTF, tmp_out, tmp_out, tmp_out, outputF, tmp_out)
    command_meta[step].append({"cl":cmd, "input_files": input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 1
    ncpu_step = int(ncpu)
    command_meta[step] = []
    input1 = command_meta[0][0]["output_files"][0]
    input2 = command_meta[0][0]["output_files"][1]
    input_files = [input1, input2]
    outputF = "%s/%s-bwa-out.bam" % (output_dir, sample_name)
    outputFrest = options.output_tophat_bam
    outputJunctionBed = "%s/%s-junctions.bed" % (output_dir, sample_name)
    #outputLog = "%s/%s-tophat.log" % (output_dir, sample_name)
    outputLog = options.output_tophat_log
    output_files = [outputF, outputJunctionBed, outputLog, outputFrest]
    cmd = " python /opt/galaxy/tools/ngs_rna/tophat2_wrapper.py --num-threads=%s --junctions-output=%s --hits-output=%s --indexes-path=%s --single-paired=paired --input1=%s --input2=%s -r 300 --mate-std-dev=20 --settings=full --read-mismatches 2 --read-edit-dist 2 --read-realign-edit-dist 1000 -a 8 -m 0 -i 70 -I 500000 -g 20 --min-segment-intron 50 --max-segment-intron 500000 --seg-mismatches=2 --seg-length=25 --library-type=fr-unstranded --max-insertion-length 3 --max-deletion-length 3 -G %s --no-coverage-search; cp ./tophat_out/align_summary.txt %s; cp %s %s" % (ncpu_step, outputJunctionBed, outputF, tophat2_reference_db, input1, input2, refGTF, outputLog, outputF, outputFrest)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files": output_files})

    #
    step = 2
    command_meta[step] = []
    inputF = command_meta[1][0]["output_file"]
    #outputF = "%s/%s-tophat-out.rd.picard.alignmentSummaryMetrics" % (output_dir, sample_name)
    outputF = options.output_picard_alignmentSummary
    outputMetricsDir = options.output_picard_alignmentSummary_directory
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/picard/picard_wrapper.py -i %s -d %s -t %s --assumesorted true -b false --adaptors \"\" --maxinsert 100000 -n \"Picard Alignment Summary Metrics\" --datatype bam -j \"%s/picard.jar CollectAlignmentSummaryMetrics\"  --tmpdir %s --ref %s" % (inputF, outputMetricsDir, outputF, JAVA_JAR_PATH, tmp_dir, picard_reference_db)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 3
    command_meta[step] = []
    inputF = command_meta[1][0]["output_file"]
    outputF = "%s/%s-tophat-out.rd.bam" % (output_dir, sample_name)
    outputMetrics = options.output_picard_markDuplicates_html
    outputMetricsDir = options.output_picard_markDuplicates_directory
    #outputMetrics = "%s/%s-SM1.dups" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF, outputMetrics]
    cmd = "python /opt/galaxy/tools/picard/picard_wrapper.py --maxjheap %sm -i %s -n Dupes_Marked --tmpdir %s -o %s --assumesorted true --readregex \"[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*\" --optdupdist 100 -j \"%s/picard.jar MarkDuplicates\" -d %s -t %s -e bam" % (jvm_heap, inputF, tmp_dir, outputF, JAVA_JAR_PATH, outputMetricsDir, outputMetrics)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    #step = 4
    #command_meta[step] = []
    inputF = command_meta[1][0]["output_file"]
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
    inputF = command_meta[1][0]["output_file"]
    #outputF = "%s/%s-tophat-out.rd.picard.collectGcBiasMetrics" % (output_dir, sample_name)
    outputF = options.output_picard_shortReadGCBias
    outputMetricsDir = options.output_picard_shortReadGCBias_directory
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/picard/picard_wrapper.py -i %s -d %s -t %s --windowsize 100 --mingenomefrac 1e-05 -n \"Short Read GC Bias Metrics\" --tmpdir %s -j \"%s/picard.jar CollectGcBiasMetrics\"  --ref %s" % (inputF, outputMetricsDir, outputF, tmp_dir, JAVA_JAR_PATH, picard_reference_db)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 5
    command_meta[step] = []
    inputF = command_meta[1][0]["output_file"]
    ribosomalRef = "/mnt/galaxyIndices/genomes/Hsapiens/hg38/annotation/GRCh38.ribosomal.tabular" 
    #outputF = "%s/%s-tophat-out.rd.picard.CollectRnaSeqMetrics_log" % (output_dir, sample_name)
    outputF = options.output_picard_collectRnaSeqSummary
    outputMetricsDir = options.output_picard_collectRnaSeqSummary_directory
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/picard/picard_wrapper.py -i %s -d %s -t %s -n \"CollectRnaSeqMetrics\" -j \"%s/picard.jar CollectRnaSeqMetrics\" --tmpdir %s --ref %s --ref_flat %s --ribosomalintervals %s --strandspecificity NONE --minimumlength 500 --rrnafragmentpercentage 0.8 --metricaccumulationlevel ALL_READS --assumesorted true" % (inputF, outputMetricsDir, outputF, JAVA_JAR_PATH, tmp_dir, picard_reference_db, refFlat, ribosomalRef)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 6
    command_meta[step] = []
    inputF = command_meta[1][0]["output_file"]
    input1 = command_meta[0][0]["output_files"][0]
    #outputF = "%s/%s-tophat-out.rd.tophat_stats" % (output_dir, sample_name)
    outputF = options.output_tophatStats
    input_files = [inputF, input1]
    output_files = [outputF]
    cmd = "perl /opt/galaxy/tools/ngs_rna/tophatstatsPE.pl %s %s > %s" % (inputF, input1, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    #step = 8
    #command_meta[step] = []
#    inputF = command_meta[1][0]["output_file"]
#    #outputF = "%s/qualimapReport.html" % (output_dir, sample_name)
#    outputF = options.output_qualimap
#    input_files = [inputF]
#    output_files = [outputF]
#    cmd = "java -Xms32m -Xmx12000M -classpath %s/qualimap.jar:%s/lib/* org.bioinfo.ngs.qc.qualimap.main.NgsSmartMain bamqc --outdir ./qualimap_tmp/ -bam %s -hm 3 -nr 1000 -nt 8  -outformat HTML -c; cp ./qualimap_tmp/qualimapReport.html %s" % (qualimap_jar, qualimap_jar, inputF, outputF)
#    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    #step = 9
    #command_meta[step] = []
    inputF = command_meta[1][0]["output_file"]
    #outputF = "%s/%s-rseqc.output.qual.boxplot.pdf" % (output_dir, sample_name)
    outputF = options.output_rseqc_readQuality
    outputR = "%s/%s-rseqc.output.qual.r" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF, outputR]
    cmd = "read_quality.py -i %s -o output -r 1000; cp %s %s" % (inputF, "output.qual.boxplot.pdf", outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    #step = 10
    #command_meta[step] = []
    inputF = command_meta[1][0]["output_file"]
    #outputF = "%s/%s-rseqc.output.NVC_plot.pdf" % (output_dir, sample_name)
    outputF = options.output_rseqc_readNVC
    outputR = "%s/%s-rseqc.output.NVC_plot.r" % (output_dir, sample_name)
    outputX = "%s/%s-rseqc.output.NVC.xls" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF, outputR, outputX]
    cmd = "read_NVC.py -i %s -o output; cp %s %s" % (inputF, "output.NVC_plot.pdf", outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    #step = 11
    #command_meta[step] = []
    inputF = command_meta[1][0]["output_file"]
    #outputF = "%s/%s-rseqc.output.GC_plot.pdf" % (output_dir, sample_name)
    outputF = options.output_rseqc_readGC
    outputR = "%s/%s-rseqc.output.GC_plot.r" % (output_dir, sample_name)
    outputX = "%s/%s-rseqc.output.GC.xls" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF, outputR, outputX]
    cmd = "read_GC.py -i %s -o output; cp %s %s" % (inputF, "output.GC_plot.pdf", outputF)
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
    parser.add_option( '--tophat2-ref', dest="tophat2_ref", help='The tophat2 reference genome to use or index' )
    parser.add_option( '--picard-ref', dest="picard_ref", help='The picard reference genome to use or index' )
    parser.add_option( '-i', '--inputDir', dest="inputDir", help='The path to fastq files to use for the mapping' )
    parser.add_option( '-g', '--gtf', dest="gtf_file", help='The gtf file' )
    parser.add_option( '--output-tophat-log', dest="output_tophat_log", help="The file to save the output log" )
    parser.add_option( '--output-tophat-bam', dest="output_tophat_bam", help="The file to save the output bam" )
    parser.add_option( '--output-picard-markDuplicates-html', dest="output_picard_markDuplicates_html", help='The file to save' )
    parser.add_option( '--output-picard-markDuplicates-directory', dest="output_picard_markDuplicates_directory", help='The dir output' )
    parser.add_option( '--output-picard-shortReadGCBias', dest="output_picard_shortReadGCBias", help='The file to save' )
    parser.add_option( '--output-picard-shortReadGCBias-directory', dest="output_picard_shortReadGCBias_directory", help='The dir output' )
    parser.add_option( '--output-picard-insertSize', dest="output_picard_insertSize", help='The file to save' )
    parser.add_option( '--output-picard-insertSize-directory', dest="output_picard_insertSize_directory", help='The dir output' )
    parser.add_option( '--output-picard-alignmentSummary', dest="output_picard_alignmentSummary", help='The file to save' )
    parser.add_option( '--output-picard-alignmentSummary-directory', dest="output_picard_alignmentSummary_directory", help='The dir output' )
    parser.add_option( '--output-picard-collectRnaSeqSummary', dest="output_picard_collectRnaSeqSummary", help='The file to save' )
    parser.add_option( '--output-picard-collectRnaSeqSummary-directory', dest="output_picard_collectRnaSeqSummary_directory", help='The file output' )
    parser.add_option( '--output-tophatStats', dest="output_tophatStats", help='The file to save' )
    parser.add_option( '--output-rseqc-readQuality', dest="output_rseqc_readQuality", help='The file output' )
    parser.add_option( '--output-rseqc-readNVC', dest="output_rseqc_readNVC", help='The file output' )
    parser.add_option( '--output-rseqc-readGC', dest="output_rseqc_readGC", help='The file output' )
    parser.add_option( '--output-rseqc-junctionSaturation', dest="output_rseqc_junctionSaturation", help='The file output' )
    parser.add_option( '--output-rseqc-junctionAnnotation-events', dest="output_rseqc_junctionAnnotation_events", help='The file output' )
    parser.add_option( '--output-rseqc-junctionAnnotation-junction', dest="output_rseqc_junctionAnnotation_junction", help='The file output' )
    parser.add_option( '--output-rseqc-readDuplication', dest="output_rseqc_readDuplication", help='The file output' )
    parser.add_option( '--output-qualimap', dest="output_qualimap", help='The file output' )
    parser.add_option( '--output-log', dest="output_log", help='The file output' )
    parser.add_option( '--username', dest="username", help='Globus username for file transfer' )
    parser.add_option( '--goauth-token', dest="goauth_token", help='Globus token' )
    parser.add_option( '--source-ep', dest="source_ep", help='Globus source endpoint' )
    parser.add_option( '--destination-ep', dest="destination_ep", help='Globus destination endpoint' )
    parser.add_option( '--destination-dir', dest="destination_dir", help='Globus destination dir' )
    parser.add_option( '--deadline', dest="deadline", help='Globus transfer deadline' )
    (options, args) = parser.parse_args()

    final_outputs = [options.output_tophat_log, options.output_picard_markDuplicates_html, options.output_picard_markDuplicates_directory, options.output_picard_shortReadGCBias, options.output_picard_shortReadGCBias_directory, options.output_picard_insertSize, options.output_picard_insertSize_directory, options.output_picard_alignmentSummary, options.output_picard_alignmentSummary_directory, options.output_picard_collectRnaSeqSummary, options.output_picard_collectRnaSeqSummary_directory, options.output_tophatStats, options.output_rseqc_readQuality, options.output_rseqc_readNVC, options.output_rseqc_readGC, options.output_rseqc_junctionSaturation, options.output_rseqc_junctionAnnotation_events, options.output_rseqc_junctionAnnotation_junction ]

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
