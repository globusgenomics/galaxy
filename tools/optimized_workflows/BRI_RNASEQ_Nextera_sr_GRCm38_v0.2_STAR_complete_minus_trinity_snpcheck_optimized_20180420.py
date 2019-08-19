# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""
Runs the BRI RNA-seq nextera pipeline in optimized mode in one Node
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
    #print cmd
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
    sample_name = options.sample_name

    #lcm = 1
    lcm = 4

    # get ncores
    ncpu = get_ncores()/lcm
    nthreads = ncpu/2

    # get memory
    jvm_heap = get_linuxRAM()/lcm

    # get the reference datasets
    picard_reference_db = options.picard_ref
    star_reference_db = options.star_ref
    #JAVA_JAR_PATH = "/mnt/galaxyTools/tools/picard/1.134"
    JAVA_JAR_PATH = "/mnt/galaxyTools/tools/picard/1.56"
    fastqc_path = "/mnt/galaxyTools/tools/FastQC/0.11.3/"
    lane_1_fastq = options.fastq_lane_1
    lane_2_fastq = options.fastq_lane_2
    lane_3_fastq = options.fastq_lane_3
    lane_4_fastq = options.fastq_lane_4
    lane_5_fastq = options.fastq_lane_5
    lane_6_fastq = options.fastq_lane_6
    lane_7_fastq = options.fastq_lane_7
    lane_8_fastq = options.fastq_lane_8

    ## Setup workflow
    command_meta = {}


    #
    step = 0
    command_meta[step] = []
    input_files = [lane_1_fastq, lane_2_fastq, lane_3_fastq, lane_4_fastq, lane_5_fastq, lane_6_fastq, lane_7_fastq, lane_8_fastq]
    outputF1 = "%s/%s-R1.fastq.gz" % (output_dir, sample_name)
    outputFrest = options.output_contatenate_fastq
    output_files = [outputF1, outputFrest]
    cmd = "python /opt/galaxy/tools/filters/catWrapper.py %s %s; cp %s %s" % (outputF1, " ".join(input_files), outputF1, outputFrest)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF1, "output_files": output_files})

    #
    step = 1
    ncpu_step = int(ncpu)
    command_meta[step] = []
    input1 = command_meta[step-1][0]["output_file"]
    input_files = [input1]
    outputF = "%s/%s-R1.fastq" % (output_dir, sample_name)
    output_files = [outputF]
    cmd = "zcat %s > %s" % (input1, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files": output_files})

    #
    step = 2
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-R1.trim.fastq" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/fastq/fastq_trimmer.py %s %s 0 1 offsets_absolute sanger exclude_zero_length" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 3
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-R1.trim.fastqmcf.fastq" % (output_dir, sample_name)
    outputLog = "%s/%s-R1.trim.fastqmcf.log" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF, outputLog]
    cmd = "fastq-mcf -l 50 -q 0 -C 600000 -o %s %s %s > %s" % (outputF, options.adapters_file, inputF, outputLog)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 4
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputHtml = "%s/%s-R1.trim.fastq.fastqc.html" % (output_dir, sample_name)
    outputHtmlRest = options.output_fastqc_html
    outputDir = "%s/%s-R1.trim.fastq.fastqc" % (output_dir, sample_name)
    outputText = "%s/%s-R1.trim.fastq.fastqc.txt" % (output_dir, sample_name)
    outputTextRest = options.output_fastqc_text
    input_files = [inputF]
    output_files = [outputHtml, outputText]
    cmd = "python /opt/galaxy/tools/rgenetics/rgFastQC.py -i %s -d %s -o %s -t %s -f fastqsanger -j FASTQC_OUT -e %s/fastqc; cp %s %s; cp %s %s" % (inputF, outputDir, outputHtml, outputText, fastqc_path, outputHtml, outputHtmlRest, outputText, outputTextRest)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    #step = 4
    #command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-R1.trim.qc.fastq" % (output_dir, sample_name)
    outputFrest = options.output_qc_fastq
    input_files = [inputF]
    output_files = [outputF, outputFrest]
    cmd = "python /opt/galaxy/tools/fastq/fastq_trimmer_by_quality.py %s %s -f sanger -s 1 -t 1 -e 53 -a min -x 0 -c '>=' -q 30.0; cp %s %s" % (inputF, outputF, outputF, outputFrest)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 5
    ncpu_step = int(ncpu)
    command_meta[step] = []
    input1 = command_meta[step-1][1]["output_files"][0]
    input_files = [input1]
    outputF = "%s/Aligned.sortedByCoord.out.bam" % (output_dir)
    outputFindex = "%s/Aligned.sortedByCoord.out.bam.bai" % (output_dir)
    outputFrest = options.output_star_bam
    output_files = [outputF, outputFindex]
    cmd = "STAR --genomeLoad NoSharedMemory --genomeDir %s --readFilesIn %s --sjdbGTFfile %s --runThreadN %s --outSAMtype BAM SortedByCoordinate --outSAMunmapped None --outSAMstrandField intronMotif; cp %s %s; samtools index %s" % (star_reference_db, input1, options.gtf_file, ncpu_step, outputF, outputFrest, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files": output_files})

    #
    step = 6
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-star-out.rd.bam" % (output_dir, sample_name)
    outputMetrics = options.output_picard_markDuplicates_html
    outputMetricsDir = options.output_picard_markDuplicates_directory
    input_files = [inputF]
    output_files = [outputF, outputMetrics]
    cmd = "python /opt/galaxy/tools/picard-pre-1.128/picard_wrapper.py --maxjheap %sm -i %s -n Dupes_Marked --tmpdir %s -o %s --assumesorted true --readregex \"[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*\" --remdups true --optdupdist 100 -j \"%s/MarkDuplicates.jar\" -d %s -t %s -e bam" % (int(jvm_heap/2), inputF, tmp_dir, outputF, JAVA_JAR_PATH, outputMetricsDir, outputMetrics)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    #step = 6
    #command_meta[step] = []
    inputF = command_meta[step-1][0]["output_files"][0]
    inputFindex = command_meta[step-1][0]["output_files"][1]
    outputF = "%s/%s-star-out.mpileup" % (output_dir, sample_name)
    outputFstdout = "%s/%s-star-out.mpileup.stdout" % (output_dir, sample_name)
    input_files = [inputF, inputFindex]
    output_files = [outputF, outputFstdout]
    cmd = "python /opt/galaxy/tools/samtools/samtools_wrapper.py -p 'samtools mpileup' --stdout %s -p '-f %s' -d \" \" \"%s\" \"bam\" \"bam_input_0\" -d \"\" \"%s\" \"bam_index\" \"bam_input_0\" -p ' -B -C 50 -d 250 -l %s -q 30 -Q 13 -g -e 20 -h 100 -L 250 -o 40 > %s'" % (outputFstdout, picard_reference_db, inputF, inputFindex, options.bed_file, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 7
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-star-out.reorder.bam" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/picard-pre-1.128/picard_wrapper.py --maxjheap %sm --input=%s --ref=%s --allow-inc-dict-concord=false --allow-contig-len-discord=false --output-format=bam --output=%s --tmpdir %s -j \"%s/ReorderSam.jar\"" % (int(jvm_heap/2), inputF, picard_reference_db, outputF, tmp_dir, JAVA_JAR_PATH)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    #step = 7
    #command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-star-out.sam" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/samtools/bam_to_sam.py --input1=%s --output1=%s --header" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    #step = 8
    #command_meta[step] = []
    inputF = command_meta[step-1][1]["output_files"][0]
    outputF = "%s/%s-star-out.snp.vcf" % (output_dir, sample_name)
    outputFrest = options.output_snpcheck_vcf
    outputFstdout = "%s/%s-star-out.snp.stdout" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF, outputFrest, outputFstdout]
    cmd = "bcftools call -O v -v -c -o %s %s 2> %s; cp %s %s" % (outputF, inputF, outputFstdout, outputF, outputFrest)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 8
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = options.output_picard_collectRnaSeqSummary
    outputMetricsDir = options.output_picard_collectRnaSeqSummary_directory
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/picard-pre-1.128/picard_wrapper.py -i %s -d %s -t %s -n \"Dupes Marked\" -j \"%s/CollectRnaSeqMetrics.jar\" --tmpdir %s --ref %s --ref_flat %s --ribosomalintervals %s --strandspecificity \"FIRST_READ_TRANSCRIPTION_STRAND\" --minimumlength 500 --rrnafragmentpercentage 0.8 --metricaccumulationlevel ALL_READS --assumesorted true" % (inputF, outputMetricsDir, outputF, JAVA_JAR_PATH, tmp_dir, picard_reference_db, options.refFlat, options.ribosomal_file)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    #step = 8
    #command_meta[step] = []
    inputF = command_meta[step-1][1]["output_file"]
    inputFastq = command_meta[1][0]["output_file"]
    outputF = options.output_tophatStats
    input_files = [inputF]
    output_files = [outputF]
    cmd = "perl /opt/galaxy/tools/ngs_rna/tophatstatsPE_bri.pl %s %s > %s" % (inputF, inputFastq, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    #step = 8
    #command_meta[step] = []
    inputF = command_meta[step-1][1]["output_file"]
    outputF = "%s/%s-star-out.sorted.sam" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/filters/sorter.old.bri.py --input=%s --out_file1=%s --column=1 --style=alpha --order=DESC" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    #step = 8
    #command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    #outputF = "%s/%s-tophat-out.rd.picard.alignmentSummaryMetrics" % (output_dir, sample_name)
    outputF = options.output_picard_alignmentSummary
    outputMetricsDir = options.output_picard_alignmentSummary_directory
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/picard-pre-1.128/picard_wrapper.py -i %s -d %s -t %s --assumesorted true -b false --adaptors \"\" --maxinsert 100000 -n \"Picard Alignment Summary Metrics\" --datatype bam -j \"%s/CollectAlignmentSummaryMetrics.jar\" --tmpdir %s --ref %s" % (inputF, outputMetricsDir, outputF, JAVA_JAR_PATH, tmp_dir, picard_reference_db)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 9
    command_meta[step] = []
    inputF = command_meta[step-1][2]["output_file"]
    outputF1 = "%s/%s-htseq-count.txt" % (output_dir, sample_name)
    outputF2 = "%s/%s-htseq-count-log.txt" % (output_dir, sample_name)
    outputF1Rest = options.output_htseq_counts
    outputF2Rest = options.output_htseq_log
    input_files = [inputF]
    output_files = [outputF]
    cmd = "htseq-count -q --mode=union --stranded=no --minaqual=0 --type=exon --idattr=gene_id %s %s | awk \'{if ($1 ~ \"no_feature|ambiguous|too_low_aQual|not_aligned|alignment_not_unique\") print $0 | \"cat 1>&2\"; else print $0}\' > temp.out.txt 2> %s && mv temp.out.txt %s; cp %s %s; cp %s %s" % (inputF, options.gtf_file, outputF2, outputF1, outputF1, outputF1Rest, outputF2, outputF2Rest )
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
    parser.add_option( '--star-ref', dest="star_ref", help='The star reference genome to use or index' )
    parser.add_option( '--picard-ref', dest="picard_ref", help='The picard reference genome to use or index' )
    parser.add_option( '', '--lane1', dest="fastq_lane_1", help='The path to fastq files to use for the mapping' )
    parser.add_option( '', '--lane2', dest="fastq_lane_2", help='The path to fastq files to use for the mapping' )
    parser.add_option( '', '--lane3', dest="fastq_lane_3", help='The path to fastq files to use for the mapping' )
    parser.add_option( '', '--lane4', dest="fastq_lane_4", help='The path to fastq files to use for the mapping' )
    parser.add_option( '', '--lane5', dest="fastq_lane_5", help='The path to fastq files to use for the mapping' )
    parser.add_option( '', '--lane6', dest="fastq_lane_6", help='The path to fastq files to use for the mapping' )
    parser.add_option( '', '--lane7', dest="fastq_lane_7", help='The path to fastq files to use for the mapping' )
    parser.add_option( '', '--lane8', dest="fastq_lane_8", help='The path to fastq files to use for the mapping' )
    parser.add_option( '', '--bed-file', dest="bed_file", help="bed file" )
    parser.add_option( '', '--gtf', dest="gtf_file", help="gtf file" )
    parser.add_option( '', '--refflat-file', dest="refFlat", help="refflat file" )
    parser.add_option( '', '--ribosomal-file', dest="ribosomal_file", help="ribosomal file" )
    parser.add_option( '', '--adapters-file', dest="adapters_file", help="adapters file" )
    parser.add_option( '--output-concatenated-fastq', dest="output_contatenate_fastq", help="the file to save")
    parser.add_option( '--output-qc-fastq', dest="output_qc_fastq", help="the file to save")
    parser.add_option( '--output-star-bam', dest="output_star_bam", help="the file to save")
    parser.add_option( '--output-picard-markDuplicates-html', dest="output_picard_markDuplicates_html", help='The file to save' )
    parser.add_option( '--output-picard-markDuplicates-directory', dest="output_picard_markDuplicates_directory", help='The dir output' )
    parser.add_option( '--output-picard-alignmentSummary', dest="output_picard_alignmentSummary", help='The file to save' )
    parser.add_option( '--output-picard-alignmentSummary-directory', dest="output_picard_alignmentSummary_directory", help='The dir output' )
    parser.add_option( '--output-picard-collectRnaSeqSummary', dest="output_picard_collectRnaSeqSummary", help='The file to save' )
    parser.add_option( '--output-picard-collectRnaSeqSummary-directory', dest="output_picard_collectRnaSeqSummary_directory", help='The file output' )
    parser.add_option( '--output-tophatStats', dest="output_tophatStats", help='The file to save' )
    parser.add_option( '--output-htseq-counts', dest="output_htseq_counts", help="File to save" )
    parser.add_option( '--output-htseq-log', dest="output_htseq_log", help="File to save" )
    parser.add_option( '--output-fastqc-html', dest="output_fastqc_html", help="File to save" )
    parser.add_option( '--output-fastqc-text', dest="output_fastqc_text", help="File to save" )
    parser.add_option( '--output-snpcheck-vcf', dest="output_snpcheck_vcf", help="File to save" )
    parser.add_option( '--output-log', dest="output_log", help='The file output' )
    parser.add_option( '--username', dest="username", help='Globus username for file transfer' )
    parser.add_option( '--goauth-token', dest="goauth_token", help='Globus token' )
    parser.add_option( '--source-ep', dest="source_ep", help='Globus source endpoint' )
    parser.add_option( '--destination-ep', dest="destination_ep", help='Globus destination endpoint' )
    parser.add_option( '--destination-dir', dest="destination_dir", help='Globus destination dir' )
    parser.add_option( '--deadline', dest="deadline", help='Globus transfer deadline' )
    (options, args) = parser.parse_args()

    final_outputs = [options.output_qc_fastq, options.output_contatenate_fastq, options.output_star_bam, options.output_picard_markDuplicates_html, options.output_picard_markDuplicates_directory, options.output_picard_alignmentSummary, options.output_picard_alignmentSummary_directory, options.output_picard_collectRnaSeqSummary, options.output_picard_collectRnaSeqSummary_directory, options.output_tophatStats, options.output_htseq_counts, options.output_htseq_log, options.output_fastqc_html, options.output_fastqc_text, options.output_snpcheck_vcf, options.output_log ]


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
