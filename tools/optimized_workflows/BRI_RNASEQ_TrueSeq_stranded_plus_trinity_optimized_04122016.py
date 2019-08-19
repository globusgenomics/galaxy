# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""
Runs the BRI RNA-seq trueseq pipeline in optimized mode in one Node
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
    tophat_reference_db = options.tophat_ref
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
    #step = 2
    #command_meta[step] = []
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
    step = 3
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-R1.trim.qc.fastq" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/fastq/fastq_trimmer_by_quality.py %s %s -f sanger -s 1 -t 1 -e 53 -a min -x 0 -c '>=' -q 30.0" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 4
    ncpu_step = int(ncpu)
    command_meta[step] = []
    input1 = command_meta[step-1][0]["output_files"][0]
    input_files = [input1]
    outputF = "%s/tophat_out/accepted_hits.bam" % output_dir
    #outputF = "%s/%s-tophat-out.bam" % (output_dir, sample_name)
    outputJunctionBed = "%s/%s-junctions.bed" % (output_dir, sample_name)
    outputFrest = options.output_tophat_bam
    output_files = [outputF, outputJunctionBed]
    cmd = "python /opt/galaxy/tools/ngs_rna/tophat_wrapper.py --num-threads=%s --junctions-output=%s --hits-output=%s --indexes-path=%s --single-paired=single --input1=%s --settings=full -a 8 -m 0 -i 70 -I 500000 -g 20 --min-segment-intron 50 --max-segment-intron 500000 --initial-read-mismatches=2 --seg-mismatches=2 --seg-length=25 --library-type=fr-firststrand --max-insertion-length 3 --max-deletion-length 3 --no-closure-search --coverage-search --min-coverage-intron 50 --max-coverage-intron 20000; cp %s %s" % (ncpu_step, outputJunctionBed, outputF, tophat_reference_db, input1, outputF, outputFrest)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files": output_files})

    #
    step = 5
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-tophat-out.rd.bam" % (output_dir, sample_name)
    outputMetrics = options.output_picard_markDuplicates_html
    outputMetricsDir = options.output_picard_markDuplicates_directory
    input_files = [inputF]
    output_files = [outputF, outputMetrics]
    cmd = "python /opt/galaxy/tools/picard-pre-1.128/picard_wrapper.py --maxjheap %sm -i %s -n Dupes_Marked --tmpdir %s -o %s --assumesorted true --readregex \"[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*\" --remdups true --optdupdist 100 -j \"%s/MarkDuplicates.jar\" -d %s -t %s -e bam" % (int(jvm_heap/2), inputF, tmp_dir, outputF, JAVA_JAR_PATH, outputMetricsDir, outputMetrics)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    #step = 5
    #command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-tophat-out.reorder.bam" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/picard-pre-1.128/picard_wrapper.py --maxjheap %sm --input=%s --ref=%s --allow-inc-dict-concord=false --allow-contig-len-discord=false --output-format=bam --output=%s --tmpdir %s -j \"%s/ReorderSam.jar\"" % (int(jvm_heap/2), inputF, picard_reference_db, outputF, tmp_dir, JAVA_JAR_PATH)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    #step = 5
    #command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-tophat-out.sam" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/samtools/bam_to_sam.py --input1=%s --output1=%s --header" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 6
    command_meta[step] = []
    inputF = command_meta[step-1][1]["output_file"]
    #outputF = "%s/%s-tophat-out.rd.picard.CollectRnaSeqMetrics_log" % (output_dir, sample_name)
    outputF = options.output_picard_collectRnaSeqSummary
    outputMetricsDir = options.output_picard_collectRnaSeqSummary_directory
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/picard-pre-1.128/picard_wrapper.py -i %s -d %s -t %s -n \"Dupes Marked\" -j \"%s/CollectRnaSeqMetrics.jar\" --tmpdir %s --ref %s --ref_flat %s --ribosomalintervals %s --strandspecificity \"FIRST_READ_TRANSCRIPTION_STRAND\" --minimumlength 500 --rrnafragmentpercentage 0.8 --metricaccumulationlevel ALL_READS --assumesorted true" % (inputF, outputMetricsDir, outputF, JAVA_JAR_PATH, tmp_dir, picard_reference_db, options.refFlat, options.ribosomal_file)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    #step = 6
    ncpu_step = int(ncpu)
    #command_meta[step] = []
    input1 = command_meta[step-3][0]["output_files"][0]
    input_files = [input1]
    outputF = "%s/trinity/Trinity.fasta" % (output_dir)
    outputFrest = options.output_trinity_transcripts
    outputLog = "%s/%s-trinity-log.txt" % (output_dir, sample_name)
    output_files = [outputF, outputLog]
    cmd = "python /opt/galaxy/tools/trinity/trinityToolWrapper.py Trinity --output %s/trinity --max_memory %sG --CPU %s --single %s --seqType fq --min_kmer_cov 1 --max_reads_per_graph 20000000 --bflyHeapSpaceMax %sG > %s; touch %s; cp %s %s" % (output_dir, int(jvm_heap/2000), ncpu_step, input1, int(jvm_heap/2000), outputLog, outputF, outputF, outputFrest )
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files": output_files})

    #
    #step = 6
    #command_meta[step] = []
    inputF = command_meta[step-1][2]["output_file"]
    inputFastq = command_meta[1][0]["output_file"]
    outputF = options.output_tophatStats
    input_files = [inputF]
    output_files = [outputF]
    cmd = "perl /opt/galaxy/tools/ngs_rna/tophatstatsPE_bri.pl %s %s > %s" % (inputF, inputFastq, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    #step = 6
    #command_meta[step] = []
    inputF = command_meta[step-1][2]["output_file"]
    outputF = "%s/%s-tophat-out.sorted.sam" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/filters/sorter.old.bri.py --input=%s --out_file1=%s --column=1 --style=alpha --order=DESC" % (inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    #step = 6
    #command_meta[step] = []
    inputF = command_meta[step-1][1]["output_file"]
    #outputF = "%s/%s-tophat-out.rd.picard.alignmentSummaryMetrics" % (output_dir, sample_name)
    outputF = options.output_picard_alignmentSummary
    outputMetricsDir = options.output_picard_alignmentSummary_directory
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/picard-pre-1.128/picard_wrapper.py -i %s -d %s -t %s --assumesorted true -b false --adaptors \"\" --maxinsert 100000 -n \"Picard Alignment Summary Metrics\" --datatype bam -j \"%s/CollectAlignmentSummaryMetrics.jar\"  --tmpdir %s --ref %s" % (inputF, outputMetricsDir, outputF, JAVA_JAR_PATH, tmp_dir, picard_reference_db)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 7
    command_meta[step] = []
    inputF = command_meta[step-1][3]["output_file"]
    outputF1 = "%s/%s-htseq-count.txt" % (output_dir, sample_name)
    outputF2 = "%s/%s-htseq-count-log.txt" % (output_dir, sample_name)
    outputF1Rest = options.output_htseq_counts
    outputF2Rest = options.output_htseq_log
    input_files = [inputF]
    output_files = [outputF]
    cmd = "htseq-count -q --mode=union --stranded=reverse --minaqual=0 --type=exon --idattr=gene_id %s %s | awk \'{if ($1 ~ \"no_feature|ambiguous|too_low_aQual|not_aligned|alignment_not_unique\") print $0 | \"cat 1>&2\"; else print $0}\' > temp.out.txt 2> %s && mv temp.out.txt %s; cp %s %s; cp %s %s" % (inputF, options.gtf_file, outputF2, outputF1, outputF1, outputF1Rest, outputF2, outputF2Rest )
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    return command_meta

def __main__():
    descr = "snpir_optimized_workflow.py: Runs the Meissnert SNPiR RNA-variant discovery pipeline in optimized mode in one Node"
    parser = optparse.OptionParser(description=descr)
    parser.add_option( '--sample-name', dest="sample_name", help='sample name' )
    parser.add_option( '--tophat-ref', dest="tophat_ref", help='The tophat reference genome to use or index' )
    parser.add_option( '--picard-ref', dest="picard_ref", help='The picard reference genome to use or index' )
    parser.add_option( '', '--lane1', dest="fastq_lane_1", help='The path to fastq files to use for the mapping' )
    parser.add_option( '', '--lane2', dest="fastq_lane_2", help='The path to fastq files to use for the mapping' )
    parser.add_option( '', '--lane3', dest="fastq_lane_3", help='The path to fastq files to use for the mapping' )
    parser.add_option( '', '--lane4', dest="fastq_lane_4", help='The path to fastq files to use for the mapping' )
    parser.add_option( '', '--lane5', dest="fastq_lane_5", help='The path to fastq files to use for the mapping' )
    parser.add_option( '', '--lane6', dest="fastq_lane_6", help='The path to fastq files to use for the mapping' )
    parser.add_option( '', '--lane7', dest="fastq_lane_7", help='The path to fastq files to use for the mapping' )
    parser.add_option( '', '--lane8', dest="fastq_lane_8", help='The path to fastq files to use for the mapping' )
    parser.add_option( '', '--gtf', dest="gtf_file", help="gtf file" )
    parser.add_option( '', '--refflat-file', dest="refFlat", help="refflat file" )
    parser.add_option( '', '--ribosomal-file', dest="ribosomal_file", help="ribosomal file" )
    parser.add_option( '--output-concatenated-fastq', dest="output_contatenate_fastq", help="the file to save")
    parser.add_option( '--output-tophat-bam', dest="output_tophat_bam", help="the file to save")
    parser.add_option( '--output-trinity-transcripts', dest="output_trinity_transcripts", help="the file to save")
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
    parser.add_option( '--output-log', dest="output_log", help='The file output' )
    parser.add_option( '--username', dest="username", help='Globus username for file transfer' )
    parser.add_option( '--goauth-token', dest="goauth_token", help='Globus token' )
    parser.add_option( '--source-ep', dest="source_ep", help='Globus source endpoint' )
    parser.add_option( '--destination-ep', dest="destination_ep", help='Globus destination endpoint' )
    parser.add_option( '--destination-dir', dest="destination_dir", help='Globus destination dir' )
    parser.add_option( '--deadline', dest="deadline", help='Globus transfer deadline' )
    (options, args) = parser.parse_args()

    final_outputs = [options.output_contatenate_fastq, options.output_tophat_bam, options.output_trinity_transcripts, options.output_picard_markDuplicates_html, options.output_picard_markDuplicates_directory, options.output_picard_alignmentSummary, options.output_picard_alignmentSummary_directory, options.output_picard_collectRnaSeqSummary, options.output_picard_collectRnaSeqSummary_directory, options.output_tophatStats, options.output_htseq_counts, options.output_htseq_log, options.output_fastqc_html, options.output_fastqc_text, options.output_log ]


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
