# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""
Runs the CGI m3p panel workflow for SV detection 
See below for options
"""

import time, optparse, os, shutil, subprocess, sys, tempfile, multiprocessing, gzip
CHUNK_SIZE = 2**20 #1mb

#def _get_sampleName(forward, reverse):
def _get_sampleName(bam):
    #common = []
    #count = 0
    #f_vals = list(forward)
    #r_vals = list(reverse)
    #for i in f_vals:
    #    if f_vals[count] == r_vals[count]:
    #        common.append(f_vals[count])
    #    else:
    #        return "".join(common)
    #    count += 1
    values = bam.split(".")
    return os.path.basename(values[0])

def time_stamp():
    return time.strftime('%d/%m/%Y %H:%M:%S\n', time.localtime(time.time()))

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

    #lcm = 1
    lcm = 4

    # get ncores
    ncpu = 32
    nthreads = ncpu/lcm

    # get memory
    jvm_heap = 256000000000/lcm

    # get the reference datasets
    #source_path = "/scratch/galaxy/test/CGI_TEST/called"
    source_path = options.source_path
    #bwa_reference_db = options.bwa_ref
    #fastq = options.fastq
    #rfastq = options.rfastq
    input_bam = options.input_bam
    upper = options.upper
    lib_size = options.lib_size
    read_length = options.read_length
    min_pairs = options.min_pairs
    rev_comp = options.rev_comp
    if options.input_normal_log:
        inputNormalLog = options.input_normal_log
    else:
        inputNormalLog = "None"

    #sample_name = _get_sampleName(os.path.basename(fastq), os.path.basename(rfastq))
    sample_name = _get_sampleName(input_bam)
    ## Setup workflow
    command_meta = {}

    #
    #step = 0
    #command_meta[step] = []
    #input_files = [fastq, rfastq]
    #outputF = "%s/%s-bwa-sorted.bam" % (output_dir, sample_name)
    #outputPrefix = "%s/%s-bwa-sorted" % (output_dir, sample_name)
    #output_files = [outputF]
    #cmd = "bwa mem %s %s %s | samtools view -Sb - | samtools sort -n -m %s - %s" % (bwa_reference_db, fastq, rfastq, str(jvm_heap), outputPrefix)
    #command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files": output_files})

    #
    step = 0
    command_meta[step] = []
    input1 = input_bam
    #input1 = command_meta[step-1][0]["output_file"]
    input_files = [input1]
    outputF = "%s/%s-bwa-sorted.sam" % (output_dir, sample_name)
    output_files = [outputF]
    cmd = "samtools view -X %s -o %s" % (input1, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files": output_files})

    #
    step = 1
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputIDF = "%s/%s-IDF-abnormal.txt" % (output_dir, sample_name)
    outputLog = "%s/%s-IDF-log.txt" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputIDF, outputLog]
    cmd = "classify_ID_sorted_SAM_BWA-MEM.pl %s %s %s %s %s %s %s" % (sample_name, inputF, outputIDF, outputLog, upper, lib_size, read_length)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputIDF, "output_files":output_files})

    #
    step = 2
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    #outputD = "%s/%s-step5_dir" % (output_dir, sample_name)
    outputD = options.output_dir
    input_files = [inputF]
    output_files = [outputD]
    cmd = "get_nr_sorted_reads_from_both_end_mapped.pl %s %s %s" % (sample_name, inputF, outputD)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputD, "output_files":output_files})

    #
    step = 3
    command_meta[step] = []
    inputD = "%s/%s" % (command_meta[step-1][0]["output_file"], sample_name)
    inputControlD = options.input_normal
    input_files = [inputD, inputControlD]
    outputD = "%s/merged_with_control" % (output_dir)
    output_files = [outputD]
    cmd = "merge_with_control.pl %s %s %s" %  (inputControlD, inputD, outputD)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputD, "output_files": output_files})

    #
    step = 4
    command_meta[step] = []
    inputD = command_meta[step-1][0]["output_file"]
    inputTumorLog = command_meta[step-3][0]["output_files"][1]
    input_files = [inputD]
    #ratio = 1  # ratio of normal SAM reads vs actual sample SAM reads
    outputD = "%s/results" % (output_dir)
    output_files = [outputD]
    cmd = "nominate_fusion_regions.pl %s %s %s %s %s %s %s %s" %  (inputD, inputNormalLog, inputTumorLog, outputD, lib_size, min_pairs, sample_name, rev_comp)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputD, "output_files": output_files})

    #
    step = 5
    command_meta[step] = []
    inputD = output_dir
    outputD = options.output_dir
    input_files = [inputD]
    #ratio = 1  # ratio of normal SAM reads vs actual sample SAM reads
    output_files = [outputD]
    cmd = "cp -r %s %s" %  (inputD, outputD)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputD, "output_files": output_files})

    return command_meta

def __main__():
    descr = "M3P workflow for tumor samples. A reference/normal sample is expected"
    parser = optparse.OptionParser(description=descr)
    parser.add_option( '--sample-name', dest="sample_name", help='sample name' )
    parser.add_option( '--input-bam', dest="input_bam", help='The BAM file' )
    parser.add_option( '--bwa-ref', dest="bwa_ref", help='The BWA reference genome path to use or index' )
    parser.add_option( '', '--fastq', dest="fastq", help='The path to forward fastq files to use for the mapping' )
    parser.add_option( '', '--rfastq', dest="rfastq", help='The path to reversefastq files to use for the mapping' )
    parser.add_option( '--upper', dest="upper", help="the distance of two ends on the same chromosome to be considered as SV")
    parser.add_option( '--lib-size', dest="lib_size", help="The size of the mate pair library, 5000 for example. usually make it a little bigger than the actual libarary size.")
    parser.add_option( '--read-length', dest="read_length", help="read length")
    parser.add_option( '--min-pairs', dest="min_pairs", help="minimum number of pairs to support a SV")
    parser.add_option( '--rev-comp', dest="rev_comp", help="if the original fastq files were reversed and complimented")
    parser.add_option( '--input-normal', dest="input_normal", help="input normal directory")
    parser.add_option( '--input-normal-log', dest="input_normal_log", default="None", help="Input normal log file")
    parser.add_option( '--output-dir', dest="output_dir", help="the dir to save")
    parser.add_option( '--output-log', dest="output_log", help='The file output' )
    parser.add_option( '--start-step', dest="start_step", default=0, help='The step to start at, default is at 0')
    parser.add_option( '--source', dest="source_path", help="source of supplementary scripts")
    (options, args) = parser.parse_args()

    final_outputs = [options.output_dir, options.output_log]

    if len(args) > 0:
        parser.error('Wrong number of arguments')

    if options.start_step:
        start_step = int(options.start_step)
    else:
        start_step = 0

    # make temp directory
    #if options.output_dir:
    #    output_dir = options.output_dir
    #    if not os.path.exists(output_dir):
    #        os.makedirs(output_dir)
    #else:
    #    output_dir = tempfile.mkdtemp(prefix="optimized-")
    #tmp_dir = tempfile.mkdtemp(prefix="%s/optimized-tmp-" % output_dir)
    output_dir = tempfile.mkdtemp(prefix="optimized-")
    tmp_dir = tempfile.mkdtemp(prefix="%s/optimized-tmp-" % output_dir)

    # create the meta dictionary of the commands to run
    wf_meta = create_wf_meta(options, output_dir, tmp_dir)

    #print wf_meta
    os.chdir(output_dir)
    for step_id, step_list in wf_meta.iteritems():
        if start_step <= step_id:
            print "START STEP %s: %s" % (step_id, time_stamp())
            step_input_files = []
            if len(step_list) == 1:
                for job in step_list:
                    print "run step %s:\n%s" % (step_id, job['cl'])
                    run_cmd ( job['cl'], tmp_dir, "running some job")
                    for ifile in job['input_files']:
                        if ifile not in step_input_files:
                            step_input_files.append(ifile)
            print "END STEP %s: %s" % (step_id, time_stamp())

    #clean up
    #shutil.rmtree(output_dir)

if __name__ == "__main__":
    __main__()
