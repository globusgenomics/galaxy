# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""
Runs the CGI m3p panel workflow for SV detection 
See below for options
"""

import time, optparse, os, shutil, subprocess, sys, tempfile, multiprocessing, gzip
CHUNK_SIZE = 2**20 #1mb

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

def create_wf_meta(options, output_dir):

    # reads
    sample_name = options.sample_name

    lcm = 1
    #lcm = 4

    # get ncores
    ncpu = 24
    nthreads = ncpu/1

    # get memory
    jvm_heap = 31000000000

    # get the reference datasets
    #source_path = "/scratch/galaxy/test/CGI_TEST/called"
    source_path = options.source_path
    bwa_reference_db = options.bwa_ref
    fastq = options.fastq
    rfastq = options.rfastq
    upper = options.upper
    lib_size = options.lib_size
    read_length = options.read_length

    ## Setup workflow
    command_meta = {}

    #
    step = 0
    command_meta[step] = []
    input_files = [fastq, rfastq]
    outputF = "%s/%s-bwa.sam" % (output_dir, sample_name)
    output_files = [outputF]
    cmd = "bwa mem %s %s %s > %s" % (bwa_reference_db, fastq, rfastq, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files": output_files})

    #
    step = 1
    command_meta[step] = []
    input1 = command_meta[step-1][0]["output_file"]
    input_files = [input1]
    outputF = "%s/%s-bwa.bam" % (output_dir, sample_name)
    output_files = [outputF]
    cmd = "samtools view -Sb %s > %s; rm -rf %s" % (input1, outputF, input1)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files": output_files})

    #
    step = 2
    command_meta[step] = []
    input1 = command_meta[step-1][0]["output_file"]
    input_files = [input1]
    outputF = "%s/%s-bwa-sorted.bam" % (output_dir, sample_name)
    outputPrefix = "%s/%s-bwa-sorted" % (output_dir, sample_name)
    output_files = [outputF]
    cmd = "samtools sort -n -m %s %s %s; rm -rf %s" % (str(jvm_heap), input1, outputPrefix, input1)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files": output_files})

    #
    step = 3
    command_meta[step] = []
    input1 = command_meta[step-1][0]["output_file"]
    input_files = [input1]
    outputF = "%s/%s-bwa-sorted.sam" % (output_dir, sample_name)
    output_files = [outputF]
    cmd = "samtools view -X %s -o %s; rm -rf %s" % (input1, outputF, input1)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files": output_files})

    #
    step = 4
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputIDF = "%s/%s-IDF-abnormal.txt" % (output_dir, sample_name)
    outputLog = "%s/%s-IDF-log.txt" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputIDF, outputLog]
    cmd = "perl %s/classify_ID_sorted_SAM_BWA-MEM.pl %s %s %s %s %s %s %s" % (source_path, sample_name, inputF, outputIDF, outputLog, upper, lib_size, read_length)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputIDF, "output_files":output_files})

    #
    step = 5
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    #outputD = "%s/%s-step5_dir" % (output_dir, sample_name)
    outputD = options.output_dir
    input_files = [inputF]
    output_files = [outputD]
    cmd = "perl %s/get_nr_sorted_reads_from_both_end_mapped.pl %s %s %s" % (source_path, sample_name, inputF, outputD)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputD, "output_files":output_files})

    return command_meta

def __main__():
    descr = "M3P workflow for normal samples. Generates output in output_dir"
    parser = optparse.OptionParser(description=descr)
    parser.add_option( '--sample-name', dest="sample_name", help='sample name' )
    parser.add_option( '--bwa-ref', dest="bwa_ref", help='The BWA reference genome path to use or index' )
    parser.add_option( '', '--fastq', dest="fastq", help='The path to forward fastq files to use for the mapping' )
    parser.add_option( '', '--rfastq', dest="rfastq", help='The path to reversefastq files to use for the mapping' )
    parser.add_option( '--upper', dest="upper", help="the distance of two ends on the same chromosome to be considered as SV")
    parser.add_option( '--lib-size', dest="lib_size", help="library size")
    parser.add_option( '--read-length', dest="read_length", help="read length")
    parser.add_option( '--output-dir', dest="output_dir", help='The dir to save' )
    parser.add_option( '--output-idf-log', dest="output_idf_log", help="log file for IDF" )
    parser.add_option( '--output-log', dest="output_log", help='The file output' )
    parser.add_option( '--start-step', dest="start_step", help='The step to start at, default is at 0')
    parser.add_option( '--source', dest="source_path", help="source of supplementary scripts")
    (options, args) = parser.parse_args()

    final_outputs = [options.output_dir, options.output_log, options.output_idf_log]

    if len(args) > 0:
        parser.error('Wrong number of arguments')

    if options.start_step:
        start_step = int(options.start_step)
    else:
        start_step = 0

    # make temp directory 
    if options.output_dir:
        output_dir = options.output_dir
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    else:
        output_dir = tempfile.mkdtemp(prefix="/scratch/galaxy/test/CGI_TEST/tmp/optimized-")
        #output_dir = tempfile.mkdtemp(prefix="optimized-")
    tmp_dir = tempfile.mkdtemp(prefix="%s/optimized-tmp-" % output_dir)

    # create the meta dictionary of the commands to run
    wf_meta = create_wf_meta(options, output_dir)

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
