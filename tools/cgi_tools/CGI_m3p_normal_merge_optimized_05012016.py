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

def create_wf_meta(options, output_dir, tmp_dir):

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
    ## Setup workflow
    command_meta = {}

    #
    step = 0
    command_meta[step] = []
    inputDs = options.inputSamples
    input_files = [inputDs]
    outputD = options.output_dir
    outputF = "%s/merged_all_sample.txt" % output_dir
    output_files = [outputF]
    cmd = "perl %s/merge_normal_sample.pl %s/merged %s" % (source_path, outputD, " ".join(inputDs))
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files": output_files})

    return command_meta

def __main__():
    descr = "M3P workflow for tumor samples. A reference/normal sample is expected"
    parser = optparse.OptionParser(description=descr)
    parser.add_option( '--sample-name', dest="sample_name", help='sample name' )
    parser.add_option( '', '--sample-dir', dest="inputSamples", action='append', help='The path to sample directory outputs from previous steps' )
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
    if options.output_dir:
        output_dir = options.output_dir
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    else:
        output_dir = tempfile.mkdtemp(prefix="/scratch/galaxy/test/CGI_TEST/tmp/optimized-")
        #output_dir = tempfile.mkdtemp(prefix="optimized-")
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
