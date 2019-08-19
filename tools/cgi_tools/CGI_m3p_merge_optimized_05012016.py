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

def create_wf_meta(options, output_dir, tmp_dir, inputDs):

    # reads

    lcm = 1
    #lcm = 4

    # get ncores
    ncpu = 24
    nthreads = ncpu/1

    # get memory
    jvm_heap = 31000000000

    # get the reference datasets
    lib_size = options.lib_size
    ## Setup workflow
    command_meta = {}

    #
    step = 0
    command_meta[step] = []
    input_files = [inputDs]
    outputD = options.output_dir
    outputF = "%s/merged_all_sample.txt" % output_dir
    output_files = [outputF]
    cmd = "merge_multiple_sample_results_mark_same_SV.pl %s %s" % (outputD, " ".join(inputDs))
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files": output_files})

    #
    step = 1
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    input_files = [inputF]
    outputD = "%s/final" % (output_dir)
    output_files = [outputD]
    cmd = "split_results_and_recall_regions_NEW.pl %s %s %s" % (outputD, inputF, lib_size)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputD, "output_files": output_files})

    return command_meta

def __main__():
    descr = "M3P workflow for tumor samples. A reference/normal sample is expected"
    parser = optparse.OptionParser(description=descr)
    parser.add_option( '', '--sample-dir', dest="inputSamples", action='append', help='The path to sample directory outputs from previous steps' )
    parser.add_option( '--lib-size', dest="lib_size", help="The size of the mate pair library, 5000 for example. usually make it a little bigger than the actual libarary size.")
    parser.add_option( '--output-dir', dest="output_dir", help="the dir to save")
    parser.add_option( '--output-log', dest="output_log", help='The file output' )
    parser.add_option( '--start-step', dest="start_step", default=0, help='The step to start at, default is at 0')
    (options, args) = parser.parse_args()

    final_outputs = [options.output_dir, options.output_log]

    #if len(args) > 0:
    #    parser.error('Wrong number of arguments')

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
    wf_meta = create_wf_meta(options, output_dir, tmp_dir, args)

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
