# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""
Runs the CGI m3p panel workflow for SV detection 
See below for options
"""

import time, optparse, os, shutil, subprocess, sys, tempfile, gzip
CHUNK_SIZE = 2**20 #1mb

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

def _input_pairs(input_dir):
    print input_dir
    forward = []
    reverse = []
    count = 0
    for infile in sorted(os.listdir(input_dir)):
        # assume files are fastq
        if infile.endswith(".fastq.gz"):
            if "_R1_" in infile or "_1_" in infile or "_1.fastq.gz" in infile:
                forward.append(infile)
            elif "_R2_" in infile or "_2_" in infile or "_2.fastq.gz" in infile:
                reverse.append(infile)
    return (forward, reverse)

def _get_sampleName(forward, reverse):
    common = []
    count = 0
    f_vals = list(forward)
    r_vals = list(reverse)
    for i in f_vals:
        if f_vals[count] == r_vals[count]:
            common.append(f_vals[count])
        else:
            return "".join(common)        
        count += 1

def _time_stamp():
    return time.strftime('%d/%m/%Y %H:%M:%S\n', time.localtime(time.time()))

def __main__():
    descr = "M3P workflow for normal samples. A reference/normal sample is expected"
    parser = optparse.OptionParser(description=descr)
    parser.add_option( '--project', dest="project_name", help='project name' )
    parser.add_option( '--bwa-ref', dest="bwa_ref", help='The BWA reference genome path to use or index' )
    parser.add_option( '--input-dir', dest="input_dir", help="input directory with fastq samples")
    parser.add_option( '--upper', dest="upper", help="the distance of two ends on the same chromosome to be considered as SV")
    parser.add_option( '--lib-size', dest="lib_size", help="The size of the mate pair library, 5000 for example. usually make it a little bigger than the actual libarary size.")
    parser.add_option( '--read-length', dest="read_length", help="read length")
    parser.add_option( '--min-pairs', dest="min_pairs", help="minimum number of pairs to support a SV")
    parser.add_option( '--rev-comp', dest="rev_comp", help="if the original fastq files were reversed and complimented")
    parser.add_option( '--output-dir', dest="output_dir", help="the dir to save")
    parser.add_option( '--output-log', dest="output_log", help='The file output' )
    parser.add_option( '--script-source', dest="script_source", help="Location of scripts" )
    (options, args) = parser.parse_args()

    if len(args) > 0:
        parser.error('Wrong number of arguments')

    # get the sequence pairs
    forward, reverse = _input_pairs(options.input_dir)

    # make temp directory
    if options.output_dir:
        output_dir = options.output_dir
        project_dir = "%s/%s" % (output_dir, options.project_name)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        if not os.path.exists(project_dir):
            os.makedirs(project_dir)
    tmp_dir = tempfile.mkdtemp(prefix="%s/optimized-tmp-" % output_dir)

    count = 0
    sample_dirs = []
    for seq in forward:
        fastq = forward[count]
        rfastq = reverse[count]
        sample_name = _get_sampleName(fastq, rfastq)
        output_dir = "%s/%s" % (project_dir, sample_name)
        output_log = "%s/%s.log.txt" % (output_dir, sample_name)
        output_idf_log = "%s/%s.idf-log.txt" % (output_dir, sample_name)
        sample_dirs.append(output_dir)
        # submit the sample
        cmd = "python %s/CGI_m3p_normal_optimized_05012016.py --source %s/called --sample-name %s --bwa-ref %s --fastq %s --rfastq %s --upper %s --lib-size %s --read-length %s --output-dir %s --output-log %s --output-idf-log %s" % (options.script_source, options.script_source, sample_name, options.bwa_ref, "%s/%s" % (options.input_dir, fastq), "%s/%s" % (options.input_dir, rfastq), options.upper, options.lib_size, options.read_length, output_dir, output_log, output_idf_log)
        print "START STEP %s: %s" % (sample_name, _time_stamp())
        run_cmd ( cmd, tmp_dir, "running sample %s" % sample_name)
        print "END STEP %s: %s" % (sample_name, _time_stamp())
        count += 1

    # merge step
    samples = ""
    for sample_dir in sample_dirs:
        samples += " --sample-dir %s " % sample_dir
    output_log = "%s/%s.log.txt" % (output_dir, options.project_name)
    cmd = "python %s/CGI_m3p_normal_merge_optimized_05012016.py --source %s/called --sample-name %s --output-dir %s --output-log %s %s" % (options.script_source, options.script_source, options.project_name, project_dir, output_log, samples)
    print "START STEP %s: %s" % ("Merge", _time_stamp())
    run_cmd ( cmd, tmp_dir, "running sample %s" % sample_name)
    print "END STEP %s: %s" % ("Merge", _time_stamp())

    #clean up
    #shutil.rmtree(output_dir)

if __name__ == "__main__":
    __main__()
