# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""
Runs the FermiKit Assembly program 
See below for options
"""

import time, optparse, os, shutil, subprocess, sys, tempfile, multiprocessing, gzip
CHUNK_SIZE = 2**20 #1mb

def time_stamp():
    return time.strftime('%d/%m/%Y %H:%M:%S\n', time.localtime(time.time()))

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

def run_cmd ( cmd , wd_tmpdir, descriptor):
    #print cmd
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

def __main__():
    descr = ""
    parser = optparse.OptionParser(description=descr)
    parser.add_option( '-s', '--genome-size', dest="genomeSize", help='Approximate genome size' )
    parser.add_option( '-f', '--fastq', dest="fastq", help='The (forward) fastq file to use for the mapping' )
    parser.add_option( '-F', '--rfastq', dest="rfastq", help='The reverse fastq file to use for mapping if paired-end data' )
    parser.add_option( '-o', '--output-assembly', dest="output_assembly", help="The file to save the output" )
    parser.add_option( '-l', '--read-length', dest="readLength", help='The primary read length' )
    parser.add_option( '-t', '--trim-adapters', dest="trimAdapters", action="store_true", help='trim adapters from reads' )
    (options, args) = parser.parse_args()

    if len(args) > 0:
        parser.error('Wrong number of arguments')

    # make temp directory 
    output_dir = tempfile.mkdtemp(prefix="optimized-")
    tmp_dir = tempfile.mkdtemp(prefix="optimized-tmp-")

    # create the meta dictionary of the commands to run
    wf_meta = {}
    output_prefix = "%s/prefix" % output_dir
    make_file = "%s.mak" % output_prefix
    cmd = "fermi2.pl unitig -s%s -l%s -t%s -p %s " % (options.genomeSize, options.readLength, get_ncores(), output_prefix)

    input_cmd = None
    if options.rfastq is not None:
        input_cmd = "seqtk mergepe %s %s" % (options.fastq, options.rfastq)
        if options.trimAdapters is not None:
            input_cmd += " | trimadapt-mt -p4"
    else:
        if options.trimAdapters is not None:
            input_cmd = "trimadapt-mt %s" %  options.fastq
        else:
            input_cmd = options.fastq


    cmd += " \"%s\" > %s" % (input_cmd, make_file)
    wf_meta['1'] = [{'cl' : cmd}]

    cmd = "make -f %s" % make_file
    wf_meta['2'] = [{'cl' : cmd}]

    output_gz = "%s/prefix.mag.gz" % (output_dir)
    cmd = "mv %s %s" % (output_gz, options.output_assembly)
    wf_meta['3'] = [{'cl' : cmd}]
  
    #print wf_meta
    for step_id, step_list in sorted(wf_meta.iteritems()):
        print "START STEP %s: %s" % (step_id, time_stamp())
        step_input_files = []
        for job in step_list:
            print "run step %s:\n%s" % (step_id, job['cl'])
            run_cmd ( job['cl'], tmp_dir, "running some job")
        print "END STEP %s: %s" % (step_id, time_stamp())

    #clean up
    shutil.rmtree(output_dir)

if __name__ == "__main__":
    __main__()
