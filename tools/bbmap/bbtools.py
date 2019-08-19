# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""
Wrappers for the BBmap set tools
"""

import time, optparse, os, shutil, subprocess, sys, tempfile, multiprocessing, gzip
CHUNK_SIZE = 2**20 #1mb

def time_stamp():
    return time.strftime('%d/%m/%Y %H:%M:%S\n', time.localtime(time.time()))

def get_seq_link(fastq_file, direction, directory):
    ftype = file_type(fastq_file)
    link_path = None
    if ftype == "gz":
        link_path = "%s/in_%s.fastq.gz" % (directory, direction)
    else:
        link_path = "%s/in_%s.fastq" % (directory, direction)

    os.symlink(fastq_file, link_path)
    return link_path

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

def __main__():
    descr = "snpir_optimized_workflow.py: Runs the Meissnert SNPiR RNA-variant discovery pipeline in optimized mode in one Node"
    parser = optparse.OptionParser(description=descr)
    parser.add_option( '-p', '--pass_through', dest='pass_through_options', action='append', type="string", help='These options are passed through directly to GATK, without any modification.' )
    parser.add_option( '--command', dest="command", help='The bbtool command to use, bbmap, bbmerge, bbduk' )
    parser.add_option( '-f', '--in', dest="in1", help='The (forward) fastq file to use for the mapping' )
    parser.add_option( '-F', '--in2', dest="in2", help='The reverse fastq file to use for mapping if paired-end data' )
    parser.add_option( '', '--log', dest="output_log", help="The log file" )
    (options, args) = parser.parse_args()

    if len(args) > 0:
        parser.error('Wrong number of arguments')

    # make temp directory 
    tmp_dir = tempfile.mkdtemp(prefix="bbtool-tmp-")

    # create links for the inputs
    # check if input is in gz format or not
    in1_link = get_seq_link(options.in1, "1", tmp_dir) 
    if options.in2:
        in2_link = get_seq_link(options.in2, "2", tmp_dir)

    # prepare the command
    cmd = "%s in1=%s " % (options.command, in1_link)
    if options.in2:
        cmd += "in2=%s " % in2_link

    if options.pass_through_options:
        cmd += ' '.join( options.pass_through_options )
    else:
        cmd += ''

    cmd += "> %s" % options.output_log 

    print "run: %s" % (cmd)
    run_cmd ( cmd, tmp_dir, "running some job")

    #clean up
    shutil.rmtree(tmp_dir)

if __name__ == "__main__":
    __main__()
