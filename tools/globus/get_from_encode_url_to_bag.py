#!/usr/bin/env python

"""
"""

import optparse, os, subprocess, sys, tempfile, shutil, requests
CHUNK_SIZE = 2**20 #1mb

def run_cmd ( cmd, descriptor):
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

def stop_err( msg ):
    sys.stderr.write( '%s\n' % msg )
    sys.exit()

def __main__():
    parser = optparse.OptionParser()
    parser.add_option( '-u', '--url', dest='source',  help='where to get data' )
    parser.add_option( '-o', '--output', dest="output", help="output file" )
    parser.add_option( '-d', '--output-dir', dest="output_dir", help="output directory" )
    (options, args) = parser.parse_args()
    #print args

    if not os.path.exists(options.output_dir):
        os.makedirs(options.output_dir)
    inurls = options.source.split(",")
    fw = open(options.output, "w")

    index = None
    if len(inurls) > 1:
        index = 1
    else:
        index = ""
    standard_name = "sample"

    for i in inurls:
        i = i.replace("__at__", "@")
        orig_filename = i.split("/")[-1]
        if len(inurls) > 1:
            filename = "standard_name_%s.fastq.gz" % str(index)
            index += 1
        else:
            filename = "standard_name.fastq.gz"
        localname = os.path.join(options.output_dir, filename)
        cmd = "wget -P %s %s" % (options.output_dir, i) 
        run_cmd(cmd, "running wget: %s" % i)
        os.rename("%s/%s" % (options.output_dir, orig_filename), "%s/%s" % (options.output_dir, filename))
        fw.write("%s\n" % i)
    fw.close()

if __name__ == "__main__" : __main__()
