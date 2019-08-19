# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""
Runs the ITMI QC pipeline in optimized mode in one Node
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

def transfer_data(options, tmpdir, input_destination):
    python_script = "/opt/galaxy/tools/globus/transfer.py"
    transfer_log = "%s/transfer.log" % tmpdir
    cmd = "%s %s -k \"\" -c \"\" -o %s --goauth-token \"%s\" --source-ep=\"%s\" --source-path=\"%s\" --destination-ep=\"%s\" --destination-path=\"%s\" --deadline=\"1200\" --final %s --final_extra %s --type %s" % (python_script, options.username, options.output_transfer_log, options.goauth_token, options.source_ep, options.source_path, options.destination_ep, options.destination_path, transfer_log, input_destination, options.path_type)
    print cmd
    target_dir = tmpdir
    stdout = tempfile.NamedTemporaryFile( prefix="transfer-stdout-", dir=tmpdir )
    stderr = tempfile.NamedTemporaryFile( prefix="transfer-stderr-", dir=tmpdir )
    process = subprocess.Popen(cmd, shell=True, stderr=stderr, stdout=stdout, cwd=target_dir)
    rval = process.wait()
    if rval:
        stderr_transfer = sys.stderr
    else:
        stderr_transfer = sys.stdout
    stderr.flush()
    stderr.seek(0)
    while True:
        chunk = stderr.read( CHUNK_SIZE )
        if chunk:
            stderr_transfer.write( chunk )
        else:
            break
    stderr.close()

def transfer_multiple_data(options, datasets, tmpdir, input_destination):
    python_script = "/opt/galaxy/tools/globus/transfer_multiple.py"
    transfer_log = "%s/transfer.log" % tmpdir
    dataset_line = " "
    for dataset in datasets:
        dataset_line += " --dataset %s %s %s" % (dataset[0], dataset[1], dataset[2])
    cmd = "%s %s -k \"\" -c \"\" -o %s --goauth-token \"%s\" --source-ep=\"%s\" --destination-ep=\"%s\" %s --deadline=\"1200\"" % (python_script, options.username, options.output_transfer_log, options.goauth_token, options.source_ep, options.destination_ep, dataset_line)
    print cmd
    target_dir = tmpdir
    stdout = tempfile.NamedTemporaryFile( prefix="transfer-stdout-", dir=tmpdir )
    stderr = tempfile.NamedTemporaryFile( prefix="transfer-stderr-", dir=tmpdir )
    process = subprocess.Popen(cmd, shell=True, stderr=stderr, stdout=stdout, cwd=target_dir)
    rval = process.wait()
    if rval:
        stderr_transfer = sys.stderr
    else:
        stderr_transfer = sys.stdout
    stderr.flush()
    stderr.seek(0)
    while True:
        chunk = stderr.read( CHUNK_SIZE )
        if chunk:
            stderr_transfer.write( chunk )
        else:
            break
    stderr.close()


def _cleanup_dir(directory):
    for dir_file in sorted(os.listdir(directory)):
        infile = "%s/%s" % (directory, dir_file)
        print infile
        if os.path.isdir(infile):
            shutil.rmtree(infile)
            print "remove dir"
        elif not infile.endswith("fastq.gz") and os.path.isfile(infile):
            print "deleting file"
            os.remove(infile)


def _s3_root_path():
    return "s3://sulakhe/scratch/galaxy"

def _s3_put_cmd():
    return "/home/galaxy/s3cmd-1.6.1/s3cmd put"

