#!/usr/bin/env python

"""
A wrapper script for running the PEMAPPER on the cloud
Including file transfer which lives on the user's endpoint 
"""

import tarfile, time, sys, optparse, os, tempfile, glob, subprocess, shutil
from collections import namedtuple

CHUNK_SIZE = 2**20 #1mb

_ntuple_diskusage = namedtuple('usage', 'total used free')

def disk_usage(path):
    """Return disk usage statistics about the given path.

    Returned valus is a named tuple with attributes 'total', 'used' and
    'free', which are the amount of total, used and free space, in bytes.
    """
    st = os.statvfs(path)
    free = st.f_bavail * st.f_frsize
    total = st.f_blocks * st.f_frsize
    used = (st.f_blocks - st.f_bfree) * st.f_frsize
    return _ntuple_diskusage(total, used, free)

def cleanup_before_exit( tmp_dir ):
    if tmp_dir and os.path.exists( tmp_dir ):
        shutil.rmtree( tmp_dir )

def _transfer_data(options, tmpdir, input_destination):
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

def _create_input_configs(input_dir):
    print input_dir
    forward = []
    reverse = []
    count = 0
    for infile in sorted(os.listdir(input_dir)):
        # assume files are fastq
        if infile.endswith(".fastq.gz"):
            if "_R1_" in infile or "_1_" in infile or "_1.fastq.gz" in infile:
                forward.append(infile)
                #link_name = "input%s.fastq.gz" % (count)
                #forward.append("input%s.fastq.gz" % (count))
            elif "_R2_" in infile or "_2_" in infile or "_2.fastq.gz" in infile:
                reverse.append(infile)
                #link_name = "input%s.fastq.gz" % (count)
                #reverse.append("input%s.fastq.gz" % (count))
            #os.symlink(infile, link_name)
            count += 1
    f_config = "file1.txt"
    r_config = "file2.txt" 

    return_files = []
    if len(forward) > 0:
        fh1 = open("%s/%s" % (input_dir, f_config), "w")
        for i in forward:
            fh1.write("%s\n" % i)
        fh1.close()
        return_files.append(f_config)
    if len(reverse) > 0:
        fh1 = open("%s/%s" % (input_dir, r_config), "w")
        for i in reverse:
            fh1.write("%s\n" % i)
        fh1.close()
        return_files.append(r_config)

    return return_files

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

def __main__():
    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option("-u", "--username", dest="username", help="username")
    parser.add_option("-c", "--cert", dest="cert_file", help="client cert file", metavar="CERT_FILE")
    parser.add_option("-k", "--key", dest="key_file", help="client key file", metavar="KEY_FILE")
    parser.add_option("-b", "--base-url", dest="base_url", help="alternate base URL", metavar="URL")
    parser.add_option("-o", "--out-transfer-log", dest="output_transfer_log", help="write log output to PATH", metavar="PATH")
    parser.add_option("--source-ep", dest="source_ep", help="Endpoint to transfer from")
    parser.add_option("--source-path", dest="source_path", help="Source endpoint filepath to transfer")
    parser.add_option("--extra-source-path", dest="extra_source_path", help="Source endpoint filepath to transfer for BAM Index file")
    parser.add_option("--destination-ep", dest="destination_ep", help="Endpoint to transfer to")
    parser.add_option("--destination-path", dest="destination_path", help="Destination endpoint filepath to transfer")
    parser.add_option("-d", "--deadline", dest="deadline", help="Deadline for transfer in minutes.")
    parser.add_option("-g", "--galaxy-dataset-id", dest="galaxy_dataset_id", help="Galaxy Dataset Id For This Transfer")
    parser.add_option('-a', '--goauth-token', dest="goauth_token", help="Use the Globus Access Token as the authentication method for this transfer")
    parser.add_option("--type", dest="path_type", help="directory or file")
    parser.add_option("--pemapper-ref", dest="pemapper_ref", help="pemapper ref") 
    parser.add_option( '-p', '--pass_through', dest='pass_through_options', action='append', type="string", help='These options are passed through directly to GATK, without any modification.' )
    parser.add_option("--out-summary", dest="output_summary", help="summary output")
    parser.add_option("--out-pileup", dest="output_pileup", help="pileup output")
    parser.add_option("--out-indel", dest="output_indel", help="indel output")
    parser.add_option("--out-mfiles", dest="output_mfiles", help="mfiles output")
    parser.add_option("--out-config1", dest="output_config_forward", help="config1 output")
    parser.add_option("--out-config2", dest="output_config_reverse", help="config2 output")
    (options, args) = parser.parse_args()
  
    # worker node tmp directory where all will be processed 
    tmp_dir = tempfile.mkdtemp( prefix='tmp-tool-' )
    #tmp_dir = tempfile.mkdtemp( prefix='/scratch/galaxy/tmp/tmp_dirtmp-tool-' )
    worker_input_dir = "%s/input" % tmp_dir
    #os.mkdir(worker_input_dir)
    _transfer_data(options, tmp_dir, worker_input_dir)
    _cleanup_dir(worker_input_dir)

    # create input files from the input_dir
    os.chdir(worker_input_dir)
    input_configs = _create_input_configs(worker_input_dir)
    paired_single = "sa"
    shutil.copy("%s/%s" % (worker_input_dir, input_configs[0]), options.output_config_forward)
    if len(input_configs) > 1:
        paired_single = "pa"
        shutil.copy("%s/%s" % (worker_input_dir,input_configs[1]), options.output_config_reverse)
        

    pass_through = ""
    if options.pass_through_options:
        pass_through = ' '.join( options.pass_through_options )

    pemapper_cmd = "pemapper sampleName %s %s %s %s" % (options.pemapper_ref, paired_single, " ".join(input_configs), pass_through)
    print pemapper_cmd

    stdout = tempfile.NamedTemporaryFile( prefix="pemapper-stdout-", dir=tmp_dir )
    stderr = tempfile.NamedTemporaryFile( prefix="pemapper-stderr-", dir=tmp_dir ) 
    proc = subprocess.Popen( args=pemapper_cmd, stdout=stdout, stderr=stderr, shell=True, cwd=worker_input_dir )
    return_code = proc.wait()
    
    if return_code:
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

    # copy files to final output locations
    shutil.copy("%s/sampleName.summary.txt" % worker_input_dir, options.output_summary)
    shutil.copy("%s/sampleName.pileup.gz" % worker_input_dir, options.output_pileup)
    shutil.copy("%s/sampleName.indel.txt.gz" % worker_input_dir, options.output_indel)
    shutil.copy("%s/transfer.log" % tmp_dir, options.output_transfer_log)
    fh = open(options.output_transfer_log, "a")
    fh.write("\n%s\n" % pemapper_cmd)
    fh.close()

    # tar and zip the mfiles
    mfiles = glob.glob("%s/*.mfile" % worker_input_dir)
    tarfile_name = "%s/mfiles.tar.gz" % tmp_dir
    tar = tarfile.open(tarfile_name, "w:gz")
    for name in mfiles:
        tar.add(name)
    tar.close()
    shutil.copy(tarfile_name, options.output_mfiles)

    #cleanup_before_exit( tmp_dir )

if __name__=="__main__": __main__()
