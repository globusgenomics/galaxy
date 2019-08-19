#!/usr/bin/python

#!/usr/bin/python

import optparse, os, shutil, sys, tempfile, glob, shlex, tarfile, pysam
from subprocess import *
import subprocess

CHUNK_SIZE = 2**20 #1mb

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
    parser.add_option( '', '--reference_name', dest='refName', action='store', type="string", help='' )
    parser.add_option( '', '--reference', dest='ref', action='store', type="string", help='' )
    parser.add_option( '', '--bam', dest='inputbam', action='store', type="string", help='' )
    #parser.add_option( '', '--summary', dest='summaryF', action='store', type="string", default=None, help='')
    parser.add_option( '', '--log', dest='logF', action='store', type="string", default=None, help='')
    parser.add_option( '', '--output', dest='outputF', action='store', type="string", default=None, help='')
#    parser.add_option( '', '--output-dir', dest='output_dir', action='store', type="string", default=None, help='If specified, the output directory for extra files.' )
    parser.add_option( '', '--checkChr', dest='checkchr', action='store', type="string", default=None, help='If specified, the output directory for extra files.' )
    (options, args) = parser.parse_args()

    # worker node tmp directory where all will be processed
    tmp_dir = tempfile.mkdtemp( prefix='tmp-tool-' )
    worker_input_dir = "%s/input" % tmp_dir
    _transfer_data(options, tmp_dir, worker_input_dir)
    #_cleanup_dir(worker_input_dir)

    # worker node tmp directory where all will be processed
#    if not os.path.exists(options.output_dir):
#        os.mkdir(options.output_dir)
        #subprocess.call(['chmod', '-R', '0777', options.output_dir])
        #os.chmod(options.output_dir, 0o777)

    # create input files from the input_dir
    tmp_bam = glob.glob("%s/*.bam" % worker_input_dir)
    input_bam = ''.join(tmp_bam)
    pysam.index(input_bam)

    if options.checkchr == "false":
        chrList="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"
    else:
        chrList="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"

#    cnvnator_cmd = "docker run -v %s:%s mustxyk/ubuntu-cnvnator /root/run_cnvnator.sh %s %s %s/cnvnator.root \"%s\" 500 %s/cnvnator.log %s/cnvnator.out %s" % (options.output_dir, options.output_dir, options.refName, options.ref, options.output_dir, chrList, options.output_dir, options.output_dir, input_bam)
    cnvnator_cmd = "docker run -v %s:%s -v /tmp:/tmp mustxyk/ubuntu-cnvnator /root/run_cnvnator.sh %s %s %s/cnvnator.root \"%s\" 500 %s/cnvnator.log %s/cnvnator.out %s" % (worker_input_dir, worker_input_dir, options.refName, options.ref, worker_input_dir, chrList, worker_input_dir, worker_input_dir, input_bam)
    print cnvnator_cmd

    stdout = tempfile.NamedTemporaryFile( prefix="cnvnator-stdout-", dir=worker_input_dir)
    stderr = tempfile.NamedTemporaryFile( prefix="cnvnator-stderr-", dir=worker_input_dir )
    proc = subprocess.Popen( args=cnvnator_cmd, stdout=stdout, stderr=stderr, shell=True, cwd=worker_input_dir )
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
    stdout.close()
    # copy files to final output locations
    #shutil.copy('%s/cnvnator.root' % options.output_dir, options.summaryF)
    shutil.copy('%s/cnvnator.log' % worker_input_dir, options.logF)
    shutil.copy('%s/cnvnator.out' % worker_input_dir, options.outputF)
    #cleanup_before_exit( options.output_dir )

if __name__=="__main__": __main__()
