#!/usr/bin/python

import optparse, os, shutil, sys, tempfile, glob, shlex, vcf, pysam, tarfile
from subprocess import *
import subprocess, json

CHUNK_SIZE = 2**20 #1mb

def cleanup_before_exit( tmp_dir ):
    if tmp_dir and os.path.exists( tmp_dir ):
        shutil.rmtree( tmp_dir )

def str_to_dict(tmp_inputF):
    inputF = tmp_inputF.replace("'", "\"")
    inputF = inputF.replace(".dat","_files")
    data = json.loads(inputF)
    newdict={}
    for k,v in [(key,d[key]) for d in data for key in d]:
        if k not in newdict: newdict[k]=v
        else: newdict[k].append(v)
    return newdict

def __main__():
    parser = optparse.OptionParser()
    parser.add_option('-f','--input', dest="inputF", action='store', type="string", default=None, help='')
    parser.add_option( '', '--output', dest='outputF', action='store', type="string", default=None, help='')
    parser.add_option( '', '--out-dir', dest='output_dir', action='store', type="string", default=None, help='If specified, the output directory for extra files.' )

    (options, args) = parser.parse_args()
    new_dict = str_to_dict(options.inputF)
    #print(options.output_dir)
    if not os.path.exists(options.output_dir): 
        os.mkdir(options.output_dir)
    input_dir = "%s/input" % options.output_dir
    if not os.path.exists(input_dir):
        os.mkdir(input_dir)

    # copy all input data to a input directory
    for key, value in new_dict.iteritems():
        for i in range(0,len(value)):
            path=input_dir+'/'+os.path.basename(value[i])
            shutil.copytree(value[i],path)
    #print "new_dict: %s\n" % new_dict
    new_dict1 = json.dumps(new_dict)
    #print "new_dict1: %s\n" % new_dict1
    # test = json.dumps(new_dict)
    # create a Rscript command with dictionary, input folder,  and output
    sleuth_cmd = "Rscript /opt/galaxy/tools/sleuth/sleuth.R -o %s -f '%s' -d %s" % (options.outputF, new_dict1, input_dir) 
    print "sleuth cmd: %s " % sleuth_cmd

    stdout = tempfile.NamedTemporaryFile( prefix="sleuth-stdout-", dir=options.output_dir)
    stderr = tempfile.NamedTemporaryFile( prefix="sleuth-stderr-", dir=options.output_dir)
    proc = subprocess.Popen( args=sleuth_cmd, stdout=stdout, stderr=stderr, shell=True, cwd=options.output_dir )
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
#    shutil.copy('%s/DU_K562_HINT.bed' % output_dir, options.outputF)
#    cleanup_before_exit( options.output_dir )

if __name__=="__main__": __main__()

