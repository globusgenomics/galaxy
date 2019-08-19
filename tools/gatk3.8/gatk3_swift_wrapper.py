#!/usr/bin/env python
#Dan Blankenberg

"""
A wrapper script for running the GenomeAnalysisTK.jar using Swift.
"""

import sys, optparse, os, tempfile, subprocess, shutil, glob
from binascii import unhexlify
from string import Template
import time  #liubo added
CHUNK_SIZE = 2**20 #1mb

## GLOBALS    
swift_bin = '/mnt/galaxyTools/tools/swift/swift-0.94.1/bin/swift'
sites_file = '/opt/galaxy/tools/gatk3/gatk3_hc_sites.xml'
tc_file = '/opt/galaxy/tools/swift/tc.data'
swift_file = '/opt/galaxy/tools/gatk3/haplotypeCaller.swift'

def cleanup_before_exit( tmp_dir ):
    if tmp_dir and os.path.exists( tmp_dir ):
        shutil.rmtree( tmp_dir )

def __main__():
    print "Start time:"
    print time.strftime('%d/%m/%Y %H:%M:%S\n', time.localtime(time.time()))

    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '', '--input_dir', dest='input_dir', action='store', type="string", help='Location of input BAM files' )
    parser.add_option( '', '--output_log', dest='output_log', action='store', type="string", help='Location of output log files' )
    parser.add_option( '', '--output_dir', dest='output_dir', action='store', type="string", help='Location of output VCF files' )
    parser.add_option( '', '--gatk_options', dest='gatk_options', action='append', type="string" )
    (options, args) = parser.parse_args()

    tmp_dir = tempfile.mkdtemp(dir=os.getcwd())   
    completed_samples = 0
    input_samples = {}
    input_samples_count = 0
    for i in os.listdir(options.input_dir):
        if i.endswith("bam"):
            sample_name = os.path.basename(i)
            bam_file = "%s/%s" % (options.input_dir, i)
            index_name = "%s/%s.bai" % (options.input_dir, i)
            vcf_name = "%s/%s.vcf" % (options.output_dir, sample_name)
            vcf_idx_name = "%s/%s.vcf.idx" % (options.output_dir, sample_name)
            input_samples_count += 1
            input_samples[sample_name] = {'index':index_name, 'bam':bam_file, 'vcf':vcf_name, 'vcf_idx':vcf_idx_name }

    flag = 0
    if flag==0:
        # figure out which samples need to be run
        incomplete = []
        for sample in input_samples.keys():
            if os.path.exists(input_samples[sample]['vcf_idx']):
                completed_samples += 1
            else:
                incomplete.append(input_samples[sample]['bam'])

        if completed_samples == input_samples_count:
            ## all samples are complete
            flag = 1
        else:
            ## samples still need to be analyzed
            ## link bam files to temporary directory for input dataset
            print "Samples Complete: %s" % (completed_samples)
            print "Samples Incomplete: %s" % (input_samples_count-completed_samples)
            print "Submitting Incomplete Samples: %s" % (incomplete)
            tmpdir = tempfile.mkdtemp(dir=os.getcwd())
            print tmpdir
            for bam in incomplete:
                sample_name = os.path.basename(bam)
                os.symlink(bam, "%s/%s" % (tmpdir, sample_name))

            gatk_cmd = " ".join(options.gatk_options)
            gatk_cmd = gatk_cmd.replace('"', '\\"')
            swift_params = list()
            swift_params.append('-bam_dir=' + options.input_dir) 
            #swift_params.append('-bam_dir=' + tmpdir)
            swift_params.append('-outputdir=' + options.output_dir)

            if not os.path.isdir(options.output_dir):
                os.mkdir(options.output_dir)

            ## check ruby version
            ruby_cmd = "ruby -v"
            proc = subprocess.Popen( args=ruby_cmd, shell=True )
            returncode = proc.wait()

            ## construct the swift command
            swift_cmd = "%s -sites.file %s -tc.file %s %s " %   (swift_bin, sites_file, tc_file, swift_file)
            cmd = "%s %s %s 2>&1" % (swift_cmd, ' '.join(swift_params), '-gatk_cmd=\"'+gatk_cmd+'\"')
            print cmd
            proc = subprocess.Popen( args=cmd, shell=True, cwd=tmp_dir )
            returncode = proc.wait()

            stdout = open( options.output_log, 'w' )
            stdout.write("BAM\t%s\nVCF\t%s\nREF\t%s\n" % (options.input_dir, options.output_dir, "FASTA"))
            stdout.close

            swift_log_files = glob.glob("%s/*.log" % tmp_dir)
            cmdSummary = "export PYTHONPATH=/mnt/galaxyTools/tools/pymodules/python2.7/lib/python:/opt/galaxy/lib:\$PYTHONPATH;. /mnt/galaxyTools/tools/pymodules/python2.7/env.sh;/opt/galaxy/tools/swift/parse_swift_log.py "
            for logF in swift_log_files:
                if "swift.log" in logF:
                    continue
                cmdSummary += " -l %s " % logF
            cmdSummary += " -o %s" % options.output_log

            return_code = None
            stderr = tempfile.NamedTemporaryFile( prefix="TOOL-stderr-", dir=tmp_dir )
            stdout = tempfile.NamedTemporaryFile( prefix="TOOL-stdout-", dir=tmp_dir )
            if return_code is None or not return_code:
                proc = subprocess.Popen( args=cmdSummary, stdout=stdout, stderr=stderr, shell=True, cwd=tmp_dir )
                return_code = proc.wait()
                if return_code:
                    stderr_target = sys.stderr
                else:
                    if stdout:
                        stderr_target = stdout
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

            print "\n\nEnd time:"
            print time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))

if __name__=="__main__": __main__()
