#!/usr/bin/env python

"""
A wrapper script for running swift on globus-galaxy
"""

## MODULES
import sys, optparse, os, tempfile, subprocess, shutil

## GLOBALS    
swift_bin = '/mnt/galaxyTools/tools/swift/swift-0.94.1/bin/swift'
sites_file = '/opt/galaxy/tools/swift/sites.xml'
tc_file = '/opt/galaxy/tools/swift/tc.data'
swift_file = '/opt/galaxy/tools/swift/sv_sort_vcf.swift'
vcftools_bin = '/mnt/galaxyTools/tools/vcftools/vcftools_0.1.11/bin'
vcflib_bin = '/mnt/galaxyTools/tools/vcflib/03.08.2016/bin'

#Parse Command Line     
parser = optparse.OptionParser()
parser.add_option('--vcf-dir', dest='vcf_dir', type="string", help='Input directory where all VCF file are located' )
parser.add_option('--config', dest='config_file', type="string", help='Config file with location of VCF files' )
parser.add_option('--output-dir', dest='outputdir', type="string", help='Output directory to store recoded VCF files' )
parser.add_option( '-d', '--dummy', dest='dummy_out', type="string")

(opts, args) = parser.parse_args()

## construct a command line template and have swift fill in the "parallelizalbe"
##  parameters using regexps.
if not os.path.isdir(opts.outputdir):
    os.mkdir(opts.outputdir)

tmp_dir = tempfile.mkdtemp( dir=opts.outputdir, prefix='tmp-TOOL-' )

## parse mandatory commands
mapper_exec = None
if opts.vcf_dir:
    vcf_directory = opts.vcf_dir
elif opts.config_file:
    vcf_directory = opts.config_file.split('.')[0]+'_files'
    print vcf_directory
   # open file and extract variables
#   fh = open(opts.config_file)
#   for line in fh:
#      if "VCF\t" in line:
#         line = line.rstrip('\n')
#         vcf_directory = line.split("\t")[1]
#         mapper_exec = "/opt/galaxy/tools/swift/vcf_mapper.sh"
#      elif "INDEL\t" in line:
#         line = line.rstrip('\n')
#         vcf_directory = line.split("\t")[1]
#         mapper_exec = "/opt/galaxy/tools/swift/indel_mapper.sh"
else:
   sys.exit()
mapper_exec = "/opt/galaxy/tools/swift/vcf_mapper.sh"
#print opts.vcf_dir
#print mapper_exec
swift_params = list()
swift_params.append('-vcf_dir=' + vcf_directory)
print "vcf_directory:%s" % vcf_directory
swift_params.append('-outputdir=' + opts.outputdir)

if mapper_exec:
    swift_params.append('-mapper_exec=' + mapper_exec)

## These are parsed explicitly for vcftools 

vcfsort_cmd = "%s/vcf-sort -c --vcf-sort_file" % vcftools_bin

## construct the swift command
swift_cmd = "%s -sites.file %s -tc.file %s %s " %   (swift_bin, sites_file, tc_file, swift_file)
cmd = "%s %s %s 2>&1" % (swift_cmd, ' '.join(swift_params), '-sort_cmd=\''+vcfsort_cmd+'\'')
print cmd
proc = subprocess.Popen( args=cmd, shell=True, cwd=tmp_dir)
returncode = proc.wait()


stdout = open( opts.dummy_out, 'w' )
vcf_output_dir = opts.dummy_out.split('.')[0]+'_files'
stdout.write("Input directory: %s\nOutput directory: %s\n" % (vcf_directory, vcf_output_dir))
stdout.close()
