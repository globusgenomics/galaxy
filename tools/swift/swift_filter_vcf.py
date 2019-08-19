#!/usr/bin/env python
# vasa

"""
A wrapper script for running swift on globus-galaxy
"""

## MODULES
import sys, optparse, os, tempfile, subprocess, shutil

## GLOBALS    
swift_bin = '/mnt/galaxyTools/tools/swift/swift-0.94.1/bin/swift'
sites_file = '/opt/galaxy/tools/swift/sites.xml'
tc_file = '/opt/galaxy/tools/swift/tc.data'
swift_file = '/opt/galaxy/tools/swift/filter_vcf.swift'
vcftools_bin = '/mnt/galaxyTools/tools/vcftools/vcftools_0.1.11/bin'
vcflib_bin = '/mnt/galaxyTools/tools/vcflib/10.27.2016/bin'

#Parse Command Line     
parser = optparse.OptionParser()
parser.add_option('--vcf-dir', dest='vcf_dir', type="string", help='Input directory where all VCF file are located' )
parser.add_option('--config', dest='config_file', type="string", help='Config file with location of VCF files' )
parser.add_option('--output-dir', dest='outputdir', type="string", help='Output directory to store recoded VCF files' )
parser.add_option('--bed', dest='bedfile', type="string", help='Bed file' )

parser.add_option('--minQ', dest='minq', type="string", help='Minimum Quality Threshold.')
parser.add_option('--filters', dest='filters', action="append", type="string", help='')
parser.add_option( '-d', '--dummy', dest='dummy_out', type="string")

(opts, args) = parser.parse_args()

## construct a command line template and have swift fill in the "parallelizalbe"
##  parameters using regexps.
if not os.path.isdir(opts.outputdir):
    os.mkdir(opts.outputdir)

## parse mandatory commands
mapper_exec = None
vcf_directory = opts.vcf_dir
if opts.config_file:
   # open file and extract variables
   fh = open(opts.config_file)
   for line in fh:
      if "VCF\t" in line:
         mapper_exec = "/opt/galaxy/tools/swift/vcf_mapper.sh"
      elif "INDEL\t" in line:
         mapper_exec = "/opt/galaxy/tools/swift/indel_mapper.sh"

swift_params = list()
swift_params.append('-vcf_dir=' + vcf_directory)
swift_params.append('-outputdir=' + opts.outputdir)
swift_params.append('-mapper_exec=' + mapper_exec)

## These are parsed explicitly for vcftools 
vcftools_params = list()
if opts.minq: vcftools_params.append(" --minQ " + opts.minq)
if opts.bedfile != "None": vcftools_params.append(" --bed " + opts.bedfile)
vcftools_params.append(" --recode ")

## initial vcftools command template
vcftools_cmd = "%s/vcftools --vcf --out %s " % (vcftools_bin, ' '.join(vcftools_params))


## Construct the vcffixup command
vcffixup_cmd = "%s/vcffixup " % (vcflib_bin)

## Construct the vcffilter command
# replace __gt__ or __lt__ with > or <
new_filters = list()
for filter_string in opts.filters:
   if '__gt__' in filter_string:
      new_filters.append(filter_string.replace('__gt__', '>'))
   elif '__lt__' in filter_string:
      new_filters.append(filter_string.replace('__lt__', '<'))
   else:
      new_filters.append(filter_string)
vcffilter_cmd = "%s/vcffilter -f \"%s\" --vcffilter_file " % (vcflib_bin, " & ".join(new_filters))

filter_cmd = "%s; %s; %s" % (vcftools_cmd, vcffixup_cmd, vcffilter_cmd)

## construct the swift command
swift_cmd = "%s -sites.file %s -tc.file %s %s " %   (swift_bin, sites_file, tc_file, swift_file)
cmd = "%s %s %s 2>&1" % (swift_cmd, ' '.join(swift_params), '-filter_cmd=\''+filter_cmd+'\'')
print cmd
proc = subprocess.Popen( args=cmd, shell=True )
returncode = proc.wait()


stdout = open( opts.dummy_out, 'w' )
stdout.write("INPUT\t%s\nVCF\t%s\n" % (vcf_directory, opts.outputdir))
stdout.close

"""
swift_log_files = glob.glob("%s/*.log" % tmp_dir)
cmdSummary = "/opt/galaxy/tools/swift/parse_swift_log.py "
for logF in swift_log_files:
    if "swift.log" in logF:
        continue
    cmdSummary += " -l %s " % logF
cmdSummary += " -o %s" % options.swift_log

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
"""
