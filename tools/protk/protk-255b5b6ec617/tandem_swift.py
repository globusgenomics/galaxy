#!/usr/bin/env python

"""
A wrapper script for running swift on globus-galaxy
"""

## MODULES
import sys, optparse, os, tempfile, subprocess, shutil

## GLOBALS    
swift_bin = '/mnt/galaxyTools/tools/swift/swift-0.94.1/bin/swift'
tandem_bin = '/mnt/galaxyTools/tools/ruby/1.9.3-p488/bin/ruby /mnt/galaxyTools/tools/protkgem/protkgem-1.0.0/protk-1.2.3/bin/tandem_search.rb'
sites_file = '/opt/galaxy/tools/protk/protk-255b5b6ec617/sites.xml'
tc_file = '/opt/galaxy/tools/protk/protk-255b5b6ec617/tc.data'
swift_file = '/opt/galaxy/tools/protk/protk-255b5b6ec617/tandem_search.swift'


#Parse Command Line     
parser = optparse.OptionParser()
parser.add_option('--input-dir', dest='input_dir', type="string", help='Input directory where all mzml file are located' )
parser.add_option('--outdir', dest='outdir', type="string", help='Output directory to store output files' )
parser.add_option('-o', '--output', dest='output_file', type="string", help='Output file to store output files' )
parser.add_option( '-p', '--pass_through', dest='pass_through_options', action='append', type="string", help='These options are passed through directly to tandem, without any modification.' )
(options, args) = parser.parse_args()

## construct a command line template and have swift fill in the "parallelizalbe"
##  parameters using regexps.

# get the path to the mzml directory
mzmlfiles = options.input_dir.split()

# create the output directory
if not os.path.exists(options.outdir):
    os.mkdir(options.outdir)

## parse mandatory commands
swift_params = list()
swift_params.append('-mzml_dir=' + os.path.dirname(mzmlfiles[0]))
swift_params.append('-outputdir=' + options.outdir)


## parse the optional commands
## These are parsed explicitly for clarity -- atlas command flags are opaque!
if options.pass_through_options:
    tool_cmd = 'export HOME=/home/galaxy; source /mnt/galaxyTools/tools/protkgem/protkgem-1.0.0/env.sh; export RUBYLIB="$RUBYLIB:/mnt/galaxyTools/tools/protkgem/protkgem-1.0.0/ftools-0.0.0/lib:/mnt/galaxyTools/tools/protkgem/protkgem-1.0.0/protk-1.2.3/lib:/mnt/galaxyTools/tools/protkgem/protkgem-1.0.0/open4-1.3.0/lib:/mnt/galaxyTools/tools/protkgem/protkgem-1.0.0/bio-1.4.3/lib:/mnt/galaxyTools/tools/protkgem/protkgem-1.0.0/libxml-ruby-2.6.0/lib:/mnt/galaxyTools/tools/protkgem/protkgem-1.0.0/rest-client-1.6.7/lib:/mnt/galaxyTools/tools/protkgem/protkgem-1.0.0/mascot-dat-0.3.1/lib:/mnt/galaxyTools/tools/protkgem/protkgem-1.0.0/net-ftp-list-3.2.5/lib:/mnt/galaxyTools/tools/protkgem/protkgem-1.0.0/ruby-ole-1.2.11.6/lib:/mnt/galaxyTools/tools/protkgem/protkgem-1.0.0/mime-types-1.23/lib:/mnt/galaxyTools/tools/protkgem/protkgem-1.0.0/spreadsheet-0.8.5/lib"; ' + tandem_bin + ' ' + ' '.join( options.pass_through_options )
else:
    tool_cmd = ''

## initial atlas command template
tool_cmd += " OUTPUT INPUT"

## check ruby version
ruby_cmd = "ruby -v"
proc = subprocess.Popen( args=ruby_cmd, shell=True )
returncode = proc.wait()

## construct the swift command
swift_cmd = "%s -sites.file %s -tc.file %s %s " %   (swift_bin, sites_file, tc_file, swift_file)
cmd = "%s %s %s 2>&1" % (swift_cmd, ' '.join(swift_params), '-tool_cmd=\"'+tool_cmd+'\"')
print cmd
proc = subprocess.Popen( args=cmd, shell=True )
returncode = proc.wait()

stdout = open( options.output_file, 'w' )
stdout.write("MZML\t%s\nOUTPUT\t%s\n" % (os.path.dirname(mzmlfiles[0]), options.outdir))
stdout.close()

