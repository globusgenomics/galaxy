#!/usr/bin/env python

"""
A wrapper script for running swift on globus-galaxy
"""

## MODULES
import sys, optparse, os, tempfile, subprocess, shutil

## GLOBALS    
swift_bin = '/mnt/galaxyTools/tools/swift/swift-0.94.1/bin/swift'
replace_bin = "/opt/galaxy/tools/protk/tpp_prophets-c04896f31ff7/replace_spectra_bug.sh"
search_perl_bin = "perl /opt/galaxy/tools/protk/tpp_prophets-c04896f31ff7/search.pl"
sites_file = '/opt/galaxy/tools/protk/protk-255b5b6ec617/sites.xml'
tc_file = '/opt/galaxy/tools/protk/protk-255b5b6ec617/tc.data'
swift_file = '/opt/galaxy/tools/protk/tpp_prophets-c04896f31ff7/spectra_search.swift'


#Parse Command Line     
parser = optparse.OptionParser()
parser.add_option('--input-dir', dest='input_dir', type="string", help='Input directory where all mzml file are located' )
parser.add_option('--output-dir', dest='outdir', type="string", help='Output directory to store output files' )
parser.add_option('-o', '--output', dest='output_file', type="string", help='Output file to store output files' )
parser.add_option( '', '--splibType', dest="splibType", type="string", help="splibType" )
parser.add_option( '', '--LibraryFile', dest="LibraryFile", type="string", help="LibraryFile" )
parser.add_option( '', '--search_database', dest="search_database", type="string", help="search_database" )
parser.add_option( '', '--queryType', dest="queryType", type="string", help="queryType" )
(options, args) = parser.parse_args()

## construct a command line template and have swift fill in the "parallelizalbe"
##  parameters using regexps.

# get the path to the mzml directory
mzmlfiles = options.input_dir.split()

# create the output directory
config_files = []
config_dir = None
if not os.path.exists(options.outdir):
    os.mkdir(options.outdir)
    # create the config_files for each mzml file in a specified directory
    config_dir = "%s/config" % options.outdir
    os.mkdir(config_dir)
    for mzml in mzmlfiles:
        name = os.path.basename(mzml)
        ext = name.split(".")[-1]
        config_name = "%s/%s.txt" % (config_dir, name)
        config_files.append(config_name)
        fh = open(config_name, "w")
        fh.write("*+* searchParams.splibType=\"%s\" *=-*\n" % options.splibType)
        fh.write("*+* searchParams.LibraryFile=\"%s\"*-*\n" % options.LibraryFile)
        fh.write("*+* searchParams.search_database=\"%s\"*-*\n" % options.search_database)
        fh.write("*+* searchParams.queryType=\"%s\" *=-*\n" % options.queryType)
        fh.write("*+* searchParams.inputType=\"%s\" *-*\n" % ext)
        fh.write("*+* searchParams.queryFile=\"%s\" *-*\n" % mzml)
        fh.write("*+* searchParams.queryFileName=\"%s\" *-*\n" % name)
        fh.close()
else:
    config_dir = "%s/config" % options.outdir

## parse mandatory commands
swift_params = list()
swift_params.append('-mzml_dir=' + config_dir)
swift_params.append('-outputdir=' + options.outdir)


## parse the optional commands
## These are parsed explicitly for clarity -- atlas command flags are opaque!
tool_cmd = 'export HOME=/home/galaxy; source /mnt/galaxyTools/tools/protkgem/protkgem-1.0.0/env.sh; source /mnt/galaxyTools/tools/tpp/TPP-4.6.2/env.sh; export RUBYLIB="$RUBYLIB:/mnt/galaxyTools/tools/protkgem/protkgem-1.0.0/ftools-0.0.0/lib:/mnt/galaxyTools/tools/protkgem/protkgem-1.0.0/protk-1.2.3/lib:/mnt/galaxyTools/tools/protkgem/protkgem-1.0.0/open4-1.3.0/lib:/mnt/galaxyTools/tools/protkgem/protkgem-1.0.0/bio-1.4.3/lib:/mnt/galaxyTools/tools/protkgem/protkgem-1.0.0/libxml-ruby-2.6.0/lib:/mnt/galaxyTools/tools/protkgem/protkgem-1.0.0/rest-client-1.6.7/lib:/mnt/galaxyTools/tools/protkgem/protkgem-1.0.0/mascot-dat-0.3.1/lib:/mnt/galaxyTools/tools/protkgem/protkgem-1.0.0/net-ftp-list-3.2.5/lib:/mnt/galaxyTools/tools/protkgem/protkgem-1.0.0/ruby-ole-1.2.11.6/lib:/mnt/galaxyTools/tools/protkgem/protkgem-1.0.0/mime-types-1.23/lib:/mnt/galaxyTools/tools/protkgem/protkgem-1.0.0/spreadsheet-0.8.5/lib";export PATH=/mnt/galaxyTools/tools/tpp/TPP-4.6.2/bin:$PATH; %s INPUT OUTPUT; %s OUTPUT' % (search_perl_bin, replace_bin)

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

