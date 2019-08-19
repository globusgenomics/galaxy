#!/usr/bin/env python
#
"""
A wrapper script for running swift on RepeatMasker for multiple FASTA files
Input must be an HTML file with the names of the FASTA files
"""

## MODULES
import time, sys, optparse, os, tempfile, subprocess, shutil

print "\n\nStart time:"
print time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))


## GLOBALS    
swift_bin = '/mnt/galaxyTools/tools/swift/swift-0.94.1/bin/swift'
swift_file = '/nfs/software/galaxy/tools/RepeatMasker/RepeatMasker_with_swift.swift'
repeatMasker_bin = '/mnt/galaxyTools/tools/RepeatMasker/4.0.3/RepeatMasker'

#Parse Command Line     
parser = optparse.OptionParser()
parser.add_option('--parallel', dest='num_threads', type="int", help='Parallel threads to run this in' )
parser.add_option('--nolow', dest='nolow', action="store_true", help='Does not mask low_complexity DNA or simple repeats' )
parser.add_option('--noint', dest='noint', action="store_true", help='Only masks low complex/simple repeats, no interspersed repeats' )
parser.add_option('--norna', dest='norna', action="store_true", help='Does not mask small RNA genes' )
parser.add_option('--species', dest='species', type="string", help='')
parser.add_option('--dir', dest='output_dir', type="string", help='extra files directory')
parser.add_option('--gc', dest='gc', type="int", help='Use GC depended matrices')
parser.add_option('--gccalc', dest='gccalc', action="store_true", help='Use GC depended matrices, automaticly' )
parser.add_option('--gff', dest='gff', type="string", help='gff output file' )
parser.add_option('--gff-extra', dest='gff_outdir', type="string", help='gff output file' )
parser.add_option('--html', dest='html', type="string", help='html output file' )
parser.add_option('--html-extra', dest='html_outdir', type="string", help='html output file' )
parser.add_option('--mask', dest='mask', type="string", help='mask output file' )
parser.add_option('--mask-extra', dest='mask_outdir', type="string", help='mask output file' )
parser.add_option('--summary', dest='summary', type="string", help='summary output file' )
parser.add_option('--summary-extra', dest='summary_outdir', type="string", help='summary output file' )
parser.add_option('--slow', dest='slow', action="store_true", help='slow_search' )
parser.add_option('--quick', dest='quick_search', action="store_true", help='quick_search' )
parser.add_option('--qq', dest='rush_search', action="store_true", help='rush_search' )
parser.add_option('--alu', dest='alu', action="store_true", help='Only masks Alus ' )
parser.add_option('--is_only', dest='is_only', action="store_true", help='Mask only E coli insertion elements' )
parser.add_option('--stdout', dest='job_output', help='Output file' )
parser.add_option('--query', dest='query', help='Input HTML file containing the links to FASTA input files' )
parser.add_option('--bed', dest='bed_output', help='Output file for concatenated BED file' )
(opts, args) = parser.parse_args()

## construct a command line template and have swift fill in the "parallelizalbe"
##  parameters using regexps.

## create the extra files directory which will also contain the temporary files and directories
## which will need to be deleted at the end of the script
os.mkdir('%s' % (opts.output_dir))            ## make extra files directory
wd_tmpdir = tempfile.mkdtemp(prefix='repeatMaskerTmp_', dir=opts.output_dir)  ## make tmpdir

## make swiftwork directory
swiftwork_dir = '%s/%s' % (wd_tmpdir, 'swiftwork')
os.mkdir('%s' % (swiftwork_dir))

## generate a sites.xml file for the job:
sites_file = '%s/sites.xml' % (wd_tmpdir)
f = open(sites_file,'w')
f.write('<config>\n')
f.write('   <pool handle="condor">\n')
f.write('     <execution provider="coaster" url="none" jobmanager="local:condor"/>\n')
f.write('     <gridftp url="local://localhost"/>\n')
f.write('     <workdirectory>%s</workdirectory>\n' % (swiftwork_dir))
f.write('     <profile namespace="karajan" key="jobThrottle">1000</profile>\n')
f.write('     <profile namespace="karajan" key="initialScore">10000</profile>\n')
f.write('     <profile namespace="globus" key="condor.+Tenant">"ci"</profile>\n')
f.write('     <!-- <profile namespace="globus" key="condor.+GlobusOnline">false</profile> -->\n')
f.write('     <profile namespace="globus" key="jobsPerNode">4</profile>\n')
f.write('     <profile namespace="globus" key="maxwalltime">1000:00:00</profile>\n')
f.write('     <profile namespace="globus" key="maxnodes">1</profile>\n')
f.write('     <profile namespace="globus" key="nodeGranularity">1</profile>\n')
f.write('     <profile namespace="globus" key="slots">1000</profile>\n')
f.write('   </pool>\n')
f.write('</config>\n')
f.close()

## generate a tc.data file for the job
tc_file = '%s/tc.data' % (wd_tmpdir)
f = open(tc_file,'w')
f.write('condor\tbash\t/bin/bash\n')
f.write('condor\tls\t/bin/ls\n')
f.write('condor\tRepeatMasker\t%s\n' % (repeatMasker_bin))
f.close()

## get directory where FASTA files are stored. This is the value of the html input file
## (i.e. /nfs/software/galaxy/database/files/001/dataset_32.dat => /nfs/software/galaxy/database/files/001/dataset_32_files)
dir_basename = os.path.basename(opts.query)
fasta_dir = "%s/%s_files" % (os.path.dirname(opts.query), dir_basename.split('.')[:-1][0])

swift_params = list()
swift_params.append('-input_dir=' + fasta_dir)
swift_params.append('-outputdir=' + opts.output_dir)

## These are parsed explicitly for nput parameters for RepeatMasker 
repeatMasker_params = list()
if opts.num_threads: repeatMasker_params.append(" -parallel %s" % opts.num_threads)
if opts.nolow: repeatMasker_params.append(" -nolow ")
if opts.noint: repeatMasker_params.append(" -noint ")
if opts.norna: repeatMasker_params.append(" -norna ")
if opts.species: repeatMasker_params.append(" -species " + opts.species)
if opts.output_dir: repeatMasker_params.append(" -dir " + opts.output_dir)
if opts.gc > 0: repeatMasker_params.append(" -gc %s" % opts.gc) 
if opts.gccalc: repeatMasker_params.append(" -gccalc ")
if opts.gff: repeatMasker_params.append(" -gff ")
if opts.html: repeatMasker_params.append(" -html ")
if opts.slow: repeatMasker_params.append(" -s ")
if opts.quick_search: repeatMasker_params.append(" -q ")
if opts.rush_search: repeatMasker_params.append(" -qq ")
if opts.alu: repeatMasker_params.append(" -alu ")
if opts.is_only: repeatMasker_params.append(" -is_only ")

## initial vcftools command template
app_cmd = "%s/RepeatMasker %s --input" % ( repeatMasker_bin, ' '.join(repeatMasker_params))

## construct the swift command
swift_cmd = "%s -sites.file %s -tc.file %s %s " %   (swift_bin, sites_file, tc_file, swift_file)
cmd = "%s %s %s 2>&1" % (swift_cmd, ' '.join(swift_params), '-app_cmd=\''+app_cmd+'\'')
print cmd

stderr_name = tempfile.NamedTemporaryFile( prefix = "consensus_stderr" ).name
proc = subprocess.Popen( args=cmd, shell=True, stderr=open( stderr_name, 'wb' ) )
returncode = proc.wait()
if returncode:
    for line in open( stderr_name ):
        print >> sys.stderr, line
        os.unlink( stderr_name ) #clean up
        raise Exception( "Error running RepeatMasker " )
    os.unlink( stderr_name ) #clean up


## clean up temporary directory
shutil.rmtree(wd_tmpdir)

## write final html file with links to output files
galhtmlprefix = """<?xml version="1.0" encoding="utf-8" ?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<meta name="generator" content="Galaxy tool output - see http://getgalaxy.org/" />
<title></title>
<link rel="stylesheet" href="/static/style/base.css" type="text/css" />
</head>
<body>
<div class="document">
"""
galhtmlpostfix = """</div></body></html>\n"""
##Create HTML file
f = open(opts.job_output,'w')
f.write(galhtmlprefix)
dir_contents = os.listdir("%s" % (opts.output_dir))
dir_contents.sort()
validOutFiles = list()

if opts.mask:
    f_mask = open(opts.mask, 'w')
    f_mask.write(galhtmlprefix)
    os.mkdir('%s' % (opts.mask_outdir)) 
if opts.gff:
    f_gff = open(opts.gff, 'w')
    f_gff.write(galhtmlprefix)
    os.mkdir('%s' % (opts.gff_outdir))
if opts.html:
    f_html = open(opts.html, 'w')
    f_html.write(galhtmlprefix)
    os.mkdir('%s' % (opts.html_outdir))
if opts.summary:
    f_summary = open(opts.summary, 'w')
    f_summary.write(galhtmlprefix)
    os.mkdir('%s' % (opts.summary_outdir))

for outfile in dir_contents:
    out_basename = os.path.basename(outfile)
    src = "%s/%s" % (opts.output_dir, outfile)
    if "fasta.masked" in outfile:
        if opts.mask:
            shutil.move(src, opts.mask_outdir)
            f_mask.write('<a href="%s">%s</a><br>\n' % (outfile, outfile))
        else:
            os.unlink(src)
    elif "fasta.out.gff" in outfile:
        if opts.gff:
            shutil.move(src, opts.gff_outdir)
            f_gff.write('<a href="%s">%s</a><br>\n' % (outfile, outfile))
        else:
            os.unlink(src)
    elif "fasta.out.html" in outfile:
        if opts.html:
            shutil.move(src, opts.html_outdir)
            f_html.write('<a href="%s">%s</a><br>\n' % (outfile, outfile))
        else:
            os.unlink(src)
    elif "fasta.tbl" in outfile:
        if opts.summary:
            shutil.move(src, opts.summary_outdir)
            f_summary.write('<a href="%s">%s</a><br>\n' % (outfile, outfile))
        else:
            os.unlink(src)
    elif "fasta.out" in outfile:
        f.write('<a href="%s">%s</a><br>\n' % (outfile, outfile))
        validOutFiles.append(src)
    else:
        os.unlink(src)
f.write(galhtmlpostfix)
f.close()

if opts.mask:
    f_mask.write(galhtmlpostfix)
    f_mask.close()
if opts.gff:
    f_gff.write(galhtmlpostfix)
    f_gff.close()
if opts.html:
    f_html.write(galhtmlpostfix)
    f_html.close()
if opts.summary:
    f_summary.write(galhtmlpostfix)
    f_summary.close()

## create the output BED file
bedOut = open(opts.bed_output, 'w')
for validFile in validOutFiles:
    print validFile
    with open(validFile, 'r') as f:
        lines = filter(None, (line.rstrip() for line in f))
        for line in lines:
           values = line.split()
           try:
               sw_score = int(values[0])
               perc_div = float(values[1])
               perc_del = float(values[2])
               perc_ins = float(values[3])
               query_seq = values[4]
               query_start = values[5]
               query_end = values[6]
               query_str = values[8]
               match_rep = values[9]
               rep_class = values[10]
               
               bedOut.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (query_seq, query_start, query_end, match_rep, sw_score, query_str, perc_div, perc_del, perc_ins, rep_class))
           except ValueError:
               pass  # it was a string, not an int.
bedOut.close()

print "\n\nEnd time:"
print time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))

