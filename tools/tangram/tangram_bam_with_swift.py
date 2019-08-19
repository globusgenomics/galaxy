#!/usr/bin/env python
"""
Submits tangram_bam on a complete dataset and output a sorted BAM file
usage: tangram_bam_with_swift.py [options]
"""

import glob, fnmatch, time, re, optparse, os, sys, tempfile, shutil
from subprocess import *

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

def stop_err( msg ):
    sys.stderr.write( '%s\n' % msg )
    sys.exit()

def cleanup_before_exit( tmp_dir ):
    if tmp_dir and os.path.exists( tmp_dir ):
        shutil.rmtree( tmp_dir )

def get_bam_name( bamfile ):
    ## read the BAM file header to get the sample ID name
    cmd = "samtools view -H %s" % bamfile
    #print cmd
    proc = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE )
    for line in proc.stdout:
        #line = proc.stdout.readline()
        if "SM:" in line:
            # split by tabs
            for field in line.rsplit("\t"):
                if "SM:" in field:
                    matchObj = re.search( r'SM:(.*)', field)
                    if matchObj:
                        sampleName = matchObj.group(1)
                        return sampleName


def run_cmd ( cmd , descriptor):

    stderr_name = tempfile.NamedTemporaryFile( prefix = "cmd_stderr" ).name
    proc = Popen( args=cmd, shell=True, stderr=open( stderr_name, 'wb' ) )
    exit_code = proc.wait()
    if exit_code:
        for line in open( stderr_name ):
            print >> sys.stderr, line
            os.unlink( stderr_name ) #clean up
            raise Exception( "Error running command: %s " % descriptor )
        os.unlink( stderr_name ) #clean up


def __main__():
    print "\n\nStart time:"
    print time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))

    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '-i', '--input', dest='bamFile', help='The input BAM dataset' )
    parser.add_option( '-r', '--ref', dest='refFile', help='The input reference file' )
    parser.add_option( '-o', '--output', dest='output_html', help='output dataset html file' )
    parser.add_option( '', '--extra-files-path', dest='extra_files_path', help='any extra files directory' )
    ( options, args ) = parser.parse_args()

    #print "Creating directory: " + args.files_path
    if not os.path.exists(options.extra_files_path):
        os.makedirs(options.extra_files_path)

    ## create tmp directory to be deleted at end of run
    wd_tmpdir = tempfile.mkdtemp(prefix='tangram_', dir=options.extra_files_path)  ## make tmpdir
    print wd_tmpdir
    sample_name = get_bam_name(options.bamFile)

    ## generate a config file for the job:
    tangram_bam_path = os.popen("which %s" % "tangram_bam").read().strip()

    ## create and run the command to produce the MOSAIK formatted BAM files for the dataset 
    os.chdir(wd_tmpdir)

    ## Swift job preparation
    ## generate a sites.xml file for the job:
    swiftwork_dir = '%s/%s' % (wd_tmpdir, 'swiftwork')
    os.mkdir('%s' % (swiftwork_dir))

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
    f.write('     <profile namespace="globus" key="jobsPerNode">1</profile>\n')
    f.write('     <profile namespace="globus" key="maxwalltime">1000:00:00</profile>\n')
    f.write('     <profile namespace="globus" key="maxnodes">100</profile>\n')
    f.write('     <profile namespace="globus" key="nodeGranularity">1</profile>\n')
    f.write('     <profile namespace="globus" key="slots">1000</profile>\n')
    f.write('   </pool>\n')
    f.write('</config>\n')
    f.close()

    ## generate a tc.data file for the job
    bash_bin = os.popen("which %s" % "bash").read().strip()
    samtools_path = os.popen("which %s" % "samtools").read().strip()
    tc_file = '%s/tc.data' % (wd_tmpdir)
    f = open(tc_file,'w')
    f.write('condor\tbash\t%s\n' % (bash_bin))
    f.write('condor\ttangram_bam\t%s\n' % (tangram_bam_path))
    f.write('condor\tsamtools\t%s\n' % (samtools_path))
    f.close()

    ## now create and run the swift command for running the assembly jobs
    input_bam_dir = '%s/%s_files' % (os.path.dirname(options.bamFile), ('.').join(os.path.basename(options.bamFile).split('.')[:-1]))
    app_cmd = "%s -i inputBAM -r %s -o outputBAM; %s sort -m 16G outputBAM  outputSortedBAM;rm outputBAM" % (tangram_bam_path, options.refFile, samtools_path)
    #print app_cmd

    ## create the command line
    swift_params = list()
    swift_params.append('-execution.retries=5' )
    swift_params.append('-input_dir=' + input_bam_dir)
    swift_params.append('-output_dir=' + options.extra_files_path)

    swift_file = '/nfs/software/galaxy/tools/tangram/tangram_bam.swift'
    swift_cmd = "swift -sites.file %s -tc.file %s %s " %   (sites_file, tc_file, swift_file)
    cmd = "%s %s %s 2>&1" % (swift_cmd, ' '.join(swift_params), '-app_cmd=\''+app_cmd+'\'')
    print cmd
    print "Current Work Dir: %s\n" % (os.getcwd())

    run_cmd(cmd, "Run Swift job for denovo local assembling")

    ## clean up temporary directory 
    ## To debug, comment this line so all swift files are kept alive
    shutil.rmtree(wd_tmpdir)

    # create an output file with the name of all files created and stored in the library, and name of library
    #Create HTML file
    f = open(options.output_html,'w')
    f.write(galhtmlprefix)
    dir_contents = os.listdir("%s" % (options.extra_files_path))
    dir_contents.sort()

    for outfile in dir_contents:
        src = "%s/%s" % (options.extra_files_path, outfile)
        if "sorted.bam" in outfile:
            f.write('<a href="%s">%s</a><br>\n' % (outfile, outfile))
        else:
            os.unlink(src)
    f.write(galhtmlpostfix)
    f.close()

    print "\n\nEnd time:"
    print time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))


if __name__=="__main__": __main__()

