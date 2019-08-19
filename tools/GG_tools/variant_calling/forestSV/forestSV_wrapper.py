#!/usr/bin/env python
"""
Submits forestSV on a whole Genome sample (full chromosomes) or a dataset of chromosomal BAM files of one sample 
usage: forestSV_wrapper.py [options]
            --inputBam $inputFormat.inputBam  ## if submitting one BAM file as input
            --outcall $BamCallsOutFile        ## if submitting one BAM file as input
            --outrdata $BamRdataOutFile       ## if submitting one BAM file as input
            --inputHtml $inputFormat.inputHtml              ## if submitting a dataset of BAM files
            --outHtml $outputHtmlFile                       ## if submitting a dataset of BAM files
            --files-path $outputHtmlFile.extra_files_path   ## if submitting a dataset of BAM files
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

def __main__():
    print "\n\nStart time:"
    print time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))

    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '', '--inputBam', dest='inputBam', help='The input BAM' )
    parser.add_option( '', '--inputHtml', dest='inputHtml', help='The input BAM dataset in HTML format' )
    parser.add_option( '', '--output-bed-html', dest='outputBedFile', help='the dataset output for the BED files' )
    parser.add_option( '', '--output-bed-path', dest='outputBedPath', help='the dataset output for the BED files' )
    parser.add_option( '', '--output-calls-rdata-html', dest='outputCallsRdataFile', help='the dataset output for the Rdata Calls files' )
    parser.add_option( '', '--output-calls-rdata-path', dest='outputCallsRdataPath', help='the dataset output for the Rdata Calls files' )
    parser.add_option( '', '--output-scores-rdata-html', dest='outputScoresRdataFile', help='the dataset output for the Rdata Scores files' )
    parser.add_option( '', '--output-scores-rdata-path', dest='outputScoresRdataPath', help='the dataset output for the Rdata Scores files' )
    parser.add_option( '', '--db_key', dest='db_key', help="hg18/hg19" )
    parser.add_option( '', '--output-bed-all', dest='bed_output', help='the concatenated output for all BED files' )
    ( options, args ) = parser.parse_args()

    if options.inputBam:

        # exit if input file empty
        if os.path.getsize( options.inputBam ) == 0:
            raise Exception, 'Initial BAM file empty'
 
        ## this option is not developed yet

    else:
        # input is an HTML file containing a dataset of chromosomal BAM files for a sample

        #print "Creating directory: " + args.files_path
        if not os.path.exists(options.outputBedPath):
            os.makedirs(options.outputBedPath)
        if not os.path.exists(options.outputCallsRdataPath):
            os.makedirs(options.outputCallsRdataPath)
        if not os.path.exists(options.outputScoresRdataPath):
            os.makedirs(options.outputScoresRdataPath)

        ## create tmp directory to be deleted at end of run
        ## submitting using Swift
        wd_tmpdir = tempfile.mkdtemp(prefix='forestSV_', dir=options.outputBedPath)  ## make tmpdir

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
        f.write('     <profile namespace="globus" key="jobsPerNode">32</profile>\n')
        f.write('     <profile namespace="globus" key="maxwalltime">1000:00:00</profile>\n')
        f.write('     <profile namespace="globus" key="maxnodes">1</profile>\n')
        f.write('     <profile namespace="globus" key="nodeGranularity">1</profile>\n')
        f.write('     <profile namespace="globus" key="slots">1000</profile>\n')
        f.write('   </pool>\n')
        f.write('</config>\n')
        f.close()

        ## generate a tc.data file for the job
        forestSV_bin = os.popen("which %s" % "forestSV").read().strip()
        rscript_bin = os.popen("which %s" % "Rscript").read().strip()
        bash_bin = os.popen("which %s" % "bash").read().strip()
        tc_file = '%s/tc.data' % (wd_tmpdir)
        f = open(tc_file,'w')
        f.write('condor\tbash\t%s\n' % (bash_bin))
        f.write('condor\tforestSV\t%s\n' % (forestSV_bin))
        f.write('condor\tRscript\t%s\n' % (rscript_bin))
        f.close()

        ## get input directory path for all chromosomal BAMs
        input_bam_dir = '%s/%s_files' % (os.path.dirname(options.inputHtml), ('.').join(os.path.basename(options.inputHtml).split('.')[:-1]))
        print input_bam_dir

        ## for each file in input directory, create an info.txt file in the tmp directory
        info_dir = '%s/%s' % (wd_tmpdir, 'info')
        os.mkdir(info_dir)
        bam_name = None
        for dirfile in os.listdir(input_bam_dir):
            if fnmatch.fnmatch(dirfile, '*.bam'):
                chrm = ('.').join(os.path.basename(dirfile).split('.')[:-1])
                if "MT" in chrm:
                    "print skipping chrm %s" % (chrm)
                elif "Y" in chrm:
                    "print skipping %s" % (chrm)
                else:
                    filename = "%s/%s.info.txt" % ( info_dir, chrm )
                    f = open(filename, 'w')
                    f.write('chr\tfilename\tchrlength\tbas\n')
                    f.write('%s\t%s/%s\t%s\t\n' % (chrm, input_bam_dir, dirfile, options.db_key))
                    f.close()
                    print dirfile
                    if bam_name is None:
                        bam_name = get_bam_name("%s/%s" % (input_bam_dir, dirfile))

        ## create the command line
        training_file = "/mnt/galaxyTools/tools/R/2.12/site-library/forestSV/data/rf_1KG_ILMN_BWA_HG19_v1.Rdata"
        env_file = "/mnt/galaxyTools/tools/R/2.12/env.sh"
        swift_params = list()
        swift_params.append('-input_dir=' + info_dir)

        app_cmd = "source %s; %s %s --infofile=info_txt --basename=%s/%s --forest=%s" % (env_file, rscript_bin, forestSV_bin, options.outputBedPath, bam_name, training_file)
        print app_cmd

        ## construct the swift command
        #swift_bin = '/mnt/galaxyTools/tools/swift/swift-0.94.1/bin/swift'
        swift_file = '/nfs/software/galaxy/tools/GG_tools/variant_calling/forestSV/forestSV_wrapper.swift'
        swift_cmd = "swift -sites.file %s -tc.file %s %s " %   (sites_file, tc_file, swift_file)
        cmd = "%s %s %s 2>&1" % (swift_cmd, ' '.join(swift_params), '-app_cmd=\''+app_cmd+'\'')
        print cmd

        stderr_name = tempfile.NamedTemporaryFile( prefix = "cmd_stderr" ).name
        proc = Popen( args=cmd, shell=True, stderr=open( stderr_name, 'wb' ) )
        exit_code = proc.wait()
        if exit_code:
            for line in open( stderr_name ):
                print >> sys.stderr, line
                os.unlink( stderr_name ) #clean up
                raise Exception( "Error running command " )
            os.unlink( stderr_name ) #clean up

        ## move the output files to the top level directory
        outfile_dir = "%s/%s" % (options.outputBedPath, bam_name)
        concat_bed_outfile = open(options.bed_output, 'w')
        for filename in glob.glob(r'%s/*.txt' % (outfile_dir)):
            shutil.copy(filename, options.outputBedPath)
            with open(filename) as infile:
                for outline in infile:
                    concat_bed_outfile.write(outline)
        concat_bed_outfile.close()

        for filename in glob.glob(r'%s/*_calls.Rdata' % (outfile_dir)):
            shutil.copy(filename, options.outputCallsRdataPath)

        for filename in glob.glob(r'%s/*_scores.Rdata' % (outfile_dir)):
            shutil.copy(filename, options.outputScoresRdataPath)

        ## clean up temporary directory
        shutil.rmtree(outfile_dir)
        shutil.rmtree(wd_tmpdir)

        # create an output file with the name of all files created and stored in the library, and name of library
        #Create HTML file
        outputs = {options.outputBedFile : options.outputBedPath, options.outputCallsRdataFile : options.outputCallsRdataPath, options.outputScoresRdataFile : options.outputScoresRdataPath }
        output_types = {options.outputBedFile : 'BED', options.outputCallsRdataFile : 'Rdata', options.outputScoresRdataFile : 'Rdata' }
        for fi, pa in outputs.iteritems():
            f = open(fi,'w')
            f.write(galhtmlprefix)
            f.write('SV_Caller: forestSV\n')
            f.write('SampleName: %s\n<hr>\n' % (bam_name))
            for dirfile in os.listdir(pa):
                f.write('%s&emsp;<a href="%s">%s</a><br>\n' % (output_types[fi], dirfile, dirfile))
            f.write(galhtmlpostfix)
            f.close()

    print "\n\nEnd time:"
    print time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))

if __name__=="__main__": __main__()

