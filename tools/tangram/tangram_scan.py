#!/usr/bin/env python
"""
Run tangram_scan with either swift or without
usage: tangram_scan.py [options]
"""

import subprocess
import glob, fnmatch, time, re, optparse, os, sys, tempfile, shutil

CHUNK_SIZE = 2**20 #1mb

"""
    tangram_scan.py
      --in $input
      --dir $lib_table.extra_files_path
      --lib-table $lib_table
      --hist $hist
      #if $advanced_select == "yes":
          #if $advanced.cf != 0.01:
              --cf $advanced.cf
          #end if
          #if $advanced.tr != 0.02:
              --tr $advanced.tr
          #end if
          #if $advanced.mq != 20:
              --mq $advanced.mq
          #end if
          #if $advanced.mf != 10000:
              --mf $advanced.mf
          #end if
      #end if
"""

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


def __main__():
    print "\n\nStart time:"
    print time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))

    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '', '--input', dest='bamFile', help='The input Mosaik BAM file' )
    #parser.add_option( '', '--dir', dest='workdir', help='The directory where the outptu will be temporaryly written to.' )
    parser.add_option( '', '--hist', dest='hist', help='The hist output.' )
    parser.add_option( '', '--lib-table', dest='lib_table_output', help='The output file.' )
    parser.add_option( '', '--cf-threshold', dest='cf', help='Threshold for normal read pair in the fragment length distribution.' )
    parser.add_option( '', '--tr-fragment', dest='tr', help='Trim rate for the fragment length distribution' )
    parser.add_option( '', '--mq-reads', dest='mq', help='Minimum mapping quality for a normal read pair' )
    parser.add_option( '', '--mf-normal', dest='mf', help='Minimum number of nomral fragments in a library' )
    parser.add_option( '', '--swift', dest='swift', action='store_true', help="use swift for parallelization")
    parser.add_option( '', '--input-dataset', dest='bamDataset', help="input Mosaik BAM datset" )
    parser.add_option( '', '--lib-table-extra', dest='lib_table_extra_files_path', help="" )
    parser.add_option( '', '--hist-extra', dest='hist_extra_files_path', help="" )
    ( options, args ) = parser.parse_args()

    ## create the output dir
    os.mkdir(options.hist_extra_files_path)
    os.mkdir(options.lib_table_extra_files_path)

    ## Create a temporary directory where output will be written
    wd_tmpdir = tempfile.mkdtemp( prefix='tmp-tangram-scan-', dir=options.hist_extra_files_path )
    output_tmpdir = "%s/scan_output" % wd_tmpdir 
    os.mkdir(output_tmpdir)

    ## create the "in" file containing a list of BAM files as input and cmd
    cmd = ""
    if options.swift:
        bamDir = '%s/%s_files' % (os.path.dirname(options.bamDataset), ('.').join(os.path.basename(options.bamDataset).split('.')[:-1]))
        # create a unique file for each bam with it's in.txt
        for bamfile in glob.glob('%s/*.bam' % bamDir):
           bam_basename = os.path.basename(bamfile)
           inFile = open("%s/%s.in.txt" % (wd_tmpdir, bam_basename), "w")
           inFile.write("%s\n" % bamfile )
           inFile.close()
        
        ## create the swfit command line
        tangram_scan_path = os.popen("which %s" % "tangram_scan").read().strip()
        cmd = '%s -in INFILE -dir OUTDIR ' % tangram_scan_path
    else:
        inFile = open("%s/in.txt" % wd_tmpdir, "w")
        inFile.write("%s\n" % options.bamFile)
        inFile.close()

        ## create the command line
        cmd = "tangram_scan -in %s -dir %s " % ("%s/in.txt" % wd_tmpdir,  output_tmpdir)

    if options.cf:
        cmd += "-cf %s " % (options.cf)
    if options.tr:
        cmd += "-tr %s " % (options.tr)
    if options.mq:
        cmd += "-mq %s " % (options.mq)
    if options.mf:
        cmd += "-mf %s " % (options.mf)

    if options.swift:
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
        f.write('     <profile namespace="globus" key="jobsPerNode">32</profile>\n')
        f.write('     <profile namespace="globus" key="maxwalltime">1000:00:00</profile>\n')
        f.write('     <profile namespace="globus" key="maxnodes">1</profile>\n')
        f.write('     <profile namespace="globus" key="nodeGranularity">1</profile>\n')
        f.write('     <profile namespace="globus" key="slots">1000</profile>\n')
        f.write('   </pool>\n')
        f.write('</config>\n')
        f.close()

        ## generate a tc.data file for the job
        bash_bin = os.popen("which %s" % "bash").read().strip()
        tangram_scan_path = os.popen("which %s" % "tangram_scan").read().strip()
        tc_file = '%s/tc.data' % (wd_tmpdir)
        f = open(tc_file,'w')
        f.write('condor\tbash\t%s\n' % (bash_bin))
        f.write('condor\ttangram_scan\t%s\n' % (tangram_scan_path))
        f.close()

        app_cmd = cmd

        swift_params = list()
        swift_params.append('-input_dir=' + wd_tmpdir)
        swift_params.append('-output_dir=' + output_tmpdir)

        swift_file = '/nfs/software/galaxy/tools/tangram/tangram_scan.swift'
        swift_cmd = "swift -sites.file %s -tc.file %s %s " %   (sites_file, tc_file, swift_file)
        cmd = "%s %s %s 2>&1" % (swift_cmd, ' '.join(swift_params), '-app_cmd=\''+app_cmd+'\'')
        

    print cmd
    #set up stderr output options
    stderr = tempfile.NamedTemporaryFile( prefix="tangram-stderr-", dir=wd_tmpdir )

    proc = subprocess.Popen( args=cmd, stderr=stderr, shell=True, cwd=wd_tmpdir )
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

    if options.swift:
        ## go through the output directory and move files accordingly
        #Create HTML file
        f_hist = open(options.hist,'w')
        f_hist.write(galhtmlprefix)
        f_lib_tab = open(options.lib_table_output, 'w')
        f_lib_tab.write(galhtmlprefix)
 
        dir_contents = os.listdir(output_tmpdir)
        dir_contents.sort()
        for subdir in dir_contents:
            src = "%s/%s" % (output_tmpdir, subdir)
            hist_file = '%s/hist.dat' % src
            lib_tab_file = '%s/lib_table.dat' % src

            chr_prefix = ""
            chr_matchObj = re.match( r'(.*).bam.sorted', subdir)
            if chr_matchObj:
                chr_prefix = chr_matchObj.group(1)
                dest_hist_file = '%s/%s.hist.dat' % (options.hist_extra_files_path, chr_prefix)
                dest_lib_tab_file = '%s/%s.lib_table.dat' % (options.lib_table_extra_files_path, chr_prefix)           
            else:
                dest_hist_file = '%s/hist.dat' % (options.hist_extra_files_path)
                dest_lib_tab_file =  '%s/lib_table.dat' % (options.lib_table_extra_files_path)

            # move files
            shutil.copy(hist_file, dest_hist_file)
            shutil.copy(lib_tab_file, dest_lib_tab_file)

            # write the html files
            f_hist.write('<a href="%s">%s</a><br>\n' % (dest_hist_file, dest_hist_file))
            f_lib_tab.write('<a href="%s">%s</a><br>\n' % (dest_lib_tab_file, dest_lib_tab_file))
        f_hist.write(galhtmlpostfix)
        f_hist.close()
        f_lib_tab.write(galhtmlpostfix)
        f_lib_tab.close()

    else:
        ## After successful run, move output files to output file locations
        shutil.copy('%s/lib_table.dat' % (output_tmpdir), options.lib_table_output)
        shutil.copy('%s/hist.dat' % (output_tmpdir), options.hist)

    ## clean up tmp directories
    shutil.rmtree(wd_tmpdir)

    print "\n\nEnd time:"
    print time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))

if __name__=="__main__": __main__()

