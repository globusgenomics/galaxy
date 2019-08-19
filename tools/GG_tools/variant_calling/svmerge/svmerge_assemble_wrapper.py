#!/usr/bin/env python
"""
Submits svmerge assemble
usage: svmerge_assemble_wrapper.py [options]
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

def create_bam_idx(bamFile, target_dir = None):
    if target_dir is None:
        target_dir = os.getcwd()
    basename = os.path.basename(bamFile)
    link_filename = os.path.join( target_dir, "%s.fasta" % ( basename ) )
    os.symlink( bamFile, link_filename )

    index_name = tempfile.NamedTemporaryFile( prefix = "bam_index" ).name

    stderr_name = tempfile.NamedTemporaryFile( prefix = "bam_stderr" ).name
    command = 'samtools index %s %s' % ( link_filename, index_name )
    proc = subprocess.Popen( args=command, shell=True, stderr=open( stderr_name, 'wb' ) )
    exit_code = proc.wait()
    if exit_code:
        for line in open( stderr_name ):
            print >> sys.stderr, line
        os.unlink( stderr_name ) #clean up
        #cleanup_before_exit( tmp_dir )
        raise Exception( "Error indexing BAM file" )
    os.unlink( stderr_name ) #clean up
    return index_name

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
    parser.add_option( '-b', '--bam-input', dest='bamFile', help='The input BAM dataset' )
    parser.add_option( '', '--bam-index', dest='bamIndex', help='The input BAM index' )
    parser.add_option( '', '--input-calls-file', dest='inputCallsFile', help='The input Bed file with the target SV regions.' )
    parser.add_option( '', '--fastaRef', dest='refFile', help='The fasta input reference file' )
    parser.add_option( '', '--indexes-path', dest='refIndexPath', help='repository to find the reference file location' )
    #parser.add_option( '', '--joblimit', dest='joblimit', help='job limit' )
    parser.add_option( '', '--checkdone', dest='checkdone', help='check done flag' )
    parser.add_option( '', '--subseq', dest='subseq', help='subseq flag' )
    parser.add_option( '', '--output', dest='output_assembled', help='output assembled directory or file?' )
    parser.add_option( '', '--extra-files-path', dest='extra_files_path', help='any extra files directory' )
    parser.add_option( '', '--species', dest='species', help='sample species' )
    ( options, args ) = parser.parse_args()

    #print "Creating directory: " + args.files_path
    if not os.path.exists(options.extra_files_path):
        os.makedirs(options.extra_files_path)

    ## create tmp directory to be deleted at end of run
    wd_tmpdir = tempfile.mkdtemp(prefix='svmerge_', dir=options.extra_files_path)  ## make tmpdir

    # check if BAM file has index
    baiFile = None
    if options.bamIndex:
        baiFile = options.bamIndex
    else:
        baiFile = create_bam_idx(options.bamFile, target_dir = wd_tmpdir)

    sample_name = get_bam_name(options.bamFile)

    ## generate a config file for the job:
    config_file = '%s/config.txt' % (wd_tmpdir)
    svmerge_assembly_path = os.popen("which %s" % "runAssembly.pl").read().strip()
    svmerge_path = os.path.dirname(svmerge_assembly_path)
    svmerge_data_path = "%s/../data" % svmerge_path
    velveth_path = os.popen("which %s" % "velveth").read().strip()
    velvetg_path = os.popen("which %s" % "velvetg").read().strip()
    samtools_path = os.popen("which %s" % "samtools").read().strip()
    exonerate_path = os.popen("which %s" % "exonerate").read().strip()
    bedexe_path = os.popen("which %s" % "intersectBed").read().strip()

    f = open(config_file,'w')

    f.write('name=%s\n' % sample_name)
    f.write('project=%s\n' % sample_name)
    f.write('version=REL-01\n')
    f.write('svdir=%s\n' % wd_tmpdir)
    f.write('bedexe=%s\n' % bedexe_path)
    f.write('exedir=%s\n' % svmerge_path)
    f.write('projdir=%s\n' % wd_tmpdir)
    f.write('callerlist=pindel\n')
    f.write('defaultQueue=normal\n')
    f.write('chrOther=chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY\n')
    f.write('species=%s\n' % options.species)

    f.write('samtools=%s\n' % samtools_path)
    f.write('exonerateExe=%s\n' % exonerate_path)
    f.write('cov_cutoff=2\n')
    f.write('bam=%s\n' % options.bamFile)
    f.write('bai=%s\n' % baiFile)

    f.write('makeConfig=1\n')
    f.write('reffile=%s\n' % (options.refFile) )
    f.write('submatrix=%s/submat.txt\n' % (svmerge_data_path))
    f.write('checkdone=%s\n' % (options.checkdone))
    f.write('subseq=%s\n' % (options.subseq))
    f.write('outdir=%s\n' % wd_tmpdir)
    f.write('assemMin=1\n')
    f.write('assemMax=100\n')
    f.write('assemQueue=normal\n')
    f.write('assemMem=2000\n\n')
    f.write('velvet=1\n')
    f.write('velveth=%s\n' % (velveth_path))
    f.write('velvetg=%s\n' % (velvetg_path))
    f.write('hashlen=29\n')
    f.write('ins_len=220\n')
    f.write('exp_cov=35\n')
    f.write('v_cutoff=2\n')
    f.write('##Abyss parameters\n')
    f.write('abyss=0\n')
    f.write('abyss-pe=/path/to/abyss-pe\n')
    f.write('kmer=25\n')
    f.write('npairs=10\n')

    f.write('bamcheck=1\n')
    f.write('meanCov=42\n')
    f.write('zyg=het\n')
    f.write('offset=1\n')
    f.write('parseSplit=1\n')
    f.write('parseSplitLines=250\n')
    f.write('parseJoblimit=10\n')
    f.write('parseQueue=normal\n')
    f.write('parseMem=2000\n')
    f.close()

    calls_file = options.inputCallsFile

    ## create and run the command to produce the sv assembly jobs files
    os.chdir(wd_tmpdir)
    cmd1 = 'coord2config.pl -c %s -b %s' % (config_file, calls_file)
    print cmd1
    run_cmd(cmd1, "Prepare assemble jobs")

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
    tc_file = '%s/tc.data' % (wd_tmpdir)
    f = open(tc_file,'w')
    f.write('condor\tbash\t%s\n' % (bash_bin))
    f.write('condor\texonerate\t%s\n' % (exonerate_path))
    f.write('condor\tsvAssemble.pl\t%s/svAssemble.pl\n' % (svmerge_path))
    f.write('condor\tsamtools\t%s\n' % (samtools_path))
    f.write('condor\tvelveth\t%s\n' % (velveth_path))
    f.write('condor\tvelvetg\t%s\n' % (velvetg_path))
    f.close()

    ## now create and run the swift command for running the assembly jobs
    app_cmd = "%s/svAssemble.pl --svconfig %s; tar -cvzf job.tar.gz *; mv job.tar.gz %s/job.tar.gz" % (svmerge_path, samtools_path, wd_tmpdir)
    #print app_cmd

    ## create the command line
    swift_params = list()
    swift_params.append('-input_dir=' + wd_tmpdir)

    swift_file = '/nfs/software/galaxy/tools/GG_tools/variant_calling/svmerge/svmerge_assemble_wrapper.swift'
    swift_cmd = "swift -sites.file %s -tc.file %s %s " %   (sites_file, tc_file, swift_file)
    cmd = "%s %s %s 2>&1" % (swift_cmd, ' '.join(swift_params), '-app_cmd=\''+app_cmd+'\'')
    print cmd
    print "Current Work Dir: %s\n" % (os.getcwd())

    run_cmd(cmd, "Run Swift job for denovo local assembling")

    ## untar the job.tar.gz file
    cmd1 = "tar xvzf job.tar.gz"
    print cmd1
    run_cmd(cmd1, "Untar the local assembly output file")

    ## need to merge all outputs from the local denovo assembly
    ## create and run the command to produce the sv merge jobs files
    #os.chdir(wd_tmpdir)
    cmd2 = 'splitFile.sh %s 500' % (calls_file)
    print cmd2
    run_cmd(cmd2, "Split call files")

    ## Swift job preparation
    ## generate a sites.xml file for the job:
    swiftwork_dir = '%s/%s' % (wd_tmpdir, 'swiftwork_merge')
    os.mkdir('%s' % (swiftwork_dir))

    sites_file = '%s/merge_sites.xml' % (wd_tmpdir)
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

    ## now create and run the swift command for running the assembly jobs
    os.chdir(wd_tmpdir)
    app_cmd = "source %s/../env.sh; tar -xvzf tarredFile -C ./; %s/runParser.pl -s sub_merge -c %s > %s/outfile" % (svmerge_path, svmerge_path, config_file, wd_tmpdir)
    print app_cmd

    ## create the command line
    calls_file_name = os.path.basename(calls_file)
    swift_params = list()
    swift_params.append('-input_dir=' + wd_tmpdir)
    swift_params.append('-mapper_grep=' + calls_file_name)
    swift_params.append('-input_tar=job.tar.gz')

    swift_file = '/nfs/software/galaxy/tools/GG_tools/variant_calling/svmerge/svmerge_merger_wrapper.swift'
    swift_cmd = "swift -sites.file %s -tc.file %s %s " %   (sites_file, tc_file, swift_file)
    cmd = "%s %s %s 2>&1" % (swift_cmd, ' '.join(swift_params), '-app_cmd=\''+app_cmd+'\'')
    print cmd
    run_cmd(cmd, "Run swift Parser job")

    ######## Concatenate the merged jobs ( all concatenate_files )
    os.chdir(wd_tmpdir)
    cmd3 = "cat %s/outfile.* > %s/outfile" % (wd_tmpdir, wd_tmpdir)
    print cmd3
    run_cmd(cmd3, "Concatenate merged jobs")

    ########## parse the alignparse file
    cmd4 = "parseBoundary.pl -a %s/outfile -o %s/sv.final" % (wd_tmpdir, wd_tmpdir)
    print cmd4
    run_cmd(cmd4, "Parse alignparse file")

    ########## merge final calls
    cmd5 = "mergeFinalCalls.pl -f %s/sv.final.rank1.tab -c %s > %s/sv.final.rank1.merged.tab" % (wd_tmpdir, config_file, wd_tmpdir)
    print cmd5
    run_cmd(cmd5, "merge final calls")

    #### print final bed file
    cmd6 = "tab2bed.pl -f %s/sv.final.rank1.merged.tab -n %s" % (wd_tmpdir, sample_name)
    print cmd6
    run_cmd(cmd6, "Print fianl bed file")

    ## move the output files to the top level directory
    shutil.copy("%s/sv.final.rank1.merged.bed" % (wd_tmpdir), options.output_assembled)

    ## clean up temporary directory 
    ## To debug, comment this line so all swift files are kept alive
    shutil.rmtree(wd_tmpdir)

    print "\n\nEnd time:"
    print time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))


if __name__=="__main__": __main__()

