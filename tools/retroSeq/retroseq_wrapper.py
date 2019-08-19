#!/usr/bin/env python
"""
Submits retroSeq call or discover stage
"""

import fnmatch, time, re, optparse, os, sys, tempfile, shutil
from subprocess import *
import subprocess

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

def create_dict_idx_from_fasta( fasta_filename, target_dir = None ):
    if target_dir is None:
        target_dir = os.getcwd()
    basename = os.path.basename(fasta_filename)
    link_filename = os.path.join( target_dir, "%s.fasta" % ( basename ) )
    os.symlink( fasta_filename, link_filename )

    # if input reference is a fasta file, create the dict and index file needed 
    stderr_picard = tempfile.NamedTemporaryFile( prefix = "picard_stderr" ).name
    #print stderr_picard
    dict_filename = os.path.join( target_dir, "%s.dict" % ( basename ) )
    picard_jar = os.path.join( os.environ['PICARD_JAR'], "CreateSequenceDictionary.jar")

    cmd1 = "java -jar %s R= %s O= %s" % (picard_jar, link_filename, dict_filename)
    cmd2 = "samtools faidx %s" % (link_filename)
    process = subprocess.Popen(cmd1, shell=True, stderr=open( stderr_picard, 'wb' ), cwd=target_dir)
    rval = process.wait()
    if rval:
        for line in open( stderr_picard ):
            print >> sys.stderr, line
            os.unlink( stderr_picard ) #clean up
            cleanup_before_exit( tmp_dir )
            raise Exception( "Error creating reference dict file" )
        os.unlink( stderr_picard ) #clean up

    stderr_sam = tempfile.NamedTemporaryFile( prefix = "sam_stderr" ).name
    process2 = subprocess.Popen(cmd2, shell=True, stderr=open( stderr_sam, 'wb' ), cwd=target_dir)
    rval2 = process2.wait()
    if rval2:
        for line in open( stderr_sam ):
            print >> sys.stderr, line
            os.unlink( stderr_sam ) #clean up
            cleanup_before_exit( tmp_dir )
            raise Exception( "Error creating reference index file" )
        os.unlink( stderr_picard ) #clean up        

def create_bam_idx(bamFile, target_dir = None):
    if target_dir is None:
        target_dir = os.getcwd()

    index_name = "%s/%s.bai" % (target_dir, os.path.basename(bamFile))
    stderr_name = tempfile.NamedTemporaryFile( prefix = "bam_stderr" ).name
    command = 'samtools index %s %s' % ( bamFile, index_name )
    print command
    proc = subprocess.Popen( args=command, shell=True, stderr=open( stderr_name, 'wb' ) )
    exit_code = proc.wait()
    if exit_code:
        for line in open( stderr_name ):
            print >> sys.stderr, line
        os.unlink( stderr_name ) #clean up
        #cleanup_before_exit( tmp_dir )
        raise Exception( "Error indexing BAM file" )
    os.unlink( stderr_name ) #clean up

def create_sites_file(sites_file, swiftwork_dir):
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
    f.write('     <profile namespace="globus" key="maxnodes">1</profile>\n')
    f.write('     <profile namespace="globus" key="nodeGranularity">1</profile>\n')
    f.write('     <profile namespace="globus" key="slots">1000</profile>\n')
    f.write('   </pool>\n')
    f.write('</config>\n')
    f.close()

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
    parser.add_option( '', '--discover', dest='discover', action="store_true", help='Run the retroSeq discover tool' )
    parser.add_option( '', '--call', dest='call', action="store_true", help='Run the retroSeq call tool' )
    parser.add_option( '', '--swift', dest='swift', action="store_true", help='Run the retroSeq call tool using Swift' )
    parser.add_option( '', '--bam', dest='inputBam', help='The input BAM whole genome file.' )
    parser.add_option( '', '--indexBam', dest='bamFileIndex', help='The BAM index file' )
    parser.add_option( '', '--output-bed', dest='outputFile', help='Output file for a single BAM input.' )
    parser.add_option( '', '--output-vcf', dest='outputFileVcf', help='Output file for a single BAM input.' )
    parser.add_option( '', '--output-candidates', dest='outputFileCandidates', help='Output file for a single BAM input.' )
    parser.add_option( '', '--align', action="store_true", dest='align', help="Do the computational expensive exonerate PE discordant mate alignment step" )
    parser.add_option( '', '--refTEsName', dest='refTEsName', help='TE type of reference elements' )
    parser.add_option( '', '--refTEs', dest='refTEs', help='BED file of reference elements' )
    parser.add_option( '', '--eref', dest='eref', help='Tab file with list of transposon types and the corresponding fasta file of reference sequences' )
    parser.add_option( '', '--extra-files-path', dest='extra_files_path', help="Directory which will contain output data." )
    parser.add_option( '-q', '--quality', dest='quality', help="Minimum mapping quality for a read mate that anchors the insertion call" )
    parser.add_option( '', '--id', dest='id_perc', help="Minimum percent ID for a match of a read to the transposon references" )
    parser.add_option( '', '--length', dest='length', help="Minimum length of a hit to the transposon references" )
    parser.add_option( '', '--rgs', dest='rgs', help="Comma separated list of readgroups to operate on" )
    parser.add_option( '', '--exd', dest='exd', help="Fofn of BED files of regions where discordant mates falling into will be excluded" )
    parser.add_option( '', '--input', dest='discover_file', help="A single output file from the PE discover stage" )
    parser.add_option( '', '--ref', dest='ref', help="Fasta of reference genome" )
    parser.add_option( '', '--hets', dest='hets', action="store_true", help='Call heterozygous insertions. Default is homozygous' )
    parser.add_option( '', '--filter', dest='filter_file', help="Tab file with TE type and BED file of reference elements. These will be filtered out from the calling" )
    parser.add_option( '', '--region', dest='region', help="Call a particular chromosome only OR region" )
    parser.add_option( '', '--reads', dest='reads', help="Minimum number of reads required to make a call" )
    parser.add_option( '', '--depth', dest='depth', help="Max average depth of a region to be considered for calling" )
    parser.add_option( '', '--ignoreRGs', dest='ignoreRGs', help="Read group names that should be ignored" )
    ( options, args ) = parser.parse_args()

    # exit if input file empty
    if os.path.getsize( options.inputBam ) == 0:
        raise Exception, 'Initial BAM file empty'

    #print "Creating directory: " + args.files_path
    if not os.path.exists(options.extra_files_path):
        os.makedirs(options.extra_files_path)

    tmp_dir = tempfile.mkdtemp(prefix="bamIndex", dir=options.extra_files_path)

    # create link of input BAM file to the output directory
    basename = os.path.basename(options.inputBam)
    link_filename = os.path.join( tmp_dir, "%s.bam" % ( basename ) )
    os.symlink( options.inputBam, link_filename )

    # create a BAM index file for calling stage (not needed for discover stage)
    if options.call:
        if options.bamFileIndex:
            # create a link to the index file in the same directory as the bam file
            bai_link_filename = os.path.join(tmp_dir, "%s.bam.bai" % ( basename ) )
            os.symlink( options.bamFileIndex, bai_link_filename)
        else:
            create_bam_idx(link_filename, target_dir = tmp_dir)

        # does dict and fasta index files exist? if not, create link to a temp directory and create index, dict files
        fasta_path_basename = ('.').join(options.ref.split('.')[:-1])
        base = os.path.basename(options.ref)
        if not os.path.isfile("%s.dict" % fasta_path_basename):
            create_dict_idx_from_fasta( options.ref, target_dir = tmp_dir)
            index_dir = tmp_dir
            fastapath = "%s/%s.fasta" % (index_dir, base)
        else:
            index_dir = os.path.realpath(options.ref)
            fastapath = options.ref

    ## create tmp directory to be deleted at end of run
    tool_type = None
    if options.discover:
        tool_type = "discover"
    else:
        tool_type = "call"

    wd_tmpdir = tempfile.mkdtemp(prefix='retroseq_%s_' % (tool_type), dir=options.extra_files_path)  ## make tmpdir
    os.chdir(wd_tmpdir)
   
    # create the temp input file containing the path to the bam file
    te_path = '%s/TE_input.txt' % wd_tmpdir
    fh_te = open(te_path, 'w')
    if options.refTEs:
        fh_te.write("%s\t%s\n" % (options.refTEsName, options.refTEs))
    else:
        fh_te.write("%s\t%s\n" % (options.refTEsName, options.filter_file))
    fh_te.close()
 
    ## create the command
    retroseq_params = list()

    align_param = None
    if options.align:
        retroseq_params.append("-align")
    if options.refTEs:
        retroseq_params.append("-refTEs " + te_path)
    if options.eref:
        retroseq_params.append("-eref " + options.eref)
    if options.quality:
        retroseq_params.append("-q " + options.quality)
    if options.id_perc:
        retroseq_params.append("-id " + options.id_perc)
    if options.length:
        retroseq_params.append("-length " + options.length)
    if options.rgs and len(options.rgs)>0:
        retroseq_params.append("-rgs \"" + options.rgs + "\"")
    if options.exd and len(options.exd)>0:
        retroseq_params.append("-exd \"" + options.exd + "\"")
    if options.discover_file:
        retroseq_params.append("-input " + options.discover_file)
    if options.ref:
        retroseq_params.append("-ref " + fastapath)
    if options.hets:
        retroseq_params.append("-hets")
    if options.filter_file:
        retroseq_params.append("-filter " + te_path)
    if options.region:
        retroseq_params.append("-region " + options.region)
    if options.reads:
        retroseq_params.append("-reads " + options.reads)
    if options.depth:
        retroseq_params.append("-depth " + options.depth)
    if options.ignoreRGs:
        retroseq_params.append("-ignoreRGs " + options.ignoreRGs)


    if options.swift:
        ## generate a config file for the job:
        retroseq_path = os.popen("which %s" % "retroseq.pl").read().strip()
        bedtools_path = os.popen("which %s" % "bedtools").read().strip()
        bedtools_env = "/mnt/galaxyTools/tools/bedtools/2.17.0/env.sh"
        samtools_env = "/mnt/galaxyTools/tools/samtools/0.1.19/env.sh"

        app_cmd = 'source %s; source %s; %s -%s -bam %s -output OUTPUT_NAME %s -region CHR' % (bedtools_env, samtools_env, retroseq_path, tool_type, link_filename, ' '.join(retroseq_params)) 
        print app_cmd

        ## Swift job preparation
        ## generate a sites.xml file for the job:
        swiftwork_dir = '%s/%s' % (wd_tmpdir, 'swiftwork')
        os.mkdir('%s' % (swiftwork_dir))

        sites_file = '%s/sites.xml' % (wd_tmpdir)
        create_sites_file(sites_file, swiftwork_dir)

        ## generate a tc.data file for the job
        bash_bin = os.popen("which %s" % "bash").read().strip()
        samtools_path = os.popen("which %s" % "samtools").read().strip()
        tc_file = '%s/tc.data' % (wd_tmpdir)
        f = open(tc_file,'w')
        f.write('condor\tbash\t%s\n' % (bash_bin))
        f.write('condor\tretroseq.pl\t%s\n' % (retroseq_path))
        f.write('condor\tsamtools\t%s\n' % (samtools_path))
        f.write('condor\tbedtools\t%s\n' % (bedtools_path))
        f.close()
        
        ## create the command line
        swift_params = list()
        swift_params.append('-input_bam=' + link_filename)
        swift_params.append('-output_dir=' + options.extra_files_path)

        swift_file = '/nfs/software/galaxy/tools/retroSeq/retroseq_call.swift'
        swift_cmd = "swift -sites.file %s -tc.file %s %s " %   (sites_file, tc_file, swift_file)
        cmd = "%s %s %s 2>&1" % (swift_cmd, ' '.join(swift_params), '-app_cmd=\''+app_cmd+'\'')
        print cmd
        print "Current Work Dir: %s\n" % (os.getcwd())

        run_cmd(cmd, "Run Swift job for denovo local assembling")

        ## merge all chromosomal VCF, BED and CANDIDATES BED files into one per sample
        bed_files = list()
        vcf_files = list()
        candidates_bed_files = list()
        dir_contents = os.listdir(options.extra_files_path)
        for outfile in dir_contents:
            src = "%s/%s" % (options.extra_files_path, outfile)
            if "PE.vcf.candidates" in outfile:
                candidates_bed_files.append(src)
            elif "PE.vcf.bed" in outfile:
                bed_files.append(src)
            elif "PE.vcf" in outfile:
                vcf_files.append(src)
        
        #VCF
        cmd_vcf_concat = "vcf-concat %s | vcf-sort > %s" % (" ".join(vcf_files), options.outputFileVcf)
        print cmd_vcf_concat
        run_cmd(cmd_vcf_concat, "Run vcf-concat")

        #BED
        cmd_bed_concat = "cat %s | sort -u > %s" % (" ".join(bed_files), options.outputFile)
        run_cmd(cmd_bed_concat, "Run bed concat")
        
        # CANDIDATES BED
        filter_file = "%s/candidates.filter.txt" % (options.extra_files_path)
        chr_file = "%s/candidates.chr.txt" % (options.extra_files_path)
        test_file = "%s/candidates.test.txt" % (options.extra_files_path)
        cmd_candidates_filter_concat = "cat %s | grep 'FILTER' | sort -u > %s" % (" ".join(candidates_bed_files), filter_file)
        run_cmd(cmd_candidates_filter_concat, "Run candidates filter concat")
 
        cmd_candidates_chr_concat = "cat %s | grep -v 'FILTER' | grep -v 'TEST:' | sort -u > %s" % (" ".join(candidates_bed_files), chr_file)
        run_cmd(cmd_candidates_chr_concat, "Run candidates chr concat")

        cmd_candidates_test_concat = "cat %s | grep 'TEST:' | sort -u > %s" % (" ".join(candidates_bed_files), test_file)
        run_cmd(cmd_candidates_test_concat, "Run candidates test concat")

        cmd_candidates_concat = "cat %s %s %s > %s" % (filter_file, chr_file, test_file, options.outputFileCandidates)
        run_cmd(cmd_candidates_concat, "Run candidates concat")

    else:
        cmd = 'retroseq.pl -%s -bam %s -output %s %s ' % (tool_type, link_filename, options.outputFile, ' '.join(retroseq_params) )
        print cmd
        try:
            ## run the command
            os.chdir(wd_tmpdir)
            print os.getcwd()

            tmp = tempfile.NamedTemporaryFile( dir=wd_tmpdir ).name
            stderr_name = open( tmp, 'wb' )
            proc = Popen( args=cmd, shell=True, stderr=stderr_name.fileno())
            exit_code = proc.wait()
            stderr_name.close()

            stderr_name = open( tmp, 'rb' )
            stderr = ''
            buffsize = 1048576
            try:
                while True:
                    stderr += stderr_name.read( buffsize )
                    if not stderr or len( stderr ) % buffsize != 0:
                        break
            except OverflowError:
                pass
            stderr_name.close()
            if exit_code != 0:
                raise Exception, stderr
        except Exception, e:
            stop_err( 'Error running retroseq. ' + str( e ) )

        ## rename output files to what it should be:
        if options.call:
            shutil.move('%s.PE.vcf' % options.outputFile, options.outputFileVcf)
            shutil.move('%s.PE.vcf.bed' % options.outputFile, options.outputFile)
            shutil.move('%s.PE.vcf.candidates' % options.outputFile, options.outputFileCandidates)

    ## clean up temporary directory
    #shutil.rmtree(options.extra_files_path)

    print "\n\nEnd time:"
    print time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))

if __name__=="__main__": __main__()

