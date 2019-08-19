#!/usr/bin/python

import optparse, os, shutil, sys, tempfile, glob, shlex, vcf, pysam
from subprocess import *
import subprocess

CHUNK_SIZE = 2**20 #1mb

def createOutputHTML (outputF, sampleNames):
    ofh = open(outputF, "w")
    ofh.write( '<html>\n<head>\n<title>Galaxy - Breakdancer Output</title>\n</head>\n<body>\n<p/>\n<ul>\n' )
    outputDir='%s_files' % ''.join(outputF.split('.')[:-1])
    for sample in sampleNames:
        values = sample.split(" ")
        sn = values[0]
        outVCF = "%s/%s.txt" % (outputDir, sn)
        print "\noutput: %s\n" % outVCF 
#        if os.path.exists(outVCF):
        ofh.write('<li><a href="%s">%s</a></li>\n' % ( outVCF, sn ) )
    ofh.write( '</ul>\n</body>\n</html>\n' )
    ofh.close()

def open_file_from_option( filename, mode = 'rb' ):
    if filename:
        return open( filename, mode = mode )
    return None

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
    parser = optparse.OptionParser()
    parser.add_option( '', '--input_dir', dest='input_dir', action='store', type="string", help='Input directory path of BAM files' )
    parser.add_option( '', '--input_dir_file', dest='input_dir_file', action='store', type="string", help='Input directory File containing path of BAM files' )
    parser.add_option( '', '--out-dir', dest='output_dir', action='store', type="string", default=None, help='If specified, the output directory for extra files.' )
    parser.add_option( '', '--log', dest='swift_log', action='store', type="string", default=None, help='swift summary output.' )
    parser.add_option( '', '--output', dest='outputF', action='store', type="string", default=None, help='mpileup output.' )
    parser.add_option( '-p', '--pass_through', dest='pass_through_options', action='append', type="string", help='These options are passed through directly to contra, without any modification.' )
    
    (options, args) = parser.parse_args()

    swift_bin = 'swift'
    sites_file = '/opt/galaxy/tools/breakdancer/breakdancer_sites.xml'
    tc_file = '/opt/galaxy/tools/swift/tc.data'
    swift_file = '/opt/galaxy/tools/breakdancer/breakdancer.swift'
    samtools_bin = "/mnt/galaxyTools/tools/samtools/1.2/bin/samtools" # for pysam
    breakdancer_bin = "/mnt/galaxyTools/tools/breakdancer/1.4.5/bin/breakdancer-max"
    bam2cfg = "/mnt/galaxyTools/tools/breakdancer/1.4.5/lib/breakdancer-max1.4.5-unstable-66-4e44b43/bam2cfg.pl"
    breakdancer2vcf = "/mnt/galaxyTools/tools/breakdancer/1.4.5/lib/breakdancer-max1.4.5-unstable-66-4e44b43/breakdancer2vcf.py"
    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)
    output_dir = "%s/output" % options.output_dir
    inputDirectory = "%s/bams" % options.output_dir
    
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
        os.mkdir(inputDirectory)

    tmp_dir = tempfile.mkdtemp( dir=options.output_dir , prefix='tmp-TOOL-' )
    inputLinkedFiles = []
    basenames = []
    inputFiles = []
    inputIndexes = []
    if not options.input_dir_file:
        if options.input_dir_file: # path of the bam files
            infile = open(options.input_dir_file, 'r')
            inputDirLine = infile.readline()
            inputRealDirectory = inputDirLine.rstrip('\n\r')
        elif options.input_dir:
            inputRealDirectory = options.input_dir

        inputFiles = sorted(glob.glob("%s/*.bam" % inputRealDirectory ))
        inputIndexes = sorted(glob.glob("%s/*.bai" % inputRealDirectory ))
        # create a link of the BAMs inside the inputs directory
        # this is to prevent the Warning of the BAI files being older than BAMs
        #for inputF in inputFiles:
        #    #os.symlink(inputF, "%s/%s" % (inputDirectory, os.path.basename(inputF)))
        #    inputLinkedFiles.append("%s/%s" % (inputDirectory, os.path.basename(inputF)))
    else:
        inputFiles = options.input_files.strip(" ").split(" ")
        inputIndexes = sorted(glob.glob("%s/*.bai" % inputRealDirectory ))
        #for inputF in inputFiles:
        #    #os.symlink(inputF, "%s/%s" % (inputDirectory, os.path.basename(inputF)))
        #    inputLinkedFiles.append("%s/%s" % (inputDirectory, os.path.basename(inputF)) )

    # get the input BAMs into a list
    sampleNames = []


    # pass through options
    ptc = ""
    if options.pass_through_options:
        ptc = ' '.join( options.pass_through_options )

    # get the sample name and bam groups if necessary 
    # for the bam file and store in list in the same order as the Input file list
    #### Run the configManta command for each input on the head node since it's a low cost job
    for inputF in inputFiles:
        samfile = pysam.AlignmentFile(inputF, 'rb')
        sn = samfile.header['RG'][0]['SM']
   
        sampleNames.append("%s" % (sn))
        # create the output directory for the sample
        sample_outD = "%s/%s" % (output_dir, sn) 
        os.mkdir(sample_outD)

        # create symlink for input file
        os.symlink(inputF, "%s/%s.bam" % (inputDirectory, sn))
        inputLinkedFiles.append("%s/%s.bam" % (inputDirectory, sn))

    # prepare tool command
    tool_cmd = "export LD_LIBRARY_PATH=/mnt/galaxyTools/tools/samtools/1.2/lib/; export PERL5LIB=/mnt/galaxyTools/tools/perlmodules/perl5/modules/Statistics-Descriptive-3.0612:/mnt/galaxyTools/tools/perlmodules/perl5/modules/GDTextUtil-0.86/blib/lib:/mnt/galaxyTools/tools/perlmodules/perl5/modules/GD-2.46/blib/arch/auto/GD:/mnt/galaxyTools/tools/perlmodules/perl5/modules/GDGraph-1.52/blib/lib:/mnt/galaxyTools/tools/perlmodules/perl5/modules/GD-2.46/blib/lib:/mnt/galaxyTools/tools/perlmodules/perl5/modules/GDGraph-histogram-1.1/blib/lib;export PATH=%s:%s:%s:\$PATH; perl %s INPUTFILE > CONFIG1; %s CONFIG2 %s > OUTPUT1; python %s -i OUTPUT2 -o OUTPUT3" % (samtools_bin, bam2cfg, breakdancer_bin, bam2cfg, breakdancer_bin, ptc, breakdancer2vcf)


    #if no stderr file is specified, we'll use our own
    stderr = tempfile.NamedTemporaryFile( prefix="TOOL-stderr-", dir=tmp_dir )
    stdout = tempfile.NamedTemporaryFile( prefix="TOOL-stdout-", dir=tmp_dir )

    # prepare command line
    swift_params = list()
    swift_params.append('-outputdir=' + output_dir)
    swift_params.append('-samplenames=\"%s\"' %  ",".join(sampleNames))
    swift_params.append('-inputfiles=\"%s\"' %  ",".join(inputLinkedFiles))
 
    ## construct the swift command
    swift_cmd = "%s -sites.file %s -tc.file %s %s " %   (swift_bin, sites_file, tc_file, swift_file)
    cmd = "%s %s %s" % (swift_cmd, ' '.join(swift_params), '-tool_cmd=\"'+tool_cmd+'\"')
    print cmd
    return_code = None

    if return_code is None or not return_code:
        proc = subprocess.Popen( args=cmd, stdout=stdout, stderr=stderr, shell=True, cwd=tmp_dir )
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

    # copy vcf files to output directory
    for vcffile in glob.glob("%s/*/*.txt" % output_dir):
        fh_vcf_read = open(vcffile, "r")
        line1 = fh_vcf_read.readline()
        if "Error: no bams files in config file" in line1:
            filename = os.path.basename(vcffile)
            out_fake = open("%s/%s" % (options.output_dir, filename), "w")
            out_fake.write("")
            out_fake.close()
            sys.stdout.write("LOW QUALITY BAM: %s/%s.bam\n" % (inputRealDirectory, filename))
        else:
            shutil.copy(vcffile, options.output_dir)
        

    # create list output files in the HTML output
    try:
        createOutputHTML(options.outputF, sampleNames)

    except Exception, e:
        sys.stdout.write("problem while generating final output " + str(e))
    
    #try:
    #    if os.path.exists(tmp_dir):
    #        shutil.rmtree(tmp_dir)
    #    #if os.path.exists(output_dir):
    #    #    shutil.rmtree(output_dir)
    #except:
    #    pass
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

if __name__=="__main__":
	__main__()
