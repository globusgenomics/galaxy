#!/usr/bin/python

import optparse, os, shutil, sys, tempfile, glob, shlex, vcf, pysam, gzip, itertools
from subprocess import *
import subprocess

CHUNK_SIZE = 2**20 #1mb

def create_final_output_files(outputDir, finalOutputDir):
    for x in os.listdir(outputDir):
        base=x
        tmpfile="%s/%s/breakseq.vcf.gz" % (outputDir, base)
        dest_file="%s/%s.vcf" % (finalOutputDir, base)
        with gzip.open(tmpfile, 'rb') as infile:
            with open(dest_file, 'w') as outfile:
                for line in infile:
                    outfile.write(line)

def createOutputHTML (outputF, sampleNames):
    ofh = open(outputF, "w")
    ofh.write( '<html>\n<head>\n<title>Galaxy - Breakseq VCF Output</title>\n</head>\n<body>\n<p/>\n<ul>\n' )
    outputDir='%s_files' % ''.join(outputF.split('.')[:-1])
    for sample in sampleNames:
        values = sample.split(" ")
        sn = values[0]
        outVCF = "%s/%s.vcf" % (outputDir, sn)
        print "\noutVCF: %s\n" % outVCF 
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
    parser.add_option( '', '--checkChr', dest='checkchr', action='store', type="string", default=None, help='If specified, the output directory for extra files.' ) 

    (options, args) = parser.parse_args()

    swift_bin = 'swift'
    sites_file = '/opt/galaxy/tools/breakseq/breakseq_sites.xml'
    tc_file = '/opt/galaxy/tools/swift/tc.data'
    swift_file = '/opt/galaxy/tools/breakseq/breakseq.swift'
    python_bin = '/mnt/galaxyTools/tools/pymodules/python2.7/bin'
    pythonpath = "/mnt/galaxyTools/tools/pymodules/python2.7/lib/python"
    samtools_bin = "/mnt/galaxyTools/tools/samtools/1.2/bin/samtools" # for pysam
    bwa_bin = '/mnt/galaxyTools/tools/bwa/0.7.12/bwa'
    breakseq_bin = "/mnt/galaxyTools/tools/pymodules/python2.7/bin/run_breakseq2.py"

    if options.checkchr == "false":
        bplib_bin = "/mnt/galaxyTools/tools/breakseq2/2.2/lib/breakseq2_bplib_20150129.fna"
    else:
        bplib_bin = "/mnt/galaxyTools/tools/breakseq2/2.2/lib/breakseq2_bplib_20150129_new.fna"

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

        print inputRealDirectory
        inputFiles = sorted(glob.glob("%s/*.bam" % inputRealDirectory ))
        inputIndexes = sorted(glob.glob("%s/*.bai" % inputRealDirectory ))
        print inputFiles
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
    #for inputbam, inputbamindex in itertools.product(inputFiles, inputIndexes):
    for inputbam, inputbamindex in zip(inputFiles, inputIndexes):
        samfile = pysam.AlignmentFile(inputbam, 'rb')
        sn = samfile.header['RG'][0]['SM']
   
        sampleNames.append("%s" % (sn))
        # create the output directory for the sample
        sample_outD = "%s/%s" % (output_dir, sn) 
        if not os.path.exists(sample_outD):
            os.mkdir(sample_outD)

        # create symlink for input file
            os.symlink(inputbam, "%s/%s.bam" % (inputDirectory, sn))
            os.symlink(inputbamindex, "%s/%s.bam.bai" % (inputDirectory, sn))
        #pysam.index("%s/%s.bam" % (inputDirectory, sn))
            inputLinkedFiles.append("%s/%s.bam" % (inputDirectory, sn))

    #for inputF in glob.glob("%s/*.bai" % inputRealDirectory ):
    #    basename=os.path.basename(inputF).split("_")[0]
    #    os.symlink(inputF, "%s/%s.bam.bai" % (inputDirectory, basename))
    # prepare tool command
    tool_cmd = "export LD_LIBRARY_PATH=/mnt/galaxyTools/tools/samtools/1.2/lib:\$LD_LIBRARY_PATH; export PATH=%s:%s:%s:\$PATH; export PYTHONPATH=%s:\$PYTHONPATH; %s --bams INPUTFILE --samtools %s --bwa %s --bplib %s --work OUTPUTDIR %s" % (samtools_bin, breakseq_bin, python_bin, pythonpath, breakseq_bin, samtools_bin, bwa_bin, bplib_bin, ptc)


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

    # create final output with sample name
    create_final_output_files(output_dir, options.output_dir)

    # copy vcf files to output directory
#    for vcffile in glob.glob("%s/*/*.vcf" % output_dir):
#        shutil.copy(vcffile, options.output_dir)

    # create list output files in the HTML output
    try:
        createOutputHTML(options.outputF, sampleNames)

    except Exception, e:
        sys.stdout.write("problem while generating final VCF " + str(e))
    
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
