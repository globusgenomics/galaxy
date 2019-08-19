#!/usr/bin/python

import optparse, os, shutil, sys, tempfile, glob, shlex, vcf, pysam, gzip
from subprocess import *
import subprocess

CHUNK_SIZE = 2**20 #1mb

def createOutputHTML (outputF, sampleNames):
    ofh = open(outputF, "w")
    ofh.write( '<html>\n<head>\n<title>Galaxy - Manta Output</title>\n</head>\n<body>\n<p/>\n<ul>\n' )
    outputDir='%s_files' % ''.join(outputF.split('.')[:-1])
    for sample in sampleNames:
        values = sample.split(" ")
        sn = values[0]
        outVCF = "%s/%s.vcf" % (outputDir, sn)
        print "\noutVCF: %s\n" % outVCF
#        if os.path.exists(outVCF):
        ofh.write('<li><a href="%s">%s</a></li>\n' % ( outVCF, sn ) )
    ofh.write( '</ul>\n</body>\n</html>\n' )
    ofh.close()

def open_file_from_option( filename, mode = 'rb' ):
    if filename:
        return open( filename, mode = mode )
    return None

def create_final_output_files(outputDir, finalOutputDir):
    for x in os.listdir(outputDir):
        base=x
        tmpfile="%s/%s/results/variants/candidateSV.vcf.gz" % (outputDir, base)
        dest_file="%s/%s.vcf" % (finalOutputDir, base)
        with gzip.open(tmpfile, 'rb') as infile:
            with open(dest_file, 'w') as outfile:
                for line in infile:
                    outfile.write(line)

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
    parser.add_option( '-c', '--config', dest='config_file', action="store", type="string", default=None )
    parser.add_option( '-r', '--reference_path', dest='reference_path', help="reference file" )
    (options, args) = parser.parse_args()

    swift_bin = 'swift'
    sites_file = '/opt/galaxy/tools/manta/manta_sites.xml'
    tc_file = '/opt/galaxy/tools/swift/tc.data'
    swift_file = '/opt/galaxy/tools/manta/manta.swift'
    samtools_bin = "/mnt/galaxyTools/tools/samtools/1.2/bin/samtools"
    manta_bin = "/mnt/galaxyTools/tools/manta/manta-0.29.2/bin"

    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)
    output_dir = "%s/output" % options.output_dir
    config_dir = "%s/config" % options.output_dir
    inputDirectory = "%s/bams" % options.output_dir
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
        os.mkdir(config_dir)
        os.mkdir(inputDirectory)

    tmp_dir = tempfile.mkdtemp( dir=options.output_dir , prefix='tmp-TOOL-' )

    if options.input_dir_file:
        infile = open(options.input_dir_file, 'r')
        inputDirLine = infile.readline()
        inputRealDirectory = inputDirLine.rstrip('\n\r')
    else:
        inputRealDirectory = options.input_dir

    # create a link of the BAMs inside the inputs directory
    # this is to prevent the Warning of the BAI files being older than BAMs
    for inputF in glob.glob("%s/*.bam" % inputRealDirectory ):
        os.symlink(inputF, "%s/%s" % (inputDirectory, os.path.basename(inputF)))

    for inputF in glob.glob("%s/*.bai" % inputRealDirectory ):
        os.symlink(inputF, "%s/%s" % (inputDirectory, os.path.basename(inputF)))

    # get the input BAMs into a list
    inputFiles = glob.glob("%s/*.bam" % inputDirectory )
    print inputFiles
    sampleNames = []
    configFiles = []

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
   
        fh_base = open("%s/%s.txt" % (config_dir, sn), 'w')
        fh_base.write("%s\t150\t%s\n" % (inputF, sn))
        fh_base.close()
        sampleNames.append("%s" % (sn))
        configFiles.append("%s/%s.txt" % (config_dir, sn))
        # create the output directory for the sample
        sample_outD = "%s/%s" % (output_dir, sn) 
        os.mkdir(sample_outD)
        cmd = "configManta.py %s --bam %s --runDir %s --config %s" % (ptc, inputF, sample_outD, options.config_file)
        run_cmd(cmd, "Running config for sample: %s" % sn)

    # prepare tool command
    tool_cmd = "export PATH=%s:\$PATH;%s/SAMPLEID/runWorkflow.py -m local -j 2 --quiet " % (manta_bin, output_dir)

    #if no stderr file is specified, we'll use our own
    stderr = tempfile.NamedTemporaryFile( prefix="TOOL-stderr-", dir=tmp_dir )
    stdout = tempfile.NamedTemporaryFile( prefix="TOOL-stdout-", dir=tmp_dir )

    # prepare command line
    swift_params = list()
    swift_params.append('-outputdir=' + output_dir)
    swift_params.append('-samplenames=\"%s\"' %  ",".join(sampleNames))
 
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

    # create list output files in the HTML output
    try:
        #for sample in sampleNames:
        #    values = sample.split(" ")
        #    sn = values[0]
        #    sample_outputdir = "%s/%s" % (output_dir, sn)

        #    vcfFile = None
        #    gainLossFile = None
        #    for output in glob.glob("%s/table/*" % sample_outputdir):
        #        #print "OUTPUT: %s" % output
        #        if output.endswith("vcf"):
        #            vcfFile = output
        #        elif output.endswith("DetailsFILTERED.txt"):
        #            gainLossFile = output
 
        #    outputVCF = "%s/%s.vcf" % (options.output_dir, sn)
        #    modifyVCF(vcfFile, gainLossFile, outputVCF, sn)

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
