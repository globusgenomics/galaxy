#!/usr/bin/python

import optparse, os, shutil, sys, tempfile, glob, shlex, vcf, pysam
from subprocess import *
import subprocess

CHUNK_SIZE = 2**20 #1mb

def createOutputHTML (outputF, sampleNames, outputDir):
    ofh = open(outputF, "w")
    print outputF
    ofh.write( '<html>\n<head>\n<title>Galaxy - Delly Output</title>\n</head>\n<body>\n<p/>\n<ul>\n' )

    for sample in sampleNames:
        values = sample.split(" ")
        sn = values[0]
        outVCF = "%s/%s.bcf" % (outputDir, sn)
        if os.path.exists(outVCF):
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

def merge(outputs, outputFile):
#       vcfWriter = vcf.Writer(open(outputFile, 'w'), template)
        template = pysam.VariantFile(filename=outputs[0])
        vcfWriter = pysam.VariantFile(open(outputFile, 'w'), template)
        for output in outputs:
                vcfReader = pysam.VariantFile(filename=output)
                for record in vcfReader:
                        vcfWriter.write(record)
                        #vcfWriter.write_record(record)
        return 0

def __main__():
    parser = optparse.OptionParser()
    parser.add_option( '', '--input_dir', dest='input_dir', action='store', type="string", help='Input directory path of BAM files' )
    parser.add_option( '', '--input_dir_file', dest='input_dir_file', action='store', type="string", help='Input directory File containing path of BAM files' )
    parser.add_option( '', '--out-dir', dest='output_dir', action='store', type="string", default=None, help='If specified, the output directory for extra files.' )
    parser.add_option( '', '--log', dest='swift_log', action='store', type="string", default=None, help='swift summary output.' )
    parser.add_option( '', '--output', dest='outputF', action='store', type="string", default=None, help='mpileup output.' )
    parser.add_option( '-p', '--pass_through', dest='pass_through_options', action='append', type="string", help='These options are passed through directly to contra, without any modification.' )
    parser.add_option( '-v', dest='vcffile', help='input VCF/BCF file for re-genotyping', action='store', type='string')
    parser.add_option( '-t', dest='types', action='append', type="string", help='SV analysis type (DEL, DUP, INV, TRA), INS')
    (options, args) = parser.parse_args()

    swift_bin = 'swift'
    sites_file = '/opt/galaxy/tools/delly/delly2_sites.xml'
    tc_file = '/opt/galaxy/tools/swift/tc.data'
    swift_file = '/opt/galaxy/tools/delly/delly2_merge.swift'
    samtools_bin = "/mnt/galaxyTools/tools/samtools/1.2/bin/samtools" # for pysam
    delly2_bin = "/mnt/galaxyTools/tools/delly/v0.7.3/delly2"

    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)
    output_dir = "%s/output" % options.output_dir
    inputDirectory = "%s/inputs" % options.output_dir
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
        os.mkdir(inputDirectory)

    inputFiles_svTypes = {}
    for sv_type in options.types:
        inputFiles_svTypes[sv_type] = []
        if not os.path.exists("%s/%s" % (inputDirectory, sv_type)):
            os.mkdir("%s/%s" % (inputDirectory, sv_type))

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

        inputFiles = sorted(glob.glob("%s/*.bcf" % inputRealDirectory ))
    else:
        inputFiles = options.input_files.strip(" ").split(" ")

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
        for sv_type in options.types:
            if sv_type in inputF:
                inputFiles_svTypes[sv_type].append(inputF)
                continue

    # prepare tool command
    env_cmd = "export LD_LIBRARY_PATH=/mnt/galaxyTools/tools/samtools/1.2/lib/:/mnt/galaxyTools/tools/delly/v0.7.3/lib:\$LD_LIBRARY_PATH; export PATH=%s:%s:\$PATH;" % (samtools_bin, delly2_bin)
    tool_cmds = []
    for sv_type in options.types:
        sv_tool_cmd = "%s merge -t %s %s -o %s/%s.bcf -b 500 -r 0.5 %s" % (delly2_bin, sv_type, ptc, options.output_dir, sv_type, " ".join(sorted(inputFiles_svTypes[sv_type])))
        tool_cmds.append(sv_tool_cmd)
    tool_cmd = ";".join(tool_cmds)

    #if no stderr file is specified, we'll use our own
    stderr = tempfile.NamedTemporaryFile( prefix="TOOL-stderr-", dir=tmp_dir )
    stdout = tempfile.NamedTemporaryFile( prefix="TOOL-stdout-", dir=tmp_dir )

    # prepare command line
    swift_params = list()
    swift_params.append('-envcmd=\"%s\"' % env_cmd)
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
    # merge the output files
    for sv_type in options.types:
        bcfFiles = glob.glob("%s/*%s*.bcf" % (options.output_dir, sv_type))
        if len(bcfFiles) > 1:
            merge_cmd="bcftools concat %s -a -o %s/%s.bcf" % (" ".join(bcfFiles), options.output_dir, sv_type)
            subprocess.call(merge_cmd,shell=True)
           
             
    
    # create list output files in the HTML output
    try:
        createOutputHTML(options.outputF, sampleNames, options.output_dir)

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
