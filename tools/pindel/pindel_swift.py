#!/usr/bin/python

import optparse, os, shutil, sys, tempfile, glob, shlex, vcf, pysam
from subprocess import *
import subprocess

CHUNK_SIZE = 2**20 #1mb

def createOutputHTML (outputF, sampleNames, outputDir):
    ofh = open(outputF, "w")
    print outputF
    ofh.write( '<html>\n<head>\n<title>Galaxy - Contra Output</title>\n</head>\n<body>\n<p/>\n<ul>\n' )

    for sample in sampleNames:
        values = sample.split(" ")
        sn = values[0]
        outVCF = "%s/%s.vcf" % (outputDir, sn)
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

def modifyVCF(vcfFile, gainLossFile, outputVCF, sampleName):
    fh_gl = open(gainLossFile, "r")
    fh_vcf = open(vcfFile, "r")
    fh_outvcf = open(outputVCF, "w")

    gain_loss_dict = {}
    header = fh_gl.readline()
    for line in fh_gl:
       values = line.split("\t")
       name = "%s-%s" % (values[3],values[4])
       gain_loss_dict[name] = values[13]
    fh_gl.close()

    for line in fh_vcf:
        if line.startswith("##"):
            if line.startswith("##ALT"):
                line += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                line += "##FORMAT=<ID=SV,Number=1,Type=String,Description=\"SV type\">\n"

        elif line.startswith("#CHROM"):
            line = line.replace (" ", "")
            line = line.rstrip("\n")
            line += "\t%s\t%s\n" % ("FORMAT", sampleName)
        else:
            line = line.rstrip("\n")
            values = line.split("\t")
            name = "%s-%s" % (values[0],values[1])
            line += "\t%s\t%s\n" % ("GT:SV", "0/1:%s" % (gain_loss_dict[name]))
        fh_outvcf.write(line)
    fh_outvcf.close()
    fh_vcf.close()

def lineCount(inFile):
    fh = open(inFile, "r")
    counter = 0
    for line in fh:
        counter += 1
    return counter

def create_additional_bam_copy(inFile, baseline_indir):
    fh = open(inFile, "r")
    fpath = None
    for line in fh:
        fpath = line.rstrip("\n")
    fh.close()
    name = os.path.basename(fpath)
    linkname = "%s/%s.dup.bam" % (baseline_indir, name)
    os.symlink(fpath, linkname)
    
    # add link to existing file
    fh = open(inFile, "a")
    fh.write(linkname)
    fh.close()

def __main__():
    parser = optparse.OptionParser()
    parser.add_option( '', '--input_dir', dest='input_dir', action='store', type="string", help='Input directory path of BAM files' )
    parser.add_option( '', '--input_dir_file', dest='input_dir_file', action='store', type="string", help='Input directory File containing path of BAM files' )
    parser.add_option( '', '--out-dir', dest='output_dir', action='store', type="string", default=None, help='If specified, the output directory for extra files.' )
    parser.add_option( '', '--stdout', dest='stdout', action='store', type="string", default=None, help='If specified, the output of stdout will be written to this file.' )
    parser.add_option( '', '--log', dest='swift_log', action='store', type="string", default=None, help='swift summary output.' )
    parser.add_option( '', '--output', dest='outputF', action='store', type="string", default=None, help='mpileup output.' )
    parser.add_option( '-p', '--pass_through', dest='pass_through_options', action='append', type="string", help='These options are passed through directly to contra, without any modification.' )
    parser.add_option( '-t', '--target', dest='bed_file', action="store", type="string", default=None )
    parser.add_option( '-n', '--sampleName', dest='sampleName', help='sampleName' )
    parser.add_option( '-o', '--contra-vcf', dest='outputVCF', help='outpuVCF' )
    parser.add_option( '-r', '--reference_path', dest='reference_path', help="reference file" )
    (options, args) = parser.parse_args()

    swift_bin = 'swift'
    sites_file = '/opt/galaxy/tools/pindel/pindel_sites.xml'
    tc_file = '/opt/galaxy/tools/swift/tc.data'
    swift_file = '/opt/galaxy/tools/pindel/pindel.swift'
    samtools_bin = "/mnt/galaxyTools/tools/samtools/1.2/bin/samtools"
    pindel_bin = "/mnt/galaxyTools/tools/pindel/0.2.5b7/pindel"
    pindel2vcf_bin = "/mnt/galaxyTools/tools/pindel/0.2.5b7/pindel2vcf"

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

    # get the sample name and bam groups if necessary 
    # for the bam file and store in list in the same order as the Input file list
    is_cmd = "/mnt/galaxyTools/tools/breakdancer/1.4.5/lib/breakdancer-max1.4.5-unstable-66-4e44b43/bam2cfg.pl %s | awk '{ print $9; }' | cut -d \":\" -f 2" % inputFiles[0]
    insert_size = os.popen(is_cmd).read().rstrip("\n")
    for inputF in inputFiles:
        samfile = pysam.AlignmentFile(inputF, 'rb')
        sn = samfile.header['RG'][0]['SM']
   
        fh_base = open("%s/%s.txt" % (config_dir, sn), 'w')
        fh_base.write("%s\t%s\t%s\n" % (inputF, insert_size, sn))
        fh_base.close()
        sampleNames.append("%s" % (sn))
        configFiles.append("%s/%s.txt" % (config_dir, sn))
        # create the output directory for the sample
        os.mkdir("%s/%s" % (output_dir, sn))

    #### Contra part of the tool command
    ptc = ""
    if options.pass_through_options:
        ptc = ' '.join( options.pass_through_options )

#    tool_cmd = "export LD_LIBRARY_PATH=/mnt/galaxyTools/tools/samtools/1.2/lib/; export PATH=%s:%s:\$PATH;%s -o OUTPUTDIR -i INPUTFILE -L SAMPLENAME %s; %s -P OUTPUTDIR2 -r %s -R UCSC_hg19 -d 20150519 -v OUTPUT_VCF -G; awk \'!x[$0]++\' OUTPUT_VCF2 > OUTPUT_VCF3 " % (samtools_bin, pindel_bin, pindel_bin, ptc, pindel2vcf_bin, options.reference_path)
    tool_cmd = "export LD_LIBRARY_PATH=/mnt/galaxyTools/tools/samtools/1.2/lib/; export PATH=%s:%s:\$PATH;%s -o OUTPUTDIR -i INPUTFILE -L SAMPLENAME %s; %s -P OUTPUTDIR2 -r %s -R UCSC_hg19 -d 20150519 -v OUTPUT_VCF -G " % (samtools_bin, pindel_bin, pindel_bin, ptc, pindel2vcf_bin, options.reference_path)

    #if no stderr file is specified, we'll use our own
    stderr = tempfile.NamedTemporaryFile( prefix="TOOL-stderr-", dir=tmp_dir )
    stdout = tempfile.NamedTemporaryFile( prefix="TOOL-stdout-", dir=tmp_dir )


    # prepare command line
    swift_params = list()
    swift_params.append('-outputdir=' + output_dir)
    swift_params.append('-inputfiles=\"%s\"' %  ",".join(configFiles))
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
    # final output
        

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
