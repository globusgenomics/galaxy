#!/usr/bin/python

import optparse, os, shutil, sys, tempfile, glob, shlex, vcf, pysam
from subprocess import *
import subprocess, math, random

CHUNK_SIZE = 2**20 #1mb

def createOutputHTML (outputF, sampleNames):
    ofh = open(outputF, "w")
    ofh.write( '<html>\n<head>\n<title>Galaxy - Contra Output</title>\n</head>\n<body>\n<p/>\n<ul>\n' )
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

def modifyContraVCF(vcfFile, gainLossFile, outputVCF, sampleName):
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
    parser.add_option( '', '--output', dest='outputF', action='store', type="string", default=None, help='mpileup output.' )
    parser.add_option( '', '--log', dest='swift_log', action='store', type="string", default=None, help='swift summary output.' )
    parser.add_option( '-p', '--pass_through', dest='pass_through_options', action='append', type="string", help='These options are passed through directly to contra, without any modification.' )
    parser.add_option( '-t', '--target', dest='bed_file', action="store", type="string", default=None )
    parser.add_option( '-n', '--sampleName', dest='sampleName', help='sampleName' )
    parser.add_option( '-o', '--contra-vcf', dest='outputVCF', help='outpuVCF' )
    parser.add_option( '', '--trim', dest='trim', help='contra baseline: --trim' )
    parser.add_option( '', '--percent-bam-files-for-baseline', type="float", dest='percent', help='contra baseline group: percentage of BAMs to use' )
    parser.add_option( '', '--baseline-input-bam', dest='group_by_keyword', help='contra baseline group: to group or not to group' )
    parser.add_option( '', '--group-field', dest='rg_field', help='contra baseline group: RG field to use for grouping' )
    parser.add_option( '', '--keyword-separator', dest='field_separator', help='contra baseline group: RG field separator' )
    parser.add_option( '', '--keyword-field-order', dest='field_order', help='contra baseline group: RG field order' )
    (options, args) = parser.parse_args()

    swift_bin = 'swift'
    sites_file = '/opt/galaxy/tools/contra/contra_sites.xml'
    tc_file = '/opt/galaxy/tools/swift/tc.data'
    swift_file = '/opt/galaxy/tools/contra/contra.swift'
    swift_baseline_file = '/opt/galaxy/tools/contra/contra_baseline.swift'
    swift_contra_to_metasv_file= '/opt/galaxy/tools/contra/contra_to_metasv.swift'
    samtools_bin = "/mnt/galaxyTools/tools/samtools/1.1/bin/samtools"
    bedtools_bin = "/mnt/galaxyTools/tools/bedtools/2.17.0/bin"
    #R_bin = "/mnt/galaxyTools/tools/R/3.0.0"
    R_bin = "/mnt/galaxyTools/tools/R/3.2.2"
    pymodules_bin = "/mnt/galaxyTools/tools/pymodules/python2.7"
    contra_bin = "/mnt/galaxyTools/tools/contra/CONTRA.v2.0.6/bin"

    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)
    output_dir = "%s/output" % options.output_dir
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    if options.group_by_keyword is not None:
        baseline_dir = "%s/baseline" % options.output_dir
        if not os.path.exists(baseline_dir):
            os.mkdir(baseline_dir)
            os.mkdir("%s/inputs" % baseline_dir)
            os.mkdir("%s/outputs" % baseline_dir)

    tmp_dir = tempfile.mkdtemp( dir=options.output_dir , prefix='tmp-TOOL-' )

    if options.input_dir_file:
        infile = open(options.input_dir_file, 'r')
        inputDirLine = infile.readline()
        inputDirectory = inputDirLine.rstrip('\n\r')
    else:
        inputDirectory = options.input_dir

    # get the input BAMs into a list
    inputFiles = glob.glob("%s/*.bam" % inputDirectory )
    sampleNames = []
    groupNames = []
    baselineFiles = []

    # get the sample name and bam groups if necessary 
    # for the bam file and store in list in the same order as the Input file list
    for inputF in inputFiles:
        samfile = pysam.AlignmentFile(inputF, 'rb')
        sn = samfile.header['RG'][0]['SM']
   
        # get group names if needed and generate the input baseline files
        base_param = ""
        if options.group_by_keyword is not None and options.rg_field is not None:
            if options.group_by_keyword == "all":
                value = "all_bam"
            else:
                value = samfile.header['RG'][0][options.rg_field]
                if options.field_separator is not None:
                    value = value.split(options.field_separator)[int(options.field_order)-1]

            fh_base = open("%s/inputs/%s.txt" % (baseline_dir, value), 'a')
            fh_base.write("%s\n" % inputF)
            fh_base.close()
            base_param = "-c %s/outputs/%s.baseline.txt" % (baseline_dir, value)
            if value not in groupNames:
                groupNames.append(value)
            if "%s/inputs/%s.txt" % (baseline_dir, value) not in baselineFiles:
                baselineFiles.append("%s/inputs/%s.txt" % (baseline_dir, value) )
        sampleNames.append("%s %s" % (sn, base_param))

    # create the baseline files if needed
    baselineInputs = []
    if options.group_by_keyword is not None:
        # Make sure there are more than 1 BAM file for each Baseline group.
        # if there is only 1 input BAM, then make a link with a new name and add to file
        for inFile in baselineFiles:
            fileCount = lineCount(inFile)
            if fileCount < 2:
                create_additional_bam_copy(inFile, "%s/inputs" % baseline_dir)
                fileCount = 2
            ## create a new baseline file with only x% of inputs depending on the user spec.
            file_qty_to_use = math.ceil(options.percent * fileCount)
            if file_qty_to_use < 2:
                file_qty_to_use = math.ceil(1 * fileCount)
            with open(inFile, "rb") as source:
                lines = [line.rstrip("\n") for line in source]
            print "FILE QTY TO ISE: %s" % file_qty_to_use
            random_choice = random.sample(lines, int(file_qty_to_use))
            newBaseline_input = "%s.selected.txt" % inFile
            baselineInputs.append(newBaseline_input)
            with open(newBaseline_input, "wb") as sink:
                sink.write("\n".join(random_choice))

        baseline_cmd = "export R_LIBS_SITE=%s/site-library:\$R_LIBS_SITE;export PATH=%s/R-3.2.2/bin:%s:%s:%s:\$PATH; baseline.py --file-list INPUTFILE --output %s/outputs/SAMPLENAME -t %s --name baseline --trim 0.2; mv %s/outputs/SAMPLENAME/baseline.pooled2_TRIM0.2.txt %s/outputs/SAMPLENAME.baseline.txt" % (R_bin, R_bin, bedtools_bin, samtools_bin, contra_bin, baseline_dir, options.bed_file, baseline_dir, baseline_dir)

        #if no stderr file is specified, we'll use our own
        stderr = tempfile.NamedTemporaryFile( prefix="TOOL-stderr-", dir=tmp_dir )
        stdout = tempfile.NamedTemporaryFile( prefix="TOOL-stdout-", dir=tmp_dir )

        # prepare command line
        swift_params = list()
        swift_params.append('-inputfiles=\"%s\"' %  ",".join(baselineInputs))
        swift_params.append('-samplenames=\"%s\"' %  ",".join(groupNames))

        ## construct the swift command
        swift_cmd = "%s -sites.file %s -tc.file %s %s " %   (swift_bin, sites_file, tc_file, swift_baseline_file)
        cmd = "%s %s %s" % (swift_cmd, ' '.join(swift_params), '-tool_cmd=\"'+baseline_cmd+'\"')
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


    #### Contra part of the tool command

    if options.pass_through_options:
        ptc= ' '.join( options.pass_through_options )

    tool_cmd = "export R_LIBS_SITE=%s/site-library:\$R_LIBS_SITE;export PATH=%s/R-3.2.2/bin:%s:\$PATH;%s/contra.py -o OUTPUTDIR -s INPUTFILE --sampleName SAMPLENAME -t %s %s --bed" % (R_bin, R_bin, bedtools_bin, contra_bin, options.bed_file, ptc)

    #if no stderr file is specified, we'll use our own
    stderr = tempfile.NamedTemporaryFile( prefix="TOOL-stderr-", dir=tmp_dir )
    stdout = tempfile.NamedTemporaryFile( prefix="TOOL-stdout-", dir=tmp_dir )


    # prepare command line
    swift_params = list()
    swift_params.append('-outputdir=' + output_dir)
    swift_params.append('-inputfiles=\"%s\"' %  ",".join(inputFiles))
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

    # create final VCF file and list output files in the HTML output
    try:
        for sample in sampleNames:
            values = sample.split(" ")
            sn = values[0]
            sample_outputdir = "%s/%s" % (output_dir, sn)

            vcfFile = None
            gainLossFile = None
            for output in glob.glob("%s/table/*" % sample_outputdir):
                #print "OUTPUT: %s" % output
                if output.endswith("vcf"):
                    vcfFile = output
                elif output.endswith("DetailsFILTERED.txt"):
                    gainLossFile = output

            outputVCF = "%s/%s.vcf" % (output_dir, sn)
            modifyContraVCF(vcfFile, gainLossFile, outputVCF, sn)

    except Exception, e:
        sys.stdout.write("problem while generating final VCF " + str(e))

    #### covert contra to metasv format
   
    tool_cmd = "export PATH=%s:\$PATH;%s/contra_to_metasv_run.py --fin INPUTFILE --fout OUTPUTFILE" % (bedtools_bin, contra_bin)

    #if no stderr file is specified, we'll use our own
    stderr = tempfile.NamedTemporaryFile( prefix="TOOL-contra_to_metasv-stderr-", dir=tmp_dir )
    stdout = tempfile.NamedTemporaryFile( prefix="TOOL-contra_to_metasv-stdout-", dir=tmp_dir )


    # prepare command line
    swift_params = list()
    swift_params.append('-outputdir=' + output_dir)
    swift_params.append('-inputfiles=\"%s\"' %  ",".join(inputFiles))
    swift_params.append('-samplenames=\"%s\"' %  ",".join(sampleNames))



    ## construct the swift command
    swift_cmd = "%s -sites.file %s -tc.file %s %s " %   (swift_bin, sites_file, tc_file, swift_contra_to_metasv_file)
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

    createOutputHTML(options.outputF, sampleNames)
    #try:
    #    if os.path.exists(tmp_dir):
    #        shutil.rmtree(tmp_dir)
    #    if os.path.exists(output_dir):
    #        shutil.rmtree(output_dir)
    #except:
    #    pass

    # run swift log summary script
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
