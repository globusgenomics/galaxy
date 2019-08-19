#!/usr/bin/python

import optparse, os, shutil, sys, tempfile, glob, shlex, vcf, pysam
from subprocess import *
import subprocess, math, random

CHUNK_SIZE = 2**20 #1mb

def createOutputHTML (outputF, sampleNames):
    ofh = open(outputF, "w")
    print outputF
    ofh.write( '<html>\n<head>\n<title>Galaxy - CNVKIT VCF Output</title>\n</head>\n<body>\n<p/>\n<ul>\n' )
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

def lineCount(inFile):
    fh = open(inFile, "r")
    counter = 0
    for line in fh:
        counter += 1
    return counter

def __main__():
    parser = optparse.OptionParser()
    parser.add_option( '', '--input_dir', dest='input_dir', action='store', type="string", help='Input directory path of BAM files' )
    parser.add_option( '', '--input_dir_file', dest='input_dir_file', action='store', type="string", help='Input directory File containing path of BAM files' )
    parser.add_option( '', '--input_files', dest='input_files', action='store', type="string", help='Input File list containing path of BAM files' )
    parser.add_option( '', '--out-dir', dest='output_dir', action='store', type="string", default=None, help='If specified, the output directory for extra files.' )
    parser.add_option( '', '--log', dest='swift_log', action='store', type="string", default=None, help='swift summary output.' )
    parser.add_option( '', '--output', dest='outputF', action='store', type="string", default=None, help='mpileup output.' )
    parser.add_option( '-p', '--pass_through', dest='pass_through_options', action='append', type="string", help='These options are passed through directly to contra, without any modification.' )
    parser.add_option( '-c', '--config', dest='config_file', action="store", type="string", default=None )
    parser.add_option( '-t', '--target', dest='bed_file', action="store", type="string", default=None )
    parser.add_option( '-r', '--reference_path', dest='reference_path', help="reference file" )
    parser.add_option( '', '--percent-bam-files-for-baseline', type="float", dest='percent', help='contra baseline group: percentage of BAMs to use' )
    parser.add_option( '', '--baseline-input-bam', dest='group_by_keyword', help='contra baseline group: to group or not to group' )
    parser.add_option( '', '--group-field', dest='rg_field', help='contra baseline group: RG field to use for grouping' )
    parser.add_option( '', '--keyword-separator', dest='field_separator', help='contra baseline group: RG field separator' )
    parser.add_option( '', '--keyword-field-order', dest='field_order', help='contra baseline group: RG field order' )
    (options, args) = parser.parse_args()

    swift_bin = 'swift'
    sites_file = '/opt/galaxy/tools/cnvkit/cnvkit_sites.xml'
    tc_file = '/opt/galaxy/tools/swift/tc.data'
    swift_file = '/opt/galaxy/tools/cnvkit/cnvkit.swift'
    cnvkit_bin = "/mnt/galaxyTools/tools/pymodules/python2.7/bin"
    pythonpath = "/mnt/galaxyTools/tools/pymodules/python2.7/lib/python"
    r_path = "/mnt/galaxyTools/tools/R/3.2.2/bin/bin"
    r_libs = "/mnt/galaxyTools/tools/R/3.2.2/site-library"
    r_ld = "/mnt/galaxyTools/tools/R/3.2.2/ld_libs"

    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)
    output_dir = "%s/output" % options.output_dir
    inputDirectory = "%s/bams" % options.output_dir
    baseline_dir = "%s/baseline" % options.output_dir
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
        os.mkdir(inputDirectory)
        os.mkdir(baseline_dir)

    tmp_dir = tempfile.mkdtemp( dir=options.output_dir , prefix='tmp-TOOL-' )

    inputLinkedFiles = []
    basenames = []
    if not options.input_files:
        if options.input_dir_file:
            infile = open(options.input_dir_file, 'r')
            inputDirLine = infile.readline()
            inputRealDirectory = inputDirLine.rstrip('\n\r')
        elif options.input_dir:
            inputRealDirectory = options.input_dir

        inputFiles = glob.glob("%s/*.bam" % inputRealDirectory )
        # create a link of the BAMs inside the inputs directory
        # this is to prevent the Warning of the BAI files being older than BAMs
        #for inputF in inputFiles:
        #    #os.symlink(inputF, "%s/%s" % (inputDirectory, os.path.basename(inputF)))
        #    inputLinkedFiles.append("%s/%s" % (inputDirectory, os.path.basename(inputF)))
    else:
        inputFiles = options.input_files.strip(" ").split(" ")
        #for inputF in inputFiles:
        #    #os.symlink(inputF, "%s/%s" % (inputDirectory, os.path.basename(inputF)))
        #    inputLinkedFiles.append("%s/%s" % (inputDirectory, os.path.basename(inputF)) )

    # get the input BAMs into a list
    sampleNames = []
    groupNames = []
    baselineFiles = []

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

                # get group names if needed and generate the input baseline files
        base_param = ""
        if options.group_by_keyword is not None and options.rg_field is not None:
            if options.group_by_keyword == "all":
                value = "all_bam"
            else:
                value = samfile.header['RG'][0][options.rg_field]
                if options.field_separator is not None:
                    value = value.split(options.field_separator)[int(options.field_order)-1]

            fh_base = open("%s/%s.txt" % (baseline_dir, value), 'a')
            fh_base.write("%s\n" % inputF)
            fh_base.close()
            base_param = "-c %s/%s.baseline.txt" % (baseline_dir, value)
            if value not in groupNames:
                groupNames.append(value)
            if "%s/%s.txt" % (baseline_dir, value) not in baselineFiles:
                baselineFiles.append("%s/%s.txt" % (baseline_dir, value) )

    # pass through options
    ptc = ""
    if options.pass_through_options:
        ptc = ' '.join( options.pass_through_options )

    # create the baseline files if needed
    baselineInputs = []
    if options.group_by_keyword is not None:
        # Make sure there are more than 1 BAM file for each Baseline group.
        # if there is only 1 input BAM, then make a link with a new name and add to file
        for inFile in baselineFiles:
            fileCount = lineCount(inFile)
            if fileCount < 2:
                create_additional_bam_copy(inFile, baseline_dir)
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
            baselineInputs= baselineInputs+random_choice
            #baselineInputs.append(newBaseline_input)
            with open(newBaseline_input, "wb") as sink:
                sink.write("\n".join(random_choice))
        print baselineInputs
        baseline_bam = ' '.join(baselineInputs)
        # prepare baseline command

        baseline_cmd = "export R_LIBS_SITE=%s:\$R_LIBS_SITE; export LD_LIBRARY_PATH=%s:\$LD_LIBRARY_PATH; export PATH=%s:%s:\$PATH;export PYTHONPATH=%s:\$PYTHONPATH; cnvkit.py batch %s -n --output-dir %s --output-reference %s/flat_reference.cnn %s" % (r_libs, r_ld, r_path, cnvkit_bin, pythonpath, baseline_bam, baseline_dir, baseline_dir, ptc )
        print "baseline_cmd: %s " % baseline_cmd
        stderr = tempfile.NamedTemporaryFile( prefix="TOOL-baseline-stderr-", dir=tmp_dir )
        stdout = tempfile.NamedTemporaryFile( prefix="TOOL-baseline-stdout-", dir=tmp_dir ) 

        return_code = None

        if return_code is None or not return_code:
            proc = subprocess.Popen( args=baseline_cmd, stdout=stdout, stderr=stderr, shell=True, cwd=tmp_dir )
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
    # prepare tool command
    #tool_cmd = "export R_LIBS_SITE=%s:\$R_LIBS_SITE; export LD_LIBRARY_PATH=%s:\$LD_LIBRARY_PATH; export PATH=%s:%s:\$PATH;export PYTHONPATH=%s:\$PYTHONPATH; cnvkit.py batch INPUTBAM -n --output-dir OUTPUTDIR --output-reference REFNAME %s; cnvkit.py call INPUTCNS -y -m threshold -t=-1.1,-0.4,0.3,0.7 -o OUTPUTCALLFILE; cnvkit.py export vcf INPUTCALLFILE -o OUTPUTVCF; cp INT_VCF FINAL_VCF" % (r_libs, r_ld, r_path, cnvkit_bin, pythonpath, ptc)
    tool_cmd = "export R_LIBS_SITE=%s:\$R_LIBS_SITE; export LD_LIBRARY_PATH=%s:\$LD_LIBRARY_PATH; export PATH=%s:%s:\$PATH;export PYTHONPATH=%s:\$PYTHONPATH; cnvkit.py batch INPUTBAM --output-dir OUTPUTDIR -r %s/flat_reference.cnn; cnvkit.py call INPUTCNS -y -m threshold -t=-1.1,-0.4,0.3,0.7 -o OUTPUTCALLFILE; cnvkit.py export vcf INPUTCALLFILE -o OUTPUTVCF; cp INT_VCF FINAL_VCF" % (r_libs, r_ld, r_path, cnvkit_bin, pythonpath, baseline_dir)

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
    for vcffile in glob.glob("%s/*.vcf" % output_dir):
        shutil.copy(vcffile, options.output_dir)
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
