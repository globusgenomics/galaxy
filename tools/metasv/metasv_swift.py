#!/usr/bin/python

import optparse, os, shutil, sys, tempfile, glob, shlex, vcf, pysam, gzip
from subprocess import *
import subprocess

CHUNK_SIZE = 2**20 #1mb

def create_final_output_files(outputDir, finalOutputDir):
    for x in os.listdir(outputDir):
        base=x
        tmpfile="%s/%s/variants.vcf.gz" % (outputDir, base)
        dest_file="%s/%s.vcf" % (finalOutputDir, base)
        with gzip.open(tmpfile, 'rb') as infile:
            with open(dest_file, 'w') as outfile:
                for line in infile:
                    outfile.write(line)

def createOutputHTML (outputF, sampleNames):
    ofh = open(outputF, "w")
    print outputF
    ofh.write( '<html>\n<head>\n<title>Galaxy - MetaSV Output</title>\n</head>\n<body>\n<p/>\n<ul>\n' )
    outputDir='%s_files' % ''.join(outputF.split('.')[:-1])
    for sample in sampleNames:
        values = sample.split(" ")
        sn = values[0]
        outVCF = "%s/%s.vcf" % (outputDir, sn)
        #if os.path.exists(outVCF):
        ofh.write('<li><a href="%s">%s</a></li>\n' % ( outVCF, sn ) )
    ofh.write( '</ul>\n</body>\n</html>\n' )
    ofh.close()

def open_file_from_option( filename, mode = 'rb' ):
    if filename:
        return open( filename, mode = mode )
    return None

def get_samplename(vcf_path):
    sample_name = os.path.basename(vcf_path).split(".")[0]
    return sample_name


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
    parser.add_option( '', '--input_bam_dir', dest='input_bam_dir', action='store', type="string", help='Input directory path of BAM files' )
    parser.add_option( '', '--input_bam_dir_file', dest='input_bam_dir_file', action='store', type="string", help='Input directory File containing path of BAM files' )
    parser.add_option( '', '--input_bam_files', dest='input_bam_files', action='store', type="string", help='Input File list containing path of BAM files' )
    parser.add_option( '', '--input_dir', dest='input_vcf_dir', action='append', nargs=2, type="string", help='Input directory path of VCF files' )
    parser.add_option( '', '--out-dir', dest='output_dir', action='store', type="string", default=None, help='If specified, the output directory for extra files.' )
    parser.add_option( '', '--log', dest='swift_log', action='store', type="string", default=None, help='swift summary output.' )
    parser.add_option( '', '--output', dest='outputF', action='store', type="string", default=None, help='vcf output.' )
    parser.add_option( '-p', '--pass_through', dest='pass_through_options', action='append', type="string", help='These options are passed through directly to contra, without any modification.' )
    (options, args) = parser.parse_args()

    swift_bin = 'swift'
    sites_file = '/opt/galaxy/tools/metasv/metasv_sites.xml'
    tc_file = '/opt/galaxy/tools/swift/tc.data'
    swift_file = '/opt/galaxy/tools/metasv/metasv.swift'
    metasv_bin = "/mnt/galaxyTools/tools/pymodules/python2.7/bin"
    pythonpath = "/mnt/galaxyTools/tools/pymodules/python2.7/lib/python"
    bedtools_bin="/mnt/galaxyTools/tools/bedtools/2.25.0/bin"
    spades_bin = "/mnt/galaxyTools/tools/spades/3.6.2/bin"
    age_bin = "/mnt/galaxyTools/tools/age/02-22-2016"
    ld_libs = "/mnt/galaxyTools/tools/spades/3.6.2/lib:/mnt/galaxyTools/tools/pymodules/python2.7/include/python/pygsl:/mnt/galaxyTools/tools/pymodules/python2.7/lib/python/pygsl"

    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)
    output_dir = "%s/output" % options.output_dir
    inputBAMdirectory = "%s/bams" % options.output_dir
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
        os.mkdir(inputBAMdirectory)

    tmp_dir = tempfile.mkdtemp( dir=options.output_dir , prefix='tmp-TOOL-' )
    sampleNames =[]
    input_samples = {}
    basenames = {}
    tool_list= len(options.input_vcf_dir)

    for (tool_name, dir_path) in options.input_vcf_dir:
        if tool_name == "--breakdancer_native":
            inputFiles = glob.glob("%s/*.txt" % dir_path ) 
        else:
            inputFiles = glob.glob("%s/*.vcf" % dir_path )
        print ("%s:%d inputfiles" % (tool_name, len(inputFiles)))
        input_samples[tool_name] = []
        for f in sorted(inputFiles):
           input_samples[tool_name].append(f) 
           sample_name = get_samplename(f)
           basenames[sample_name] = 1
           #print "%s: %s\n" % (tool_name,sample_name)
 
    # add output dir.
    for outname in basenames:
        os.mkdir("%s/%s" % (output_dir, outname) )

    # pass through options
    ptc = ""
    if options.pass_through_options:
        ptc = ' '.join( options.pass_through_options )

    # get the sample name and bam groups if necessary 
    # for the bam file and store in list in the same order as the Input file list
    #### Run the configManta command for each input on the head node since it's a low cost job
    if options.input_bam_dir:
        inputFiles = []; inputIndexes=[]; sampleNames = []; inputLinkedFiles = [];
        inputFiles = sorted(glob.glob("%s/*.bam" % options.input_bam_dir))
        inputIndexes = sorted(glob.glob("%s/*.bai" % options.input_bam_dir))
        index_n =0

        for inputF in inputFiles:
            samfile = pysam.AlignmentFile(inputF, 'rb')
            sn = samfile.header['RG'][0]['SM']

            sampleNames.append("%s" % (sn))
        # create the output directory for the sample
            sample_outD = "%s/%s" % (output_dir, sn)
            if not os.path.exists(sample_outD):
                os.mkdir(sample_outD)

        # create symlink for input file
            os.symlink(inputF, "%s/%s.bam" % (inputBAMdirectory, sn))
            os.symlink(inputIndexes[index_n], "%s/%s.bai" % (inputBAMdirectory, sn))
            inputLinkedFiles.append("%s/%s.bam" % (inputBAMdirectory, sn))
            index_n += 1
     
        tool_cmd = "export LD_LIBRARY_PATH=%s:\$LD_LIBRARY_PATH; export PATH=%s:%s:%s:%s:\$PATH;export PYTHONPATH=%s:\$PYTHONPATH; run_metasv.py TOOLS --sample SAMPLE --bam BAM --outdir OUTDIR %s --num_threads 1" % (ld_libs, metasv_bin, spades_bin, age_bin, bedtools_bin, pythonpath, ptc)
    else: 
    # prepare tool command
        tool_cmd = "export LD_LIBRARY_PATH=%s:\$LD_LIBRARY_PATH; export PATH=%s:%s:%s:%s:\$PATH;export PYTHONPATH=%s:\$PYTHONPATH; run_metasv.py TOOLS --sample SAMPLE --outdir OUTDIR %s" % (ld_libs, metasv_bin, spades_bin, age_bin, bedtools_bin, pythonpath, ptc)

    #if no stderr file is specified, we'll use our own
    stderr = tempfile.NamedTemporaryFile( prefix="TOOL-stderr-", dir=tmp_dir )
    stdout = tempfile.NamedTemporaryFile( prefix="TOOL-stdout-", dir=tmp_dir )

    # prepare command line
    swift_params = list()
    swift_params.append('-inputbamdir=' + inputBAMdirectory)
    swift_params.append('-outputdir=' + output_dir)
    swift_params.append('-samplenames=\"%s\"' %  ",".join(sorted(basenames.keys())))
    swift_params.append('-tools=\"%s\"' % ",".join(input_samples.keys()))
    #for tool in input_samples:
    #    swift_params.append('-%s=\"%s\"' %  (tool, ",".join(input_samples[tool])))
    for tool in input_samples:
        #print input_samples[tool][0]
        swift_params.append('-%s=\"%s\"' %  (tool, os.path.dirname(input_samples[tool][0])))
 
    ## construct the swift command
    swift_cmd = "%s -sites.file %s -tc.file %s %s " %   (swift_bin, sites_file, tc_file, swift_file)
    cmd = "%s %s %s" % (swift_cmd, ' '.join(swift_params), '-tool_cmd=\"'+tool_cmd+'\"')
    print tool_cmd
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

    create_final_output_files(output_dir, options.output_dir)

    # create list output files in the HTML output
    try:
        sampleNames=sorted(basenames.keys())
        print sampleNames
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
