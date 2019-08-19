#!/usr/bin/python

import optparse, os, shutil, sys, tempfile, glob, shlex, vcf, pysam
from subprocess import *
import subprocess

CHUNK_SIZE = 2**20 #1mb

def get_sampleName(fileN):
    fh = open(fileN, "r")
    for line in fh:
        if "#CHROM" in line and "FORMAT\t" in line:
            values = line.rstrip("\n").split("\t")
            return values[-1]
    return os.path.basename(fileN).split(".")[0]

def createOutputHTML (outputF, sampleNames):
    ofh = open(outputF, "w")
    print outputF
    ofh.write( '<html>\n<head>\n<title>Galaxy - LumpyExpress Output</title>\n</head>\n<body>\n<p/>\n<ul>\n' )
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
#    outfiles = glob.glob("%s/*.vcf" % outputDir)

 #   for outF in outfiles:
 #       if os.path.exists(outF):
 #           sn = os.path.basename(outF)
 #           ofh.write('<li><a href="%s">%s</a></li>\n' % ( outF, outF ) )
 #   ofh.write( '</ul>\n</body>\n</html>\n' )
 #   ofh.close()

def __main__():
    parser = optparse.OptionParser()
    parser.add_option( '', '--input_dir', dest='input_dir', action='store', type="string", help='Input directory path of BAM files' )
    parser.add_option( '', '--input_dir_file', dest='input_dir_file', action='store', type="string", help='Input directory File containing path of BAM files' )
    parser.add_option( '', '--out-dir', dest='output_dir', action='store', type="string", default=None, help='If specified, the output directory for extra files.' )
    parser.add_option( '', '--log', dest='swift_log', action='store', type="string", default=None, help='swift summary output.' )
    parser.add_option( '', '--stdout', dest='stdout', action='store', type="string", default=None, help='If specified, the output of stdout will be written to this file.' )
    parser.add_option( '', '--output', dest='outputF', action='store', type="string", default=None, help='mpileup output.' )
    parser.add_option( '', '--bed-filter', dest='filter_bed', action='store', type="string", default=None, help='bed file' )
    parser.add_option( '-p', '--pass_through', dest='pass_through_options', action='append', type="string", help='These options are passed through directly to contra, without any modification.' )
    (options, args) = parser.parse_args()

    swift_bin = 'swift'
    sites_file = '/opt/galaxy/tools/lumpy/lumpy_sites.xml'
    tc_file = '/opt/galaxy/tools/swift/tc.data'
    swift_file_X = '/opt/galaxy/tools/lumpy/lumpyexpress_with_Xfilter.swift'
    swift_file_bed = '/opt/galaxy/tools/lumpy/lumpyexpress_with_bedfilter.swift'
    samtools_bin = "/mnt/galaxyTools/tools/samtools/1.2/bin"
    lumpy_bin = "/mnt/galaxyTools/tools/lumpy/v0.2.11/bin"
    sambamba_bin = "/mnt/galaxyTools/tools/sambamba/v0.5.4"
    samblaster_bin = "/mnt/galaxyTools/tools/samblaster/0.1.22"
    pymodules_bin = "/mnt/galaxyTools/tools/pymodules/python2.7/bin"
    pymodules_lib = "/mnt/galaxyTools/tools/pymodules/python2.7/lib/python"
    svtyper_bin = "/mnt/galaxyTools/tools/svtyper/v0.0.2"
    vcftools_bin = "/mnt/galaxyTools/tools/vcftools/vcftools_0.1.14/bin"
    vcftools_perllib = "/mnt/galaxyTools/tools/vcftools/vcftools_0.1.14/share/perl/5.18.2"

    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)
    output_dir = "%s/output" % options.output_dir
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    tmp_dir = tempfile.mkdtemp( dir=options.output_dir , prefix='tmp-TOOL-' )

    if options.input_dir_file:
        infile = open(options.input_dir_file, 'r')
        inputDirLine = infile.readline()
        inputDirectory = inputDirLine.rstrip('\n\r')
    else:
        inputDirectory = options.input_dir

    # get the input BAMs into a list
    inputFiles = glob.glob("%s/*.bam" % inputDirectory )

    # passthrough tool paramaters
    pass_through = ""
    if options.pass_through_options:
        pass_through = ' '.join( options.pass_through_options )

    #### Contra part of the tool command
    if "-x " not in pass_through and options.filter_bed:
        #tool_cmd = "export PERL5LIB=\$PERL5LIB:%s;export PYTHONPATH=%s:\$PYTHONPATH; export PATH=%s:%s:%s:%s:%s:%s:%s/\$PATH; samtools view -b -L %s PREINPUT > POSTINPUT; samtools index POSTINPUT2; lumpyexpress -o OUTPUT -B INPUTFILE -T TEMPDIR %s; mv TMP_BAI_FILE BAI_FILE; svtyper -M -B OUTPUT_BAM -S SPLITTER_BAM -i OUTPUT_VCF > OUTPUT_GT_VCF; (grep ^\"#\" INPUT_GT_VCF1 ; grep -v ^\"#\" INPUT_GT_VCF2 | LC_ALL=C sort -k1,1 -k2,2n -V) > OUTPUT_FINAL_VCF" % (vcftools_perllib, pymodules_lib, vcftools_bin, sambamba_bin, svtyper_bin, samtools_bin, lumpy_bin, samblaster_bin, pymodules_bin, options.filter_bed, pass_through )
        tool_cmd = "export PERL5LIB=\$PERL5LIB:%s;export PYTHONPATH=%s:\$PYTHONPATH; export PATH=%s:%s:%s:%s:%s:%s:%s/\$PATH; mkdir TEMPDIR; samtools view -b -F 1294 INPUTFILE_RAW > DISCORDANT_BAM_UNSORTED; samtools sort DISCORDANT_BAM_UNSORTED PREFIX_DISCORDANT; samtools view -h INPUTFILE_RAW | /mnt/galaxyTools/tools/lumpy/v0.2.13/scripts/extractSplitReads_BwaMem -i stdin| samtools view -Sb - > SPLITTER_BAM_UNSORTED; samtools sort SPLITTER_BAM_UNSORTED PREFIX_SPLITTER; samtools index SPLITTER_BAM_SORTED; lumpyexpress -o OUTPUT -B INPUTFILE_RAW -S SPLITTER_BAM_SORTED -D DISCORDANT_BAM_SORTED -T TEMPDIR %s; mv TMP_BAI_FILE BAI_FILE; svtyper -M -B OUTPUT_BAM -S SPLITTER_BAM_SORTED -i OUTPUT_VCF > OUTPUT_GT_VCF; (grep ^\"#\" INPUT_GT_VCF1 ; grep -v ^\"#\" INPUT_GT_VCF2 | LC_ALL=C sort -k1,1 -k2,2n -V) > OUTPUT_FINAL_VCF;  vcftools --vcf OUTPUT_FINAL_VCF2 --out OUTPUT_FINAL_VCF3 --bed %s --recode --recode-INFO-all --non-ref-ac 1 " % (vcftools_perllib, pymodules_lib, vcftools_bin, sambamba_bin, svtyper_bin, samtools_bin, lumpy_bin, samblaster_bin, pymodules_bin, pass_through, options.filter_bed )

        swift_file = swift_file_bed
    else:
        tool_cmd = "export PERL5LIB=\$PERL5LIB:%s;export PYTHONPATH=%s:\$PYTHONPATH; export PATH=%s:%s:%s:%s:%s:%s:%s/\$PATH; mkdir TEMPDIR; samtools view -b -F 1294 INPUTFILE_RAW > DISCORDANT_BAM_UNSORTED; samtools sort DISCORDANT_BAM_UNSORTED PREFIX_DISCORDANT; samtools view -h INPUTFILE_RAW | /mnt/galaxyTools/tools/lumpy/v0.2.13/scripts/extractSplitReads_BwaMem -i stdin| samtools view -Sb - > SPLITTER_BAM_UNSORTED; samtools sort SPLITTER_BAM_UNSORTED PREFIX_SPLITTER; samtools index SPLITTER_BAM_SORTED; lumpyexpress -o OUTPUT -B INPUTFILE_RAW -S SPLITTER_BAM_SORTED -D DISCORDANT_BAM_SORTED -T TEMPDIR %s; mv TMP_BAI_FILE BAI_FILE; svtyper -M -B OUTPUT_BAM -S SPLITTER_BAM_SORTED -i OUTPUT_VCF > OUTPUT_GT_VCF; (grep ^\"#\" INPUT_GT_VCF1 ; grep -v ^\"#\" INPUT_GT_VCF2 | LC_ALL=C sort -k1,1 -k2,2n -V) > OUTPUT_FINAL_VCF" % (vcftools_perllib, pymodules_lib, vcftools_bin, sambamba_bin, svtyper_bin, samtools_bin, lumpy_bin, samblaster_bin, pymodules_bin, pass_through )
        swift_file = swift_file_X

    #if no stderr file is specified, we'll use our own
    stderr = tempfile.NamedTemporaryFile( prefix="TOOL-stderr-", dir=tmp_dir )
    stdout = tempfile.NamedTemporaryFile( prefix="TOOL-stdout-", dir=tmp_dir )

    # prepare command line
    swift_params = list()
    swift_params.append('-outputdir=' + output_dir)
    swift_params.append('-inputfiles=\"%s\"' %  ",".join(inputFiles))

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
    sampleNames=[]
    ## Fix the header and rename the output file by taking the sample name
    for fileF in glob.glob("%s/*.vcf" % options.output_dir):
        sampleName = get_sampleName(fileF)
        sampleNames.append("%s" % sampleName)
        new_fh = open("%s/%s.vcf" % (options.output_dir, sampleName), "w")
        fh = open(fileF, "r")
        for line in fh:
            if "##reference=" in line:
                #new_fh.write("##reference=hg19\n")
                pass
            else:
                new_fh.write(line)
        os.remove(fileF)

    createOutputHTML(options.outputF, sampleNames)

#    try:
#        if os.path.exists(output_dir):
#            shutil.rmtree(output_dir)
#    except:
#        pass

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
