#!/usr/bin/env python

"""
A wrapper script for running swift on globus-galaxy
"""

import optparse, os, shutil, sys, tempfile, glob, json 
from subprocess import *
import subprocess

CHUNK_SIZE = 2**20 #1mb

def createOutputHTML (outputF, sampleNames, outputDir):
    ofh = open(outputF, "w")
    ofh.write( '<html>\n<head>\n<title>Galaxy - vcffilter Output</title>\n</head>\n<body>\n<p/>\n<ul>\n' )

    for sample in sampleNames:
        outVCF = "%s/%s.vcf" % (outputDir, sample)
        ofh.write('<li><a href="%s">%s</a></li>\n' % ( outVCF, sample ) )
    ofh.write( '</ul>\n</body>\n</html>\n' )
    ofh.close()

def __main__():
    parser = optparse.OptionParser()
    parser.add_option('--vcf-dir', dest='vcf_dir', type="string", help='Input directory where all VCF file are located' )
    parser.add_option('--config', dest='config_file', type="string", help='Config file with location of VCF files' )
    parser.add_option('--output-dir', dest='output_dir', type="string", help='Output directory to store recoded VCF files' )
    parser.add_option( '', '--output', dest='outputF', type="string")
    parser.add_option('--bed', dest='bedfile', type="string", help='Bed file' )
    parser.add_option( '-p', '--pass_through', dest='pass_through_options', action='append', type="string", help='These options are passed through directly to contra, without any modification.' )
    (opts, args) = parser.parse_args()

    swift_bin = '/mnt/galaxyTools/tools/swift/swift-0.94.1/bin/swift'
    sites_file = '/opt/galaxy/tools/swift/sites.xml'
    tc_file = '/opt/galaxy/tools/swift/tc.data'
    swift_file = '/opt/galaxy/tools/swift/filter_sv_vcf.swift'
    vcftools_bin = '/mnt/galaxyTools/tools/vcftools/vcftools_0.1.14/bin'
    vcflib_bin = '/mnt/galaxyTools/tools/vcflib/10.27.2016/bin'

    if not os.path.exists(opts.output_dir):
        os.mkdir(opts.output_dir)
    output_dir = "%s/output" % opts.output_dir
    inputDirectory = "%s/vcfs" % opts.output_dir
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
        os.mkdir(inputDirectory)

    tmp_dir = tempfile.mkdtemp( dir=opts.output_dir, prefix='tmp-TOOL-' )

## parse mandatory commands

    if opts.vcf_dir:
        vcf_directory = opts.vcf_dir
    elif opts.config_file:
        vcf_directory = opts.config_file.split('.')[0]+'_files'
        #print vcf_directory
    else:
       sys.exit()

    sampleNames=[]
    configFiles=[]

    for inputF in glob.glob("%s/*.vcf" % vcf_directory):
        os.symlink(inputF, "%s/%s" % (inputDirectory, os.path.basename(inputF)))
        configFiles.append("%s/%s" % (inputDirectory, os.path.basename(inputF)))
        sampleNames.append(os.path.basename(inputF).split('.')[0])
        
    ptc = ""
    if opts.pass_through_options:
        ptc = ' '.join( opts.pass_through_options )
        #print "ptc1: %s" % ptc
        if '__gt__' in ptc:
            ptc = ptc.replace('__gt__', '>')
        elif '__lt__' in ptc:
            ptc = ptc.replace('__lt__', '<')
        else:
            print "ptc: %s" % ptc
        ptc = json.dumps(ptc)
        ptc = ptc[1:-1]
    #print ptc
    if opts.bedfile is not None:
        if ptc != "":
            tool_cmd = "export PATH=%s:%s:\$PATH; %s/vcftools --vcf INPUTFILE --out TMPFILE --bed %s --recode --recode-INFO-all; %s/vcffilter %s TMPFILE > OUTPUTFILE" % (vcftools_bin, vcflib_bin, vcftools_bin, opts.bedfile, vcflib_bin, ptc)
        else:
            tool_cmd = "export PATH=%s:%s:\$PATH; %s/vcftools --vcf INPUTFILE --out TMPFILE --bed %s --recode --recode-INFO-all; cat TMPFILE > OUTPUTFILE" % (vcftools_bin, vcflib_bin, vcftools_bin, opts.bedfile)
    else:
        tool_cmd = "export PATH=%s:\$PATH;%s/vcffilter %s INPUTFILE > OUTPUTFILE" % (vcflib_bin, vcflib_bin, ptc)
    
    #print "tool cmd: %s " % tool_cmd

    #if no stderr file is specified, we'll use our own
    stderr = tempfile.NamedTemporaryFile( prefix="TOOL-stderr-", dir=tmp_dir )
    stdout = tempfile.NamedTemporaryFile( prefix="TOOL-stdout-", dir=tmp_dir )

    swift_params = list()
    swift_params.append('-output_dir=' +  output_dir)
    swift_params.append('-inputfiles=\"%s\"' % ",".join(configFiles))
    swift_params.append('-samplenames=\"%s\"' %  ",".join(sampleNames))


## construct the swift command
    swift_cmd = "%s -sites.file %s -tc.file %s %s " %   (swift_bin, sites_file, tc_file, swift_file)
    #cmd = "%s %s %s 2>&1" % (swift_cmd, ' '.join(swift_params), '-tool_cmd=\"'+tool_cmd+'\"')
    cmd = "%s %s %s " % (swift_cmd, ' '.join(swift_params), '-tool_cmd=\"'+tool_cmd+'\"')
    print "cmd: %s" % cmd
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

    try:
        tmp_dir=opts.outputF.split('.')[0]+'_files'
        createOutputHTML(opts.outputF, sampleNames, tmp_dir)
    except Exception, e:
        sys.stdout.write("problem while generating final VCF " + str(e))

if __name__=="__main__":
    __main__()
