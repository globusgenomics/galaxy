#!/usr/bin/python

import optparse, os, shutil, sys, tempfile, glob, shlex, vcf, pysam
from subprocess import *
import subprocess

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

def __main__():
    parser = optparse.OptionParser()
    parser.add_option( '-p', '--pass_through', dest='pass_through_options', action='append', type="string", help='These options are passed through directly to contra, without any modification.' )
    parser.add_option( '-n', '--sampleName', dest='sampleName', help='sampleName' )
    parser.add_option( '-o', '--contra-vcf', dest='outputVCF', help='outpuVCF' )
    (options, args) = parser.parse_args()

    tempDir = tempfile.mkdtemp();
    output_dir = "%s/contra-output" % tempDir

    # get the sample name for the bam file and store in dictionary
    #samfile = pysam.AlignmentFile(linkname,'rb')
    #sn = samfile.header['RG'][0]['SM']
    #filebase = '.'.join(os.path.basename(linkname).split('.')[:-1]) 
    #sample_names[filebase] = sn

    cmd = "contra.py -o %s --sampleName %s " % (output_dir, options.sampleName)

    if options.pass_through_options:
        cmd += ' '.join( options.pass_through_options )

    try:
        print cmd
        run_cmd( cmd, "runing contra " )
    except Exception, e:
        sys.stdout.write("problem while runing contra " + str(e))
        sys.exit()

    try:
        vcfFile = None
        gainLossFile = None
        for output in glob.glob("%s/table/*" % output_dir):
            print "OUTPUT: %s" % output
            if output.endswith("vcf"):
                vcfFile = output
            elif output.endswith("DetailsFILTERED.txt"):
                gainLossFile = output

        modifyContraVCF(vcfFile, gainLossFile, options.outputVCF, options.sampleName)

    except Exception, e:
        sys.stdout.write("problem while generating final VCF " + str(e))

#    finally:
#        try:
#            if os.path.exists(tempDir):
#                shutil.rmtree(tempDir)
#        except:
#            pass

if __name__=="__main__":
	__main__()
