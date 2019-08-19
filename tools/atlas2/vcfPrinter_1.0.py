#!/usr/bin/env python
import sys, os, commands, string, subprocess, shutil, optparse, tempfile

""" version 1.0 """

def stop_err( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit()

def __main__():
    galaxyhome=os.environ.get('GALAXY_HOME')

    '''
    n = len(sys.argv)
    for i in range(0,n):
    print sys.argv[i]
    '''
    
    #Parse Command Line     
    parser = optparse.OptionParser()
    parser.add_option( '--inputVCF', dest='inputVCF', action="append", type="string", help='Input VCF file' )
    parser.add_option( '--inputBAM', dest='inputBAM', action="append", type="string", help='Input BAM file' )
    parser.add_option( '--config', dest='config', type="string", help='Configuration file containing path of input VCF and BAM files and Fasta reference file' )
    parser.add_option( '--bamdir', dest='bamdir', type="string", help='Input directory where all BAM file are located' )
    parser.add_option( '--vcfdir', dest='vcfdir', type="string", help='Input directory where all VCF file are located' )
    parser.add_option( '--outputVCF', dest='outputVCF', type="string", help='OutputVCF files' )
    parser.add_option( '--reference', dest='reference', type="string", help='Reference file' )
    parser.add_option( '--refconfig', dest='referenceConfig', type="string", help='Reference Config file' )
    parser.add_option( '--ownfile', dest='ownfile', type="string", help='Reference fasta file from history' )
    
    parser.add_option( '-n', dest='n_param', action="store_true", help='skip pileup' )
    parser.add_option( '-p', dest='p_param', action="store_true", help='Reference fasta file from history' )
    
    (options, args) = parser.parse_args()
    
    
    if options.inputVCF:
        #vcfPrinter.rb requires that the input VCF file and BAM file have the same file name and .vcf/.bam postfix.
        # Copy input VCF file to input_FOLDER/file.vcf
        new_VCF = os.path.join( output_FOLDER, "file.vcf" )
        shutil.copyfile( input_VCF, new_VCF)
        # Copy input BAM file to input_FOLDER/file.bam
        new_BAM = os.path.join( output_FOLDER, "file.bam" )
        shutil.copyfile( input_BAM, new_BAM)

        #get additional input VCF/BAM files
        additional_VCF = ""
        additional_BAM = ""

        i=7
        while (i<len(sys.argv)):
            add_VCF = os.path.join( output_FOLDER, "file%s.vcf" % i)
            shutil.copyfile( sys.argv[i], add_VCF)
            additional_VCF +=  " " + add_VCF
        
            add_BAM = os.path.join( output_FOLDER, "file%s.bam" % i)
            shutil.copyfile( sys.argv[i+1], add_BAM)
            additional_BAM +=  " " + add_BAM
        
            i=i+2

            parameter = " -i " + new_VCF + additional_VCF + " -b " + new_BAM + additional_BAM + input_FASTA + n + p + " -o " + outputVCF
            #print parameter

            ##uncomment when running on VM
            os.chdir(galaxyhome + "/tools/atlas2/")     #uncomment when running on cluster
            #os.chdir("/media/Work/galaxy-proteonics/tools/atlas2/")    #for local running

            command="ruby ./vcfPrinter/vcfPrinter.rb " + parameter
            print command
            proc = subprocess.Popen( args=command, shell=True, stderr=subprocess.PIPE )
            returncode = proc.wait()
            sys.exit()
        
    #### config option
    elif options.config:
        # open file and extract variables
        fh = open(options.config)
        for line in fh:
            if "BAM\t" in line:
                bamdir = line.split("\t")[1]
            elif "VCF\t" in line:
                vcfdir = line.split("\t")[1] + "/*.vcf"
            elif "REF\t" in line:
                ref_path = line.split("\t")[1] + "/*.bam"
                
    elif options.bamdir and options.vcfdir:
        bamdir = options.bamdir + "/*.bam"
        vcfdir = options.vcfdir + "/*.vcf"
    
    if options.n_param:
        n = "-n"
    else:
        n = ""

    if options.p_param:
        p = "-p"
    else:
        p = ""

    if options.reference != "None":
        myref = options.reference
    elif options.ownfile != "None":
        myref = options.ownfile
    else:
        myref = ref_path

    parameters = "-i \"%s\" -o %s -b \"%s\" -r %s %s %s" % (vcfdir, options.outputVCF, bamdir, myref, n, p)
    #cmd = "ruby ./vcfPrinter/vcfPrinter.rb " + parameters
    cmd = "ruby /opt/galaxy/tools/atlas2//vcfPrinter/vcfPrinter.rb " + parameters

    # run:
    try:
        tmp_err = tempfile.NamedTemporaryFile().name
        tmp_stderr = open( tmp_err, 'wb' )
        print cmd
        proc = subprocess.Popen( args=cmd, shell=True, stderr=tmp_stderr, cwd=".")
        returncode = proc.wait()
        tmp_stderr.close()

        # get stderr, allowing for case where it's very large
        tmp_stderr = open( tmp_err, 'rb')
        stderr = ''
        buffsize = 1048576
        try:
            while True:
                stderr += tmp_stderr.read( buffsize )
                if not stderr or len( stderr ) % buffsize != 0:
                    break
        except OverflowError:
            pass
        tmp_stderr.close()

        if returncode != 0:
            raise Exception, stderr

    except Exception, e:
        stop_err( 'Error in vcfPrinter:\n' + str( e ) ) 


if __name__=="__main__": __main__()

