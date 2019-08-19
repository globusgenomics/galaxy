#!/usr/bin/env python
import sys, os, commands, string, subprocess, shutil, optparse, tempfile

""" version 1.1 """

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
    parser.add_option( '--inputPileup', dest='inputPileup', action="append", type="string", help='Input pileup/mpileup file' )
    parser.add_option( '--config', dest='config', type="string", help='Configuration file containing path of input VCF and pileup/mpileup files' )
    parser.add_option( '--vcfdir', dest='vcfdir', type="string", help='Input directory where all VCF file are located' )
    parser.add_option( '--pileupdir', dest='pileupdir', type="string", help='Input directory where all pileup files are located' )
    parser.add_option( '--outputVCF', dest='outputVCF', type="string", help='OutputVCF files' )
    parser.add_option( '-s', dest='region', help='region' )
    parser.add_option( '-p', dest='p_param', action="store_true", help='pass variants' )
    (options, args) = parser.parse_args()

    vcfdir_replace = None
    output_FOLDER = options.outputVCF.split('.')[0]+'_files'
    if not os.path.exists(output_FOLDER):
        os.makedirs(output_FOLDER)
    else:
        pass

    parameters = ""
    
    if options.inputVCF:
        if range(len(options.inputVCF)) > 0:
            #vcfPrinter.rb requires that the input VCF file and pileup files have the same file name and .vcf/.bam postfix.
            # Copy input VCF file to input_FOLDER/file.vcf
            vcfs = ""

            for i in range(len(options.inputVCF)):
                add_VCF = os.path.join( output_FOLDER, "file%s.vcf" % i)
                shutil.copyfile( options.inputVCF[i], add_VCF)
                vcfs +=  " " + add_VCF

            parameters += " -i " + vcfs
        
            if options.inputPileup:
                if range(len(options.inputPileup)) > 0:
                    #vcfPrinter.rb requires that the input VCF file and BAM file have the same file name and .vcf/.bam postfix.
                    # Copy input pileup files to input_FOLDER/file.pileup if available
                    pileups = ""

                    for i in range(len(options.inputPileup)):
                        add_pileup = os.path.join( output_FOLDER, "file%s.pileup" % i)
                        shutil.copyfile( options.inputPileup[i], add_pileup)
                        pileups +=  " " + add_pileup

                    parameters += " -l " + pileups
        vcfdir_replace = add_VCF
    #### config option
    elif options.config and options.vcfdir:
        # open file and extract variables
        fh = open(options.config)
        vcfdir_replace = options.vcfdir
        for line in fh:
            if "VCF\t" in line:
                vcfdir = options.vcfdir + "/*.vcf"
            elif "INDEL\t" in line:
                vcfdir = options.vcfdir + "/*.indel"

        parameters += " -i \"" + vcfdir + "\""

        if options.pileupdir:
            pileupdir = options.pileupdir + "/*.pileup"
            parameters += " -l \"" + pileupdir + "\""

    elif options.vcfdir:
        vcfdir = options.vcfdir + "/*.vcf"
        vcfdir_replace = options.vcfdir
        parameters += " -i \"" + vcfdir + "\""

        if options.pileupdir:
            pileupdir = options.pileupdir + "/*.pileup"
            parameters += " -l \"" + pileupdir + "\""

    if options.region:
        parameters += " -s " + options.region
        
    if options.p_param:
        parameters += " -p"

    ## set the parameters
    parameters += " -o " + options.outputVCF

    #cmd = "ruby ./vcfPrinter/vcfPrinter.rb " + parameters
    tmpoutput = "%s.tmp" % options.outputVCF

    cmd = "ruby /opt/galaxy/tools/atlas2/vcfPrinter/vcfPrinter.rb %s; ruby /opt/galaxy/tools/atlas2/atlas2-code-195-trunk/utils/add_allele_freq.rb %s -clean  > %s" % (parameters,  options.outputVCF, tmpoutput)
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

    shutil.move(tmpoutput, options.outputVCF)

    ## Rename header with the actual filename which should have (i.e. N100.indel should be N100)
    # get all the columns in the output VCF
    f = open(options.outputVCF, "r")
    header = None
    for line in f:
        if line.startswith("#CHROM\t"):
            header = line.rstrip("\n")
            break
    f.close()

    # Starting in the 10th column are the sample columns we should change names for
    cols = header.split("\t")
    new_header = cols[0:9]
    samples = cols[9:]
    for sample in samples:
        filename = "%s/%s" % (vcfdir_replace, sample)
        if os.path.isfile(filename):
            # get the column name for the sample
            f = open(filename, "r")
            head = None
            for line in f:
                if line.startswith("#CHROM\t"):
                    head = line.rstrip("\n")
                    break
            f.close()
            samplename = head.split("\t")[-1]
            # replace sample with samplename
            new_header.append(samplename)
    # replace header in the file
    f = open(options.outputVCF, "r")
    tmpoutput = "%s.tmp" % options.outputVCF
    fo = open(tmpoutput, "w")
    header = None
    for line in f:
        if line.startswith("#CHROM\t"):
            line = "\t".join(new_header) + "\n" 
        fo.write(line)
    f.close()
    fo.close()
    shutil.move(tmpoutput, options.outputVCF)


    # fix the lines where the reference allele shows up in the alternate allele column (i.e. chr2	47680085	.	GCA	G,GCA)
    # and fix the sample columns appropriately (i.e. 2/2 to 0/0 or for non-repetetive alleles modify the allele number) 
    f = open(options.outputVCF, "r")
    tmpoutput = "%s.tmp" % options.outputVCF
    fo = open(tmpoutput, "w")
    for line in f:
        line = line.rstrip()
        if not line.startswith("#"):
            cols = line.split("\t")
            if "," in cols[4]:
                values_alt_col = cols[4].split(",")
                if cols[3] in values_alt_col:
                    number = 0
                    genotype_number = {cols[3] : number}
                    new_alt = []
                    for alt in values_alt_col:
                        if alt != cols[3]:
                            number += 1
                            genotype_number[alt] = number
                            new_alt.append(alt)
                        else:
                            genotype_number[alt] = "0"

                    cols[4] = ",".join(new_alt)
                    for i in xrange(9, len(cols)):
                        geno_values = cols[i].split(":")
                        if "/" in geno_values[0]:
                            (alt_geno, ref_geno) = geno_values[0].split("/")
                            if geno_values[-1] in genotype_number and alt_geno != ".":
                                if geno_values[-1] == cols[3]:
                                    ref_geno = 0
                                if ref_geno > genotype_number[geno_values[-1]] and ref_geno == alt_geno:
                                    ref_geno = genotype_number[geno_values[-1]]
                                geno_values[0] = "%s/%s" % (genotype_number[geno_values[-1]], ref_geno)
                        cols[i] = ":".join(geno_values)
                    line = "\t".join(cols)
        fo.write(line + "\n")    
    f.close()
    fo.close()
    #shutil.move(tmpoutput, options.outputVCF)

    # re-run add_allele_freq.rb on new VCF
    cmd = "ruby /opt/galaxy/tools/atlas2/atlas2-code-195-trunk/utils/add_allele_freq.rb %s -clean  > %s" % (tmpoutput, options.outputVCF)  
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
        
    os.remove(tmpoutput)

if __name__=="__main__": __main__()
