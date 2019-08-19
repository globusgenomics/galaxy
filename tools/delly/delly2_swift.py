#!/usr/bin/python

import optparse, os, shutil, sys, tempfile, glob, shlex, vcf, pysam
from subprocess import *
import subprocess

CHUNK_SIZE = 2**20 #1mb

def map_svTypes_to_sitesFiles(types, vcfdir):
    maps = {}
    for sv_type in types:
        files = glob.glob("%s/%s.bcf" % (vcfdir, sv_type))
        maps[sv_type] = files[0]
    return maps

def regenotyping_map_svTypes_to_sitesFiles(vcfdir):
    files=glob.glob("%s/*.bcf" % vcfdir)
    print "VCFDIR: %s" % vcfdir
    print "FILES: %s" % files
    maps={}
    for j in range(0,len(files)):
        filename=os.path.basename(files[j])
        #if len(filename)==7:
        if "_TRA" in filename or "_INV" in filename or "_DEL" in filename or "_INS" in filename or "_DUP" in filename:
            #sv_type=filename[0:3]
            if "_TRA" in filename:
                sv_type = "TRA"
            elif "_INV" in filename:
                sv_type = "INV"
            elif "_DEL" in filename:
                sv_type = "DEL"
            elif "_INS" in filename:
                sv_type = "INS"
            elif "_DUP" in filename:
                sv_type = "DUP"


            sv_file="%s/%s" % (vcfdir, filename)
            maps[sv_type]=sv_file
    return maps

def createOutputHTML (outputF, sampleNames):
    ofh = open(outputF, "w")
    ofh.write( '<html>\n<head>\n<title>Galaxy - Delly2 VCF Output</title>\n</head>\n<body>\n<p/>\n<ul>\n' )
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
    parser.add_option( '-t', dest='types', help='SV analysis type (DEL, DUP, INV, TRA), INS', action='store', nargs=5, type="string" )
    (options, args) = parser.parse_args()

    swift_bin = 'swift'
    sites_file = '/opt/galaxy/tools/delly/delly2_sites.xml'
    tc_file = '/opt/galaxy/tools/swift/tc.data'
    swift_file = '/opt/galaxy/tools/delly/delly2.swift'
    samtools_bin = "/mnt/galaxyTools/tools/samtools/1.2/bin/samtools" # for pysam
    delly2_bin = "/mnt/galaxyTools/tools/delly/v0.7.3"
    delly2 = "/mnt/galaxyTools/tools/delly/v0.7.3/delly2"
    bedtools_bin = "/mnt/galaxyTools/tools/bedtools/2.17.0/bin"
    swift_delly_to_metasv_file= '/opt/galaxy/tools/delly/delly_to_metasv.swift'    
 
    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)
    output_dir = "%s/output" % options.output_dir
    inputDirectory = "%s/bams" % options.output_dir
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
        os.mkdir(inputDirectory)

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

        inputFiles = sorted(glob.glob("%s/*.bam" % inputRealDirectory ))
        inputIndexes = sorted(glob.glob("%s/*.bai" % inputRealDirectory ))
        # create a link of the BAMs inside the inputs directory
        # this is to prevent the Warning of the BAI files being older than BAMs
        #for inputF in inputFiles:
        #    #os.symlink(inputF, "%s/%s" % (inputDirectory, os.path.basename(inputF)))
        #    inputLinkedFiles.append("%s/%s" % (inputDirectory, os.path.basename(inputF)))
    else:
        inputFiles = options.input_files.strip(" ").split(" ")
        inputIndexes = sorted(glob.glob("%s/*.bai" % inputRealDirectory ))
        #for inputF in inputFiles:
        #    #os.symlink(inputF, "%s/%s" % (inputDirectory, os.path.basename(inputF)))
        #    inputLinkedFiles.append("%s/%s" % (inputDirectory, os.path.basename(inputF)) )

    # get the input BAMs into a list
    sampleNames = []

    # pass through options
    ptc = ""
    if options.pass_through_options:
        ptc = ' '.join( options.pass_through_options )

    # get the sample name and bam groups if necessary 
    # for the bam file and store in list in the same order as the Input file list
    #### Run the configManta command for each input on the head node since it's a low cost job
    index_n = 0
    for inputF in inputFiles:
        samfile = pysam.AlignmentFile(inputF, 'rb')
        sn = samfile.header['RG'][0]['SM']
   
        sampleNames.append("%s" % (sn))
        # create the output directory for the sample
        sample_outD = "%s/%s" % (output_dir, sn) 
        os.mkdir(sample_outD)

        # create symlink for input file
        os.symlink(inputF, "%s/%s.bam" % (inputDirectory, sn))
        os.symlink(inputIndexes[index_n], "%s/%s.bai" % (inputDirectory, sn))
        inputLinkedFiles.append("%s/%s.bam" % (inputDirectory, sn))
        index_n += 1

    # prepare tool command
    print "HERE: %s" % options.vcffile
    if options.vcffile == "None":
        env_cmd = "export LD_LIBRARY_PATH=/mnt/galaxyTools/tools/samtools/1.2/lib/:/mnt/galaxyTools/tools/delly/v0.7.3/lib:\$LD_LIBRARY_PATH; export PATH=%s:%s:\$PATH;" % (samtools_bin, delly2_bin)
        tool_cmd = " %s call -o OUTPUTDIR -t %s %s INPUTFILE; %s call -o OUTPUTDIR -t %s %s INPUTFILE; %s call -o OUTPUTDIR -t %s %s INPUTFILE; %s call -o OUTPUTDIR -t %s %s INPUTFILE; %s call -o OUTPUTDIR -t %s %s INPUTFILE" % (delly2, options.types[0], ptc, delly2, options.types[1], ptc, delly2, options.types[2], ptc, delly2, options.types[3], ptc, delly2, options.types[4], ptc)

    else:
        #svlist=regenotyping_map_svTypes_to_sitesFiles(options.vcffile)
        #svlist=glob.glob("%s/*.bcf" % options.vcffile)
        #types={}
        #for j in range(0,len(svlist)):
        #    if len(os.path.basename(svlist[j]))==7:
        #        svName=os.path.basename(svlist[j])[0:3]
        #        fileName="%s/%s" % (options.vcffile, os.path.basename(svlist[j]))
        #        types[svName]=fileName
        env_cmd = "export LD_LIBRARY_PATH=/mnt/galaxyTools/tools/samtools/1.2/lib/:/mnt/galaxyTools/tools/delly/v0.7.3/lib:\$LD_LIBRARY_PATH; export PATH=%s:%s:\$PATH;" % (samtools_bin, delly2_bin)
        svlist=regenotyping_map_svTypes_to_sitesFiles(options.vcffile)
        print "SVLIST: %s" % svlist
        tmp=[]
        svtypes=[] 
        for key in svlist:
            print "KEY: %s" % key
            svtypes.append(key)
            tmp.append(" %s call -o OUTPUTDIR -t %s %s -v %s INPUTFILE;" % (delly2, key, ptc, svlist[key]))
        tool_cmd =''.join(tmp)     
        print "tool_cmd:%s" % tool_cmd
        print svtypes
        print tmp

    #if no stderr file is specified, we'll use our own
    stderr = tempfile.NamedTemporaryFile( prefix="TOOL-stderr-", dir=tmp_dir )
    stdout = tempfile.NamedTemporaryFile( prefix="TOOL-stdout-", dir=tmp_dir )

    # prepare command line
    swift_params = list()
    swift_params.append('-envcmd=\"%s\"' % env_cmd)
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

    # create list output files in the HTML output
    fileList = os.listdir(output_dir)
    fileList = [output_dir+"/"+filename for filename in fileList]
    for f in fileList:
        if os.path.isfile(f):
            shutil.copy(f, options.output_dir)

    # final output
    if options.vcffile != "None":
        fullFileList = glob.glob("%s/*.bcf" % output_dir)
        fullFileList.sort()
        typeLen=len(svtypes)
        for i in range(0, len(fullFileList)-1, typeLen):
            fileList=fullFileList[i:i+typeLen-1]
            baseName=fileList[0].split('/')[-1].split('_')[0]
            final_merge_cmd = "bcftools concat %s -a -O v -o %s/%s.vcf" % (" ".join(fileList), output_dir, baseName)
            print "\nfinal_cmd=%s" % final_merge_cmd
            subprocess.call(final_merge_cmd, shell=True)

        #### covert delly to metasv format

        tool_cmd = "export PATH=%s:\$PATH;%s/delly_to_metasv_run.py --fin INPUTFILE --fout OUTPUTFILE" % (bedtools_bin, delly2_bin)

        #if no stderr file is specified, we'll use our own
        stderr = tempfile.NamedTemporaryFile( prefix="TOOL-delly_to_metasv-stderr-", dir=tmp_dir )
        stdout = tempfile.NamedTemporaryFile( prefix="TOOL-delly_to_metasv-stdout-", dir=tmp_dir )

        # prepare command line
        swift_params = list()
        swift_params.append('-outputdir=' + output_dir)
        swift_params.append('-inputfiles=\"%s\"' %  ",".join(inputFiles))
        swift_params.append('-samplenames=\"%s\"' %  ",".join(sampleNames))

        ## construct the swift command
        swift_cmd = "%s -sites.file %s -tc.file %s %s " %   (swift_bin, sites_file, tc_file, swift_delly_to_metasv_file)
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
