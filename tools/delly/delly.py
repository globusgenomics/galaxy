#!/usr/bin/python

import argparse, os, shutil, subprocess, sys, tempfile, shlex, vcf, pysam

parser = argparse.ArgumentParser(description='')
parser.add_argument ( '-b', dest='bam', help='the bam file', required=True )	
parser.add_argument ( '-t', dest='types', help='SV analysis type (DEL, DUP, INV, TRA)', nargs='+', required=True )
parser.add_argument ( '-q', dest='map_qual', help='min. paired-end mapping quality', default='1' )
parser.add_argument ( '-u', dest='min_map_qual', help='min. mapping quality for genotyping', default='1' )
parser.add_argument ( '-s', dest='mad_cutoff', help='insert size cutoff, median+s*MAD (deletions only)', default=9 )
parser.add_argument ( '-g', dest='genome', help='genome fasta file' )
parser.add_argument ( '-m', dest='min_flank', help='minimum flanking sequence size', default=13 )
parser.add_argument ( '-e', dest='epsilon', help='allowed epsilon deviation of PE vs. SR deletion', default=0.1 )
parser.add_argument ( '-o', dest='output', help='output file' )
parser.add_argument ( '-x', "--exclude", dest="exclude", help="exclude BED file" )
parser.add_argument ( '--output-tmp', dest='temp_dir', help='output temp dir' )

#dellyPath=os.path.dirname(os.path.realpath(__file__))+"/delly_v0.5.5"
#os.system("chmod +x "+str(dellyPath))
dellyPath="delly_v0.6.3_parallel"

def execute( cmd, output="" ):
	tmp_dir = tempfile.mkdtemp()
	try: 
		err = open(tmp_dir+"/errorLog", 'a')
		if output != "":
			out = open(output, 'w')
		else:
			out = subprocess.PIPE
		process = subprocess.Popen( args=shlex.split(cmd), stdout=out, stderr=err )
		process.wait()
		err.close()
		if out != subprocess.PIPE:
			out.close()
	except Exception, e:
		sys.stderr.write("problem doing : %s\n" %(cmd))
		sys.stderr.write( '%s\n\n' % str(e) )


def delly(args, tempDir):
	tempOutputs=[];
	for typ in args.types:
		output=str(tempDir)+"/"+str(typ)
		tempOutputs.append(output)
		cmd = "%s %s -t %s -o %s -q %s -s %s -m %s -e %s -u %s " % (dellyPath, args.bam, typ, output, args.map_qual, args.mad_cutoff, args.min_flank, args.epsilon, args.min_map_qual)
		if (args.genome!="" and args.genome!=None):
			cmd += " -g %s" % (args.genome)
                if (args.exclude != "" and args.exclude!=None):
                        cmd += " -x %s" % (args.exclude)
		print ("command = "+cmd)
		execute( cmd )
	return tempOutputs


def merge(outputs, outputFile):
	template = vcf.Reader(filename=outputs[0])
	vcfWriter = vcf.Writer(open(outputFile, 'w'), template)
	for output in outputs:
		vcfReader = vcf.Reader(filename=output)
		for record in vcfReader:
			vcfWriter.write_record(record)
	return 0


def getVersion(program):
	try:
		tmp = tempfile.NamedTemporaryFile().name
		tmp_stdout = open( tmp, 'wb' )
		proc = subprocess.Popen( args=program, shell=True, stdout=tmp_stdout )
		tmp_stdout.close()
		returncode = proc.wait()
		stdout = None
		for line in open( tmp_stdout.name, 'rb' ):
			if line.lower().find( 'version' ) >= 0:
				stdout = line.strip()
				break
		if stdout:
			sys.stdout.write( '%s\n' % stdout )
	except:
		sys.stdout.write( 'Could not determine %s version\n' % (program) )

def fix_header(outputF, sample_names, tempDir):
    inF = open(outputF, "r")
    outF = open("%s/tmp.vcf" % tempDir, "w")
    for line in inF:
        if line.startswith("#CHROM"):
            line = line.rstrip()
            cols = line.split("\t")
            new_cols = []
            for col in cols:
                if col in sample_names:
                    col = sample_names[col]
                new_cols.append(col)
            line = "\t".join(new_cols) + "\n"
        outF.write(line)
    inF.close()
    outF.close()
    shutil.move("%s/tmp.vcf" % tempDir, outputF)

def __main__():
	print(os.path.dirname(os.path.realpath(__file__)))
	args = parser.parse_args()

        if args.temp_dir is not None:
            tempDir = args.temp_dir
            if not os.path.exists(tempDir):
                os.mkdir(tempDir)
        else:
	    tempDir = tempfile.mkdtemp();

	getVersion(dellyPath)

        ## In our implementation we are assuming BAI files are already generated. 
        ## Will need to modify for future implementations
        input_dir = "%s/input" % tempDir
        if os.path.exists(input_dir):
            shutil.rmtree(input_dir)
        os.mkdir(input_dir)

	try:
	#    execute("samtools index " + args.bam)
            bamfiles = args.bam.split()
            linked_bamfiles = []
            sample_names = {}
            for bamfile in bamfiles:
                name = os.path.basename(bamfile)
                linkname = "%s/%s" % (input_dir, name)
                os.symlink(bamfile, linkname)
                linked_bamfiles.append(linkname)
                print "FOUND BAM FILE: %s" % linkname

                # create or link the bai index
                file_base = ('.').join(bamfile.split('.')[:-1])
                bai_file = None
                if os.path.exists("%s.bai" % file_base):
                   bai_file = "%s.bai" % file_base
                elif os.path.exists("%s.bai" % bamfile):
                   bai_file = "%s.bai" % bamfile
                else:
                   bai_file = None

                if os.path.exists(bai_file):
                    bai_link_name = "%s/%s.bai" % (input_dir, name)
                    os.symlink(bai_file, bai_link_name)
                else:
                    execute("samtools index " + linkname)

                # get the sample name for the bam file and store in dictionary
                samfile = pysam.AlignmentFile(linkname,'rb')
                sn = samfile.header['RG'][0]['SM']
                filebase = '.'.join(os.path.basename(linkname).split('.')[:-1]) 
                sample_names[filebase] = sn

            args.bam = " ".join(linked_bamfiles)

	except Exception, e:
		print("problem while indexing bam file " + str(e))

	try:
		tempOutputs = delly(args, tempDir)
	except Exception, e:
		sys.stdout.write("problem while runing delly " + str(e))

	try:	
		if args.output:
			outputs=[]
			for output in tempOutputs:
				if os.path.exists(output):
					if os.path.getsize(output)>0:
						outputs.append(output)
			if len(outputs)>0:
				merge(outputs, args.output)
                                fix_header(args.output, sample_names, tempDir)
	except Exception, e:
		sys.stdout.write("problem while merging files " + str(e))

	finally:
		try:
			if os.path.exists(tempDir):
				shutil.rmtree(tempDir)
                                #print "TempDir: %s" % tempDir
			#os.system("rm "+str(args.bam)+".bai")
		except:
			pass



if __name__=="__main__":
	__main__()
