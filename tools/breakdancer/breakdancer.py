#!/usr/bin/python

import argparse, os, shutil, subprocess, sys, tempfile, shlex, time

parser = argparse.ArgumentParser(description='')
# required
parser.add_argument('-i1', dest='inputBamFile', required=True, help='the bam input file')
parser.add_argument('-o1', dest='outputRawFile', required=True, help='the raw output file')
# optional
parser.add_argument('-o2', dest='outputVcfFile', required=False, help='the vcf output file')
parser.add_argument('-o', dest='chromosome', required=False, help='operate on a single chromosome')
parser.add_argument('-s', dest='minLength', type=int, required=False, help='minimum length of a region', default=7)
parser.add_argument('-c', dest='cutoff', type=int, required=False, help='cutoff in unit of standard deviation', default=3)
parser.add_argument('-m', dest='maxSvSize', type=int, required=False, help='maximum SV size', default=1000000000)
parser.add_argument('-q', dest='minMapQuality', type=int, required=False, help='minimum alternative mapping quality', default=35)
parser.add_argument('-r', dest='minReadDepth', type=int, required=False, help='minimum number of read pairs required to establish a connection', default=2)
parser.add_argument('-x', dest='maxHaploidCov', type=int, required=False, help='maximum threshold of haploid sequence coverage for regions to be ignored', default=1000)
parser.add_argument('-b', dest='bufferSize', type=int, required=False, help='buffer size for building connection', default=100)
parser.add_argument('-t', dest='onlyTrans', action='store_true', help='only detect transchromosomal rearrangement', default=False)
parser.add_argument('-d', dest='prefix', required=False, help='prefix of fastq files that SV supporting reads will be saved by library')
parser.add_argument('-g', dest='bedFormat', required=False, help='dump SVs and supporting reads in BED format for GBrowse')
parser.add_argument('-l', dest='matePair', required=False, help='analyze Illumina long insert (mate-pair) library')
parser.add_argument('-a', dest='sortByLibrary', action='store_true', help='print out copy number and support reads per library rather than per bam', default=False)
# parser.add_argument('-h', dest='AFColumn', action='store_true', help='print out Allele Frequency column', default=False)
parser.add_argument('-y', dest='scoreFilter', type=int, required=False, help='output score filter', default=30)



# binPath = os.environ['BREAKDANCER_BIN']

# bam2cfgPath = binPath+"/bam2cfg.pl"
# breakdancer2vcfPath = binPath+"/breakdancer2vcf.py"


# def bam2cfg(args, tempDir):
# 	config = tempDir+"/breakdancer_config"
# 	cmd = 'perl %s %s' % (bam2cfgPath, args.inputBamFile)
# 	execute(cmd, output=config)
# 	return config


# def breakdancer(args, config):
# 	cmd = 'breakdancer-max %s' % (config)
# 	execute(cmd, output=args.outputRawFile)


# def breakdancer2vcf(args):
# 	cmd = "python %s -i %s -o %s" % (breakdancer2vcfPath, args.outputRawFile, args.outputVcfFile)
# 	execute(cmd)

def execute(cmd, output=None):
	# function to execute a cmd and report if an error occur
	# print(cmd)
	try: 
		process = subprocess.Popen(args=shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		process.wait()
		stdout,stderr = process.communicate()
	except Exception, e: # une erreur de ma commande : stderr
		sys.stderr.write("problem doing : %s\n%s\n" %(cmd, e))
		return
	if output:
		output = open(output, 'w')
		output.write(stdout)
		output.close()
	if stderr != '': # une erreur interne au programme : stdout (sinon, souvent des warning arrete les programmes)
		sys.stdout.write("warning or error while doing : %s\n-----\n%s-----\n\n" %(cmd, stderr))


# def execute(cmd, output=""):
# 	try: 
# 		err = open(tempDir+"/errorLog", 'a')
# 		if output != "":
# 			out = open(output, 'w')
# 		else:
# 			out = subprocess.PIPE
# 		process = subprocess.Popen(args=shlex.split(cmd), stdout=out, stderr=err)
# 		process.wait()
# 		err.close()
# 		if out != subprocess.PIPE:
# 			out.close()
# 	except Exception, e:
# 		sys.stderr.write("problem doing : %s\n" %(cmd))
# 		sys.stderr.write('%s\n\n' % str(e))


def check(output):
	if (os.path.getsize(output)>0):
		return True
	else:
		sys.stderr.write('The output file is empty : %s\n' % (output))


def getLine(file):
	try:
		f = open(file, 'r')
		lignes = f.readlines()
		n=0
		for ligne in lignes:
			if ligne.strip()[0]!="#":
				n+=1
		# n = len(lignes)
		f.close()
	except Exception, e:
		sys.stderr.write('%s\n' % str(e))
	return n


def compare(output1, output2):
	# compare le nombre de ligne entre 2 fichiers
	# pour verifier qu'aucune ligne n'a ete saute
	num_raw = getLine(output1)
	num_vcf = getLine(output2)
	if (num_raw==num_vcf):
		return True
	else:
		sys.stderr.write('Not the same number of variant between the raw file and the vcf file : %d vs %d\n' % (num_raw, num_vcf))

def getVersion(program):
	import tempfile, subprocess
	try:
		temp = tempfile.NamedTemporaryFile().name
		tempStdout = open(temp, 'wb')
		proc = subprocess.Popen(args=program, shell=True, stdout=tempStdout, stderr=tempStdout)
		tempStdout.close()
		returncode = proc.wait()
		stdout = None
		for line in open(tempStdout.name, 'rb'):
			if line.lower().find('version') >= 0:
				stdout = line.strip()
				break
		if stdout:
			sys.stdout.write('%s\n' % stdout)
	except:
		sys.stdout.write('Could not determine %s version\n' % (program))

def __main__():

	time0 = time.time()
	tempDir = tempfile.mkdtemp()
	args = parser.parse_args()
	getVersion("breakdancer-max")

	os.system('echo $PATH')
	try:
		# config = bam2cfg(args, tempDir)
		# breakdancer(args, config)
		# breakdancer2vcf(args)

		# bam2cfg
		config = tempDir+"/breakdancer_config"
		cmd = 'bam2cfg.pl %s' % (args.inputBamFile)
		execute(cmd, output=config)

		# breakdancer
		cmd = 'breakdancer-max %s' % (config)
		cmd += ' -s %d -c %d -m %d -q %d -r %d -x %d -b %d -y %d' % (args.minLength, args.cutoff, args.maxSvSize, args.minMapQuality, args.minReadDepth, args.maxHaploidCov, args.bufferSize, args.scoreFilter)
		if args.chromosome:
			cmd += ' -o %s ' % (args.chromosome)
		if args.onlyTrans:
			cmd += ' -t '
		if args.prefix:
			cmd += ' -d %s ' % (args.prefix)
		if args.bedFormat:
			cmd += ' -g %s ' % (args.bedFormat)
		if args.matePair:
			cmd += ' -l '
		if args.sortByLibrary:
			cmd += ' -a '
		# if args.AFColumn:
		# 	cmd += ' -h '
		execute(cmd, output=args.outputRawFile)

		# breakdancer2vcf
		if args.outputVcfFile:
			cmd = "breakdancer2vcf.py -i %s -o %s" % (args.outputRawFile, args.outputVcfFile)
			execute(cmd)

		#	quelques tests
		# verifier que les fichiers de sorties ne sont pas vides
		check(args.outputRawFile)
		check(args.outputVcfFile)
		# comparer le nombre de ligne entre les 2 fichiers
		# pour etre sur que toute les variations ont ete prises en compte
		compare(args.outputRawFile, args.outputVcfFile)
			
		sys.stdout.write('\nDone in %d seconds\n' %(int(time.time()-time0)))

		# if (os.path.getsize(errorFile)>0):
		# 	sys.stdout.write('At least one non fatal error was send, check the file : %s\n' %(errorFile))
		# 	try:
		# 		err = open(errorFile, 'r')
		# 		errors = err.read()
		# 		err.close()
		# 		sys.stdout.write("Errors :\n%s" %(errors))
		# 	except Exception, e:
		# 		sys.stderr.write('%s\n' % str(e))
		# else:
		sys.stdout.write('BreakDancer successful')

	except Exception, e:
		sys.stdout.write('BreakDancer fail\n')
		sys.stderr.write('%s\n' % str(e))
		sys.exit()

	finally:
		if os.path.exists(tempDir):
			shutil.rmtree(tempDir)

if __name__=="__main__":
	__main__()
