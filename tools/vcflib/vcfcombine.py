#!/usr/bin/env python

import argparse
import subprocess as sp
##import subprocess, tempfile, sys, os

def __main__():
	'''
	Combine an arbitrary number of VCF files using vcflib's vcfcombine.

	Author: vasa
	'''
	#Parse Command Line
	parser = argparse.ArgumentParser()
	parser.add_argument('--vcf', nargs='*', dest='vcfs')
	parser.add_argument('--out', dest='out')
	args = parser.parse_args()

	## construct command
	cmd = ['vcfcombine']
	for vcf in args.vcfs: cmd.append(vcf)
##        cmd.append(" > ")
##        cmd.append(args.out)

##        cmd_line = " ".join(cmd)

	## run command, redirecting stdout to file
	print cmd
	with open(args.out, 'w') as outFile:
		sp.call(cmd, stdout=outFile, shell=True)

##        stderr_name = tempfile.NamedTemporaryFile( prefix = "stderr" ).name
##        proc = subprocess.Popen( args=cmd_line, shell=True, stderr=open( stderr_name, 'wb' ) )
##        exit_code = proc.wait()
##        if exit_code:
##            for line in open( stderr_name ):
##               print >> sys.stderr, line
##                os.unlink( stderr_name ) #clean up
##                raise Exception( "Error creating consensus sequences " )
##            os.unlink( stderr_name ) #clean up


if __name__ == "__main__":
	__main__()
