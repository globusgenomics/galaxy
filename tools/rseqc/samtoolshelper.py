import sys
import subprocess as sp
import os

# Creates the sorted and indexed bam/bai files that are requried for both bam2wig and RSEQC_count
def samtools_sorted(bam):
	sortedbam = bam + ".sorted"
	indexedbam = ".".join([sortedbam,"bam.bai"])
	sp.call(['samtools', 'sort', '-m 1000000000', bam, sortedbam])
	sortedbam = sortedbam + '.bam'
	sp.call(['samtools', 'index', sortedbam, indexedbam])
	return sortedbam

def main(args):
	args[2] = samtools_sorted(args[2])
	sp.call(args)


if __name__ == "__main__":
	main(sys.argv[1:])