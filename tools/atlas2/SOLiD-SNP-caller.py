#!/usr/bin/env python

import sys, os, commands, string, subprocess

galaxyhome=os.environ.get('GALAXY_HOME')

'''
n = len(sys.argv)
for i in range(0,n):
    print sys.argv[i]
'''
'''
SOLiD-SNP-caller [in.bam] [ref.fa] [.bed region] > [output.vcf]
'''

input_BAM = " " + sys.argv[1]
input_FASTA = " " + sys.argv[2]

if sys.argv[3]!="none":
    input_BED = " " + sys.argv[3]
else:
    input_BED = ""    

output = " > " + sys.argv[4]


parameter = input_BAM + input_FASTA + input_BED + output
#print parameter


##uncomment when running on VM

os.chdir(galaxyhome + "/tools/atlas2/")     #uncomment when running on cluster
#os.chdir("/home/liubo/")    #for local running

command="./SOLiD-SNP-caller/SOLiD-SNP-caller " + parameter

#command = "SOLiD-SNP-caller [in.bam] [ref.fa] [.bed region] > [output.vcf]"


print command

subprocess.call(command,shell=True)


