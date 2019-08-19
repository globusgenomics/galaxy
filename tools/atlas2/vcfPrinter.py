#!/usr/bin/env python
import sys, os, commands, string, subprocess, shutil

galaxyhome=os.environ.get('GALAXY_HOME')

'''
n = len(sys.argv)
for i in range(0,n):
    print sys.argv[i]
'''

input_VCF = sys.argv[1]
input_BAM = sys.argv[2]
input_FASTA = " -r " + sys.argv[3]

if sys.argv[4]!="none":
    n = " " + sys.argv[4]
else:
    n = ""

if sys.argv[5]!="none":
    p = " " + sys.argv[5]
else:
    p = ""

outputVCF = sys.argv[6]

output_FOLDER = outputVCF.split('.')[0]+'_files'

if not os.path.exists(output_FOLDER):
	os.makedirs(output_FOLDER)
else:
	pass


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


