#!/usr/bin/env python

import sys, os, commands, string, subprocess, shutil

galaxyhome=os.environ.get('GALAXY_HOME')


n = len(sys.argv)
for i in range(0,n):
    print sys.argv[i]


input_BAM = " -i " + sys.argv[1]
input_FASTA = " -r " + sys.argv[2]
platform = " " + sys.argv[3]

if sys.argv[4]!="none":
    target_region = " -t " + sys.argv[4]
else:
    target_region = ""    

if sys.argv[5]=="no":    #don't output VCF
    sample_name = ""  
else:
    sample_name = " -v -n " + sys.argv[5]   #output VCF, with sample name
    output_vcf = sys.argv[16]
 

if sys.argv[6]!="none":
    cutoff = " -c " + sys.argv[6]
else:
    cutoff = ""  

if sys.argv[7]!="none":
    min_coverage = " -y " + sys.argv[7]
else:
    min_coverage = ""  

if sys.argv[8]!="none":
    prior_prob_e = " -e " + sys.argv[8]
else:
    prior_prob_e = ""  

if sys.argv[9]!="none":
    prior_prob_l = " -l " + sys.argv[9]
else:
    prior_prob_l = ""  

if sys.argv[10]!="none":
    filter_m = " -m " + sys.argv[10]
else:
    filter_m = ""  

if sys.argv[11]!="none":
    filter_g = " -g " + sys.argv[11]
else:
    filter_g = ""  

if sys.argv[12]!="none":
    filter_f = " -f " + sys.argv[12]
else:
    filter_f = ""  

if sys.argv[13]!="none":
    filter_p = " -p " + sys.argv[13]
else:
    filter_p = ""  

input_BAI = sys.argv[14]

output = " -o " + sys.argv[15]


parameter = input_BAM + input_FASTA + output + platform + target_region + sample_name + cutoff + min_coverage + prior_prob_e + prior_prob_l + filter_m + filter_g + filter_f + filter_p
#print parameter

#Atlas-SNP2 needs an index file (.bai) that located in the same directory with the input BAM file. 
#The bai file will be copied to the same directory of input BAM file and renamed as xxx.bam.bai in background.
BAIfile = sys.argv[1] + ".bai"  #inputbamfile.bam.bai
if not os.path.exists(BAIfile):
    shutil.copyfile(input_BAI,BAIfile)

##uncomment when running on VM

os.chdir(galaxyhome + "/tools/atlas2/")     #uncomment when running on cluster
#os.chdir("/media/Work/galaxy-proteonics/tools/atlas2/")    #for local running

command="ruby ./Atlas-SNP2/Atlas-SNP2.rb " + parameter

#command = "ruby ./Atlas-SNP2/Atlas-SNP2.rb -i NA12275.bam -r ~/refs/human_g1k_v37.fasta -o"


print command

subprocess.call(command,shell=True)


#copy the generated vcf file to Galaxy's output_vcf file
if sys.argv[5]!="no":
    generate_vcf = sys.argv[15] + ".vcf"
    shutil.copyfile(generate_vcf,output_vcf)


