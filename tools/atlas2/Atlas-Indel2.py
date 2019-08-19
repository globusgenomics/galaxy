#!/usr/bin/env python

import sys, os, commands, string, subprocess

galaxyhome=os.environ.get('GALAXY_HOME')

'''
n = len(sys.argv)
for i in range(0,n):
    print sys.argv[i]
'''
'''
   $input_BAM
   $input_FASTA
   $type

   $sample
   $p-cutoff
   $p-1bp-cutoff
   $input_BED
   $orig-base-qual
   $norm-base-qual
   $min-total-depth
   $min-var-reads
   $min-var-ratio
   $strand-dir-filter
   $near-read_end_ratio
   $homo-var-cutoff
   
   $output
'''

input_BAM = " -b " + sys.argv[1]
input_FASTA = " -r " + sys.argv[2]
model_type = " " + sys.argv[3]

if sys.argv[4]!="none":
    sample = " -s " + sys.argv[4]
else:
    sample = ""    

if sys.argv[5]!="none":
    p_cutoff = " -p " + sys.argv[5]
else:
    p_cutoff = ""   

if sys.argv[6]!="none":
    p1bp_cutoff = " -P " + sys.argv[6]
else:
    p1bp_cutoff = ""  

if sys.argv[7]!="none":
    input_BED = " -B " + sys.argv[7]
else:
    input_BED = ""  

if sys.argv[8]=="Yes":
    orig_base_qual = " -O "
else:
    orig_base_qual = ""  

if sys.argv[9]=="Yes":
    norm_base_qual = " -N "
else:
    norm_base_qual = ""  

if sys.argv[10]!="none":
    mintotal_depth = " -t " + sys.argv[10]
else:
    mintotal_depth = ""  

if sys.argv[11]!="none":
    minvar_reads = " -m " + sys.argv[11]
else:
    minvar_reads = ""  

if sys.argv[12]!="none":
    minvar_ratio = " -v " + sys.argv[12]
else:
    minvar_ratio = ""  

if sys.argv[13]=="Yes":
    strand_dirfilter = " -f "
else:
    strand_dirfilter = ""  

if sys.argv[14]!="none":
    near_read_end_ratio = " -n " + sys.argv[14]
else:
    near_read_end_ratio = ""  

if sys.argv[15]!="none":
    homo_var_cutoff = " -h " + sys.argv[15]
else:
    homo_var_cutoff = ""  

if sys.argv[16]!="none":
    site_list = " -a " + sys.argv[16]
else:
    site_list = ""

output = " -o " + sys.argv[17]


parameter = input_BAM + input_FASTA + output + model_type + sample + p_cutoff + p1bp_cutoff + input_BED + orig_base_qual + norm_base_qual + mintotal_depth + minvar_reads + minvar_ratio + strand_dirfilter + near_read_end_ratio + homo_var_cutoff
#print parameter


##uncomment when running on VM

os.chdir(galaxyhome + "/tools/atlas2/")     #uncomment when running on cluster
#os.chdir("/media/Work/galaxy-proteonics/tools/atlas2/")    #for local running

command="ruby ./Atlas-Indel2/Atlas-Indel2.rb " + parameter

#command = "ruby ./Atlas-Indel2/Atlas-Indel2.rb  -b /media/Work/galaxy-proteonics/database/files/000/dataset_47.dat -r /media/Work/galaxy-proteonics/database/files/000/dataset_18.dat -o /media/Work/galaxy-proteonics/database/files/000/dataset_209.dat -I"


print command

subprocess.call(command,shell=True)


