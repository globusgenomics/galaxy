#!/usr/bin/env python

import sys, os, commands, string, subprocess, shutil

galaxyhome=os.environ.get('GALAXY_HOME')


n = len(sys.argv)
for i in range(0,n):
    print sys.argv[i]


pedfile = " -p " + sys.argv[1]
datfile = " -d " + sys.argv[2]
glfIndexFile = " -g " + sys.argv[3]

if sys.argv[4]!="none":
    cutoff = " -c " + sys.argv[4]
else:
    cutoff = ""    

if sys.argv[5]!="none":
    minMapQuality = " --minMapQuality " + sys.argv[5]
else:
    minMapQuality = ""   

if sys.argv[6]!="none":
    minDepth = " --minDepth " + sys.argv[6]
else:
    minDepth = ""    

if sys.argv[7]!="none":
    maxDepth = " --maxDepth " + sys.argv[7]
else:
    maxDepth = ""    

if sys.argv[8]!="none":
    minPercSampleWithData = " -minPercSampleWithData " + sys.argv[8]
else:
    minPercSampleWithData = ""    

if sys.argv[9]!="none":
    theta = " --theta " + sys.argv[9]
else:
    theta = ""    

if sys.argv[10]!="none":
    poly_tstv = " --poly_tstv " + sys.argv[10]
else:
    poly_tstv = ""    


if sys.argv[11]=="no":    #don't turn on de novo mutation detection
    denovo = ""
    rate_denovo = ""  
    tstv_denovo = ""  
    minLLR_denovo = ""  
else:
    denovo = " --denovo "
    if sys.argv[11]!="none":    #turn on de novo mutation
        rate_denovo = " --rate_denovo " + sys.argv[11]   
    else: 
        rate_denovo = ""
    if sys.argv[12]!="none":    #turn on de novo mutation
        tstv_denovo = " --tstv_denovo " + sys.argv[12]   
    else: 
        tstv_denovo = ""
    if sys.argv[13]!="none":    #turn on de novo mutation
        minLLR_denovo = " --minLLR_denovo " + sys.argv[13]   
    else: 
        minLLR_denovo = ""


if sys.argv[14]!="none":
    prec = " --prec " + sys.argv[14]
else:
    prec = ""  

if sys.argv[15]!="none":
    nthreads = " --nthreads " + sys.argv[15]
else:
    nthreads = ""  

if sys.argv[16]!="none":
    chr2process = " --chr2process " + sys.argv[16]
else:
    chr2process = ""  

if sys.argv[17]=="true":
    gl_off = " --gl_off " 
else:
    gl_off = ""  


output = " --vcf " + sys.argv[18]


parameter = pedfile + datfile + glfIndexFile + output + cutoff + minMapQuality + minDepth + maxDepth + minPercSampleWithData + theta + poly_tstv + denovo + rate_denovo + tstv_denovo + minLLR_denovo + prec + nthreads + chr2process + gl_off
print parameter



##uncomment when running on VM

os.chdir(galaxyhome + "/tools/polymutt/")     #uncomment when running on cluster
#os.chdir("/media/Work/galaxy-proteonics/tools/polymutt/")    #for local running

command="./polymutt.0.04/bin/polymutt " + parameter

#command = "./polymutt.0.04/bin/polymutt -p test.ped -d test.dat -g test.gif -c 0.9 --minDepth 150 --maxDepth 200 --nthreads 4  --vcf test.vcf"


print command

subprocess.call(command,shell=True)


