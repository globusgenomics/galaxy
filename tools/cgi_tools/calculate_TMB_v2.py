##===========================================================================##
## calculate_TMB_v2.py
## April, 25, 2018
## Version 2.0
## 1. # calculate the TMB given a filtered VCF and a targeted BED file
## Usage: python ./calculate_TMB.py <target_bed_file.bed> <filtered_vcf_file.vcf>
##===========================================================================##

#!/usr/bin/python
import os, sys, subprocess
#Usage: python ./calculate_TMB.py <target_bed_file.bed> <filtered_vcf_file.vcf>

# calculate the TMB given a filtered VCF and a targeted BED file
target_bed = sys.argv[1]
region_size = 0
fh = open(target_bed, "r")
for line in fh:
    values = line.split("\t")
    region_size += int(values[2]) - int(values[1])
#print region_size
region_size_mb = float(region_size) / float(1000000)
#print region_size_mb
fh.close()

vcf = sys.argv[2]
vcf_filename = vcf.split('/')[-1].split('.vcf')[0]
bed_filename = target_bed.split('/')[-1].strip('.bed')

#print vcf_filename

## Edit : April 25, 2018 add bedtools intersect to get variants restrict by BED file
## create target vcf file based on BED
target_vcf = os.path.join(os.getcwd(),vcf_filename+'_'+bed_filename+'.vcf')
bedtools_cmd = 'bedtools intersect' + ' -a ' + vcf + ' -b ' + target_bed + ' -header > ' + target_vcf
#print '==== make bedtools file ====='
subprocess.call(bedtools_cmd, shell=True)

fh_vcf = open(target_vcf, "r")
var_size = 0
for line in fh_vcf:
    values = line.split("\t")
    if line.startswith("#"):
        continue

    ## this calculate the number of mutation,
    var_size += len(values[3]) + len(values[4]) - 1
print "TMB for %s in %s is: %s total variants / %s total region (MB) = %s" % ( vcf_filename, bed_filename, str(var_size), str(region_size_mb), str(float(var_size)/float(region_size_mb)) )

