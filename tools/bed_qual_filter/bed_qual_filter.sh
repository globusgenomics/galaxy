#!/bin/bash

## input vcf file
vcf=$1
## output vcf file
out=$2
## bed file to filter on
bed=$3

## filter on min QUAL and on target
/mnt/galaxyTools/tools/vcftools_0.1.11/bin/vcftools --vcf $vcf --bed $bed --recode-to-stream --minQ 10 |\
/mnt/galaxyTools/tools/vcflib/bin/vcffixup |\
/mnt/galaxyTools/tools/vcflib/bin/vcffilter -f "AC > 0" > $out

