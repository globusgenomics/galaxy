#!/bin/bash

vcf=$1
newvcf=$2
bedFile="$3"
keepSamples="$4"
refFile="$5"

## prototype filter
echo 'vcfkeepgeno $vcf GT | vcffilter -f "QUAL > 10" | vcfintersect -b $bedFile > tmp.vcf'
vcfkeepgeno $vcf GT | vcffilter -f "QUAL > 10" | vcfintersect -b $bedFile > tmp.vcf

echo 'vt normalize tmp.vcf -r $refFile -o tmp2.vcf 2>&1'
vt normalize tmp.vcf -r $refFile -o tmp2.vcf 2>&1

echo 'cat tmp2.vcf |   grep -v "/2" | grep -v "2/" | vcfsort   > tmp3.vcf'
cat tmp2.vcf |   grep -v "/2" | grep -v "2/" | vcfsort   > tmp3.vcf

## subset samples and filter monomorphic variants
echo 'vcftools --vcf tmp3.vcf --keep $keepSamples --recode --out tmp4 2>&1'
vcftools --vcf tmp3.vcf --keep $keepSamples --recode --out tmp4 2>&1

echo 'vcffixup tmp4.recode.vcf | vcffilter -f "AC > 0" > $newvcf'
vcffixup tmp4.recode.vcf | vcffilter -f "AC > 0" > $newvcf

## cleanup log and vestigial files
#rm tmp.vcf tmp2.vcf tmp3.vcf tmp4.recode.vcf tmp4.log
