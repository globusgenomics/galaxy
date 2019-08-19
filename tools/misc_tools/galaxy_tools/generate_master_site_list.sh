#!/bin/bash

## parse command line args
config_file=$1
vcf_dir=$(grep VCF $config_file | awk '{print $2}')

output_file=$2

## pull ALL variant sites from experiment
for file in $(ls ${vcf_dir}/*.vcf); do
  echo $file
  grep -v '^#' $file | cut -f1,2 >> /tmp/atlas_tmp
done

## subset to unique sites
if [ `grep -c "chr" /tmp/atlas_tmp` -gt 1 ]; then 
    sort -V -s -k1.4,1.5 -k2 -k3 /tmp/atlas_tmp | uniq > $output_file;
else
    sort -V -s -k1 -k2 /tmp/atlas_tmp | uniq > $output_file;
fi
rm /tmp/atlas_tmp
