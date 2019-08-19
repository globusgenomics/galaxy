#!/bin/bash

#take care of the mapper args
while [ $# -gt 0 ]; do
  case $1 in
    -bam)     bam=$2;;
    *) echo "$0: bad mapper args" 1>&2
       exit 1;;
  esac
  shift 2
done

counter=0
for i in $bam/*.bam; do
  samtools view -H ${i} | grep "\@RG" | sed 's/^.*SM://g' | awk '{printf "[%s] %s\n", '$counter', $1}';
  let counter=$counter+1; 
done
