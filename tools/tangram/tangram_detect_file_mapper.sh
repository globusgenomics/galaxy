#!/bin/sh

#take care of the mapper args
while [ $# -gt 0 ]; do
  case $1 in
    -bam)     bam=$2;;
    *) echo "$0: bad mapper args" 1>&2
       exit 1;;
  esac
  shift 2
done

samtools view -H ${bam} | grep "\@SQ" | sed 's/^.*SN://g' | cut -f 1 | grep -v GL | grep -v NC_ |  awk '{printf "[%d] %s\n", NR-1, $1}'

