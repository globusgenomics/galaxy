#!/bin/bash
##<fastq1 in> <fastq1 out> <fastq2 in> <fastq2 out> <singletons reads out> <pe reads output>
## $5 and $6 are outputs
perl /usr/local/velvet/contrib/select_paired/select_paired.pl $1 $2 $3 $4 $5
perl /usr/local/velvet/shuffleSequences_fasta.pl $2 $4 $6
rm -f $2
rm -f $4 
