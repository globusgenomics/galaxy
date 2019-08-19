#!/bin/bash
##<fastq1 in> <fastq1 out> <fastq2 in> <fastq2 out> <singletons reads out> <pe reads output>
## $5 and $6 are outputs
#perl /nfs/software/galaxy/tools/sr_assembly/fastq_pe_even.pl $1 $2 $3 $4 $5
perl /nfs/software/galaxy/tools/sr_assembly/shuffleSequences_fastq.pl $1 $3 $6
rm -f $2
rm -f $4 
