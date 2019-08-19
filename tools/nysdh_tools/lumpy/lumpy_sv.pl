#!/usr/bin/perl

### script that first extract the chimeric reads from the bam file (i.e. split reads and
### insert size distribution
## then these two output are feed into lumpy to identified location of potential
## genomic insertions


use strict;

my $path = $ARGV[0];
#my $samtools = "/scratch/galaxy/test/pxl10/Projects/MTB_pipeline_upgrade/Pipeline_shared_tools/Mapping/software/samtools";
my $pairend_distro = "pairend_distro.py";
my $extractSplitReads = "extractSplitReads_BwaMem";
my $lumpy = "lumpy";



system ("samtools view $path/tmp_fixed_sorted_marked_RG_realigned.bam | tail -n+100000 | $pairend_distro -r 150 -X 4 -N 10000 -o $path/sample.pe.histo 2> /dev/null");
system ("samtools view -u -F 0x0002  $path/tmp_fixed_sorted_marked_RG_realigned.bam | samtools view -u -F 0x0100 - | samtools view -u -F 0x0004 - | samtools view -u -F 0x0008 - | samtools view -b -F 0x0400 - > $path/sample.discordant.pe.bam 2> /dev/null");
system ("samtools view -h $path/tmp_fixed_sorted_marked_RG_realigned.bam | $extractSplitReads -i stdin | samtools view -Sb - > $path/sample.sr.bam 2> /dev/null");
system ("$lumpy -mw 40 -tt 0.0 -pe bam_file:$path/sample.discordant.pe.bam,histo_file:$path/sample.pe.histo,mean:500,stdev:50,read_length:150,min_non_overlap:150,discordant_z:4,back_distance:20,weight:1,id:1,min_mapping_threshold:20 -sr bam_file:$path/sample.sr.bam,back_distance:20,weight:1,id:2,min_mapping_threshold:20 -b > $path/sample_svn.out 2> /dev/null");

exit;
