#!/usr/bin/perl -w

use strict;

my $path = $ARGV[2];
my $prefix = $ARGV[3];
my $downsampling_factor =$ARGV[4];


open (OUT,">$path/Stats.out");
unless ((-e $ARGV[0]) && (-e $ARGV[1])) {
	print OUT "Missing input read file(s) for mapping $prefix\n";
	print OUT "\nFAILED RUN\n";
	close OUT;
	exit;
}


#my $bwa = "/scratch/galaxy/test/pxl10/Projects/MTB_pipeline_upgrade/Pipeline_shared_tools/Mapping/software/bwa";
#my $samtools = "/scratch/galaxy/test/pxl10/Projects/MTB_pipeline_upgrade/Pipeline_shared_tools/Mapping/software/samtools";
my $picard = "/mnt/galaxyTools/tools/picard/1.129/picard.jar";
my $GATK = "/mnt/galaxyTools/tools/gatk3/3.3/GenomeAnalysisTK.jar";
my $bwa_ref = "/mnt/galaxyIndices/genomes/Mtuberculosis/H37Rv/bwa/REFERENCE";

## map reads
system ("bwa mem $bwa_ref -t 5 $ARGV[0] $ARGV[1] > $path/tmp.sam");
my $down=0;
if ($downsampling_factor < 1) {
	$down=1;
	system ("samtools view -s $downsampling_factor -Sb $path/tmp.sam > $path/tmp.bam");
} else {
	system ("samtools view -Sb $path/tmp.sam > $path/tmp.bam");
}


system ("samtools sort $path/tmp.bam $path/tmp_fixed_sorted");
system ("java -jar $picard MarkDuplicates I=$path/tmp_fixed_sorted.bam O=$path/tmp_fixed_sorted_marked.bam M=$path/metrics.out REMOVE_DUPLICATES=true");
system ("java -jar $picard AddOrReplaceReadGroups I=$path/tmp_fixed_sorted_marked.bam O=$path/tmp_fixed_sorted_marked_RG.bam SORT_ORDER=coordinate RGID=RAL357 RGLB=RAL357 RGPL=illumina RGPU=RAL357 RGSM=RAL357 CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT");
#unless (-e "$path/reference.dict") {
#	system ("java -jar $picard CreateSequenceDictionary R=$path/reference.fasta O=$path/reference.dict");
#}

#system ("samtools faidx $path/reference.fasta");

##GATK RealignerTargetCreator & IndelRealigner step
my $reference_fasta = "/mnt/galaxyIndices/genomes/Mtuberculosis/H37Rv/seq/H37rv.fasta";

system ("java -jar $GATK -T RealignerTargetCreator -R $reference_fasta -I $path/tmp_fixed_sorted_marked_RG.bam -o $path/indels.intervals -et NO_ET -K /opt/galaxy/tools/nysdh_tools/mapping/pascal.lapierre_health.ny.gov.key");
system ("java -jar $GATK -T IndelRealigner -R $reference_fasta -I $path/tmp_fixed_sorted_marked_RG.bam -targetIntervals $path/indels.intervals -o $path/tmp_fixed_sorted_marked_RG_realigned.bam -et NO_ET -K /opt/galaxy/tools/nysdh_tools/mapping/pascal.lapierre_health.ny.gov.key");
system ("samtools flagstat $path/tmp_fixed_sorted_marked_RG_realigned.bam > $path/flagstat.out");




## get mapping Stats ##
my $commandline = $0 . " ". (join " ", @ARGV);
if ($down == 1) {
	print OUT "\n---------------------------------------------------\n$commandline\n\nMapping Statistics (Downsampled ratio $downsampling_factor) :\n\n";
} else {
	print OUT "\n---------------------------------------------------\n$commandline\n\nMapping Statistics:\n\n";
}

open (FILE,"$path/flagstat.out");
my $flagstat_line=<FILE>;
my @split_flagstat_line = split (' ',$flagstat_line);
print OUT "Total # of mapped reads :\t$split_flagstat_line[0]\n";
$flagstat_line=<FILE>;$flagstat_line=<FILE>;
@split_flagstat_line = split ('\(',$flagstat_line);
my @split_at_percent = split ('\:',$split_flagstat_line[1]);


print OUT "Percent of mapped reads:\t$split_at_percent[0]\n";

close FILE;
close OUT;


exit;


#######################
