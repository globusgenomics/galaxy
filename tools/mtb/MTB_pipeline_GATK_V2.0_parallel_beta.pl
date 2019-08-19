#!/usr/bin/perl

##### on September 7th, edited out line $pm -> set_waitpid_blocking_sleep(0); in an attempt to solve the missing sample issue

my $version="2.0";
use strict;
use Benchmark;
use LWP::Simple;
use Parallel::ForkManager;
use File::Temp;

#$pm->set_waitpid_blocking_sleep(0);


### changelog for MTB_pipeline_GATK###
## Version 1.1 ###
## Add output for access database upload ##
## Version 1.1.1 ###
## Fix bug for when no MTB reads are found in the samples to avoid crashes (lines 149)


### changelog for MTB_pipeline_GATK_V1_parallel.pl###
### modified script from MTB_pipeline_GATK_V1.1.1 to be run by script Parallel_TB_runner_vX.pl
### Removed necessary input arguments now taken in charged by Parallel_TB_runner
### No pipeline lock necessary
### usage perl MTB_pipeline_GATK_V1_parallel.pl $Read_R1 $tmpdir
### Add step to detect africanum signature sequence

### get current date ###
my ($day, $month, $year) = (localtime)[3,4,5];
my $date = sprintf("%04d_%02d_%02d", $year+1900, $month+1, $day);

## ARGV[0] = path to file
## ARGV[1] = path where to store running files

#my $path = $ARGV[1];  ## path to the temporary directory where to store the generated files
my $path = File::Temp->newdir;
my $log_file = $ARGV[1];    ## path to log file

#my $Read1 = $ARGV[2];   ## path to input directory object
#my $sample_name = $ARGV[3]; ## sample name of the fastq files
my $output_directory = $ARGV[2];    ## path to output_directory

my $reference = "/mnt/galaxyIndices/genomes/Mtuberculosis/H37Rv/seq/H37rv.fasta"; ## location of the reference genome
my $kraken_ref = "/mnt/galaxyIndices/genomes/Mtuberculosis/Kraken/Myco";
my $MTB_mask_region_file = "/opt/galaxy/tools/nysdh_tools/SNP_caller/MTB_mask_region.txt";

## scripts and software locations
my $spoligo_script = "/opt/galaxy/tools/nysdh_tools/spoligotyper/spoligo_kmers_parrallel.pl";
#my $kraken = "/mnt/galaxyTools/tools/kraken/0.10.5-beta/kraken";
#my $kraken_report = "/mnt/galaxyTools/tools/kraken/0.10.5-beta/kraken-report";
#my $bwa = "/mnt/galaxyTools/tools/bwa/0.7.12/bwa";
my $map_read_script = "/opt/galaxy/tools/nysdh_tools/mapping/Map_reads.pl";
#my $samtools = "/scratch/galaxy/test/pxl10/Projects/MTB_pipeline_upgrade/Pipeline_shared_tools/Mapping/software/samtools";
my $lumpy_sv_script = "/opt/galaxy/tools/nysdh_tools/lumpy/lumpy_sv.pl";
my $picard = "/mnt/galaxyTools/tools/picard/1.129/picard.jar";
my $GATK = "/mnt/galaxyTools/tools/gatk3/3.3/GenomeAnalysisTK.jar";
my $annotator_script = "/opt/galaxy/tools/nysdh_tools/SNP_caller/annotator_GATK_v2.01.pl";
my $clims_script = "/opt/galaxy/tools/nysdh_tools/SNP_caller/clims_report_v1.0.pl";


## generate subdirectories
if (!-d "$path/VCF") {
	system ("mkdir $path/VCF");
}
if (!-d "$path/Reports") {
	system ("mkdir $path/Reports");
}
if (!-d "$path/Consensus") {
	system ("mkdir $path/Consensus");
}
if (!-d "$path/CLIMS") {
	system ("mkdir $path/CLIMS");
}


## get reference genome size ##
open (FILE,$reference);
my @in_file=<FILE>;
close FILE;
my $join_file = join ('',@in_file);
my @split_reference_file = split ('\n',$join_file);
shift @split_reference_file;  ##remove header line
my $join_sequence = join ('',@split_reference_file);
my $gsize=length($join_sequence); ## genome size

## set a CLIMS report where data will by append
my $clims_filename="CLIMS_report\_$date.out";
open (CLIMS,">$path/$clims_filename");
close CLIMS;


my $Read_R1 = $ARGV[0];
my $Read_R2 = $Read_R1;
$Read_R2 =~ s/\_R1/\_R2/g;

if ((! -e $Read_R1) || ($Read_R1 !~ 'R1\_001\.fastq\.gz'))  {
	die "$Read_R1 not found, no sample processed.";
}
if ((! -e $Read_R2) || ($Read_R2 !~ 'R2\_001\.fastq\.gz'))  {
	die "$Read_R2 not found, no sample processed.";
}

my $stop=0;  ## variable used as a check point to tell if one of the step along the pipeline went wrong and pipeline needs to be stop if stop ==1

#system ("cp $reference $path/reference.fasta");
#system ("chmod 775 $path/reference.fasta");


## store start time of the script
my $t0 = Benchmark->new;

########################

my @split_read_path = split ('\/',$Read_R1);
my @split_read_filename = split ('\_',$split_read_path[-1]);
my $prefix = $split_read_filename[0]; ## sample name
my $local_log = "$path/Reports/$prefix\_$date.txt";

#open (LOG,">$path/Reports/$prefix\_$date.txt");
open(LOG, ">$local_log");
print LOG "Script ./MTB_pipeline_GATK_V2.0_parallel.pl Version $version\n------------------------------------------------\n\n";
print LOG "Pipeline script $0, Sample $prefix ($Read_R1), rundate $date\n\n";


#####  RUN Kraken  #####
print "\nkraken --preload --db $kraken_ref --gzip-compressed --fastq-input $Read_R1 $Read_R2 --paired --threads 2 --output $path/kraken.out\n\n";
print "\nkraken-report -db $kraken_ref $path/kraken.out > $path/kraken_report.out\n\n";
system ("kraken --preload --db $kraken_ref --gzip-compressed --fastq-input $Read_R1 $Read_R2 --paired --threads 2 --output $path/kraken.out 2> /dev/null");
system ("kraken-report -db $kraken_ref $path/kraken.out > $path/kraken_report.out 2> /dev/null");

print "\n\nkraken completed!\n\n";
## Parse Kraken report ##
my @taxonomy=();
my $id;
my $pid;
my $downsampling_factor;
my $top_match;

open (FILE,"$path/kraken_report.out");
my $in_kraken_file=<FILE>;
$in_kraken_file =~ s/ +/ /;  ## Cant remember why I do this...
my @split_kraken_line = split ('\t',$in_kraken_file);

print "\nin_kraken_file: $in_kraken_file\n";
if ($in_kraken_file =~ 'U') {  ## get number of unclassified reads
	print LOG "Taxonomic report:\n\nUnclassified reads:\t$split_kraken_line[0]%";
	if ($split_kraken_line[0] >= 10) { ## if more than 10% of undetermined reads
		print LOG "\tWarning : Possible presence of nontuberculous mycobacterium (NTM) organisms in the sample";
	}
	print LOG "\n\n";


	my $total_reads = $split_kraken_line[1]; ## get total unclassified number of reads in sample
	$in_kraken_file=<FILE>;
	$in_kraken_file =~ s/ +/ /;
	@split_kraken_line = split ('\t',$in_kraken_file);
	$total_reads += $split_kraken_line[1]; ## add unclassified and classified number of reads in sample
	$id="";
	$pid=0;

	do { ## iterate through kraken report until reaching line with Mycobacterium
		$in_kraken_file =<FILE>;
		$in_kraken_file =~ s/ //g;
		@split_kraken_line = split ('\t',$in_kraken_file);
	}
	until (($split_kraken_line[4] == 1763) || (eof));  ## 1763 is taxonomynomic ID for Mycobacterium sp.

	$in_kraken_file =<FILE>; ## go down one taxonomic level in the report
	chomp $in_kraken_file;
	$in_kraken_file =~ s/ +/ /;
	@split_kraken_line = split ('\t',$in_kraken_file);
	$split_kraken_line[5] =~ s/ +//;
	print LOG "$split_kraken_line[5]\t$split_kraken_line[0]%\n";
	$pid=$split_kraken_line[0]; ## percent read matches this id

	$downsampling_factor = 0;  ## variable used to set read downsampling factor. If 0, no downsampling.
	if ($split_kraken_line[1] > 0) { ## fix for version 1.1.1
		$downsampling_factor = (2000000/$split_kraken_line[1])/2;
	}

	$id=$split_kraken_line[5]; ## name of the top species match in kraken
	my %hash_species=(); ## hash to store other matches


	while ($in_kraken_file=<FILE>) { ## iterate through kraken report to get other matches
		chomp $in_kraken_file;
		$in_kraken_file =~ s/  //g;
		@split_kraken_line = split ('\t',$in_kraken_file);
		if (($split_kraken_line[0] >= 2)) {  ## % ID thershold in kraken
			print LOG "$split_kraken_line[5]\t$split_kraken_line[0]%\n";
			$hash_species{$split_kraken_line[5]} = $split_kraken_line[0];
		}
	}

	foreach my $key  (sort { $hash_species{$b} <=> $hash_species{$a} } keys %hash_species ) { ## sort the hash
  		push (@taxonomy,$key); ## generate a list of all the matches
  	}

	my $taxonomy_list=join ",",@taxonomy;	## this variable might be unused now
	$top_match=$taxonomy[0];


	if (exists $taxonomy[0]) {
		if ($taxonomy[0] =~ 'Mycobacterium tuberculosis') {
			$top_match="Mycobacterium tuberculosis";
		}
	} else {
		$top_match="No species identified";
	}

### check for RD1 region      ###
### RD1 is used to type bovis vs tuberculosis or africanum ###

	if ($top_match =~ 'bovis|tuberculosis|africanum') {
		my $count_read_matches=0;
		my $query = "CACTCTGAGAGGTTGTCA"; ## signature sequence 1
		my $reverse_query = reverse $query;
		$reverse_query =~ tr/ATCG/TAGC/;
		#system ("zcat $Read_R1 $Read_R2 | agrep -1 -c $query > $path/count.out");
		system ("zcat $Read_R1 $Read_R2 > $path/a.txt; agrep -1 -c $query $path/a.txt > $path/count.out; rm $path/a.txt");
		open (FILE_COUNT,"$path/count.out");
		my $count=<FILE_COUNT>;
		close FILE_COUNT;
		chomp $count;
		$count_read_matches = $count;
		#system ("zcat $Read_R1 $Read_R2 | agrep -1 -c $reverse_query > $path/count.out");
  		system ("zcat $Read_R1 $Read_R2 > $path/b.txt; agrep -1 -c $reverse_query $path/b.txt > $path/b_tmp.out; cat $path/b_tmp.out > $path/count.out; rm $path/b.txt; rm $path/b_tmp.out");
		open (FILE_COUNT,"$path/count.out");
		$count=<FILE_COUNT>;
		close FILE_COUNT;
		chomp $count;
		$count_read_matches  += $count;

		if ($count_read_matches >= 5) { ## if number of matches >= threshold
			if ($top_match =~ 'tuberculosis') {
				$query = "AAACGTCACGAACGTAACCCCAG"; ## signature sequence 2
				$reverse_query = reverse $query;
				$reverse_query =~ tr/ATCG/TAGC/;
				system ("zcat $Read_R1 $Read_R2 | grep -E '$query|$reverse_query' -c > $path/count.out");
				#print "\n\nzcat $Read_R1 $Read_R2 > $path/c.txt; agrep -e '$query;$reverse_query' -c $path/c.txt > $path/c_tmp.out; cat $path/c_tmp.out > $path/count.out; rm $path/c.txt; rm $path/c_tmp.out\n\n";
				#system ("zcat $Read_R1 $Read_R2 > $path/c.txt; agrep -e '$query;$reverse_query' -c $path/c.txt > $path/c_tmp.out; cat $path/c_tmp.out > $path/count.out; rm $path/c.txt; rm $path/c_tmp.out");
				open (FILE_COUNT,"$path/count.out");
				$count=<FILE_COUNT>;
				close FILE_COUNT;
				chomp $count;
				if ($count >= 5) {
					print LOG "\nAfricanum signature detected\n";
					$top_match ="Mycobacterium africanum";
				}

			} elsif ($top_match =~ 'bovis') {
				print LOG "\nRD1 region detected, likely Mycobacterium bovis\n";
				$top_match="Mycobacterium bovis";
			}
		} else {
			if ($top_match =~ 'tuberculosis|africanum|bovis') {
				print LOG "\nRD1 region absent, likely Mycobacterium bovis BCG\n";
				$top_match="Mycobacterium bovis BCG";
			}
		}
	}


	print LOG "\n\nMycobacterium Species detected : $top_match\n\n";

	if ($top_match =~ 'aviumx|canettii|abscessus|fortuitum|intracellulare|indicus|parascrofulaceum|kansasii|smegmatis|xenopi|thermoresistibile|phlei|rhodesiae|colombiense|yongonense') {
		print LOG "Species other than MTB detected. Pipeline Stop\n\n";
		close LOG;
		close FILE;
		exit;
	}


} else {
	print LOG "Failed Kraken Step\n\n";
	close FILE;
	exit;
}
#system ("rm $path/kraken.out");
#system ("rm $path/kraken_report.out");

if ($id !~ 'Mycobacterium tuberculosis complex') {
	print LOG "\nPipeline stopped due to insufficient TB reads in the sample\n\n";
	close LOG;
	close FILE;
	exit;
} else {
	close LOG;
	close FILE;
}
######## End Kraken #############


######## get spoligotype ########
my $spoligo="";
my $in_spoligo_file;
my @split_spoligo_file;
print "\n\n$spoligo_script $Read_R1 $Read_R2 $path 2> /dev/null\n\n";
print "\n\ncat $path/spoligo.out >> $path/Reports/$prefix\_$date.txt\n\n";

system ("$spoligo_script $Read_R1 $Read_R2 $path 2> /dev/null");
system ("cat $path/spoligo.out >> $path/Reports/$prefix\_$date.txt");
open (FILE,"$path/spoligo.out");
$in_spoligo_file=<FILE>;$in_spoligo_file=<FILE>;$in_spoligo_file=<FILE>;$in_spoligo_file=<FILE>;$in_spoligo_file=<FILE>;
if ($in_spoligo_file =~ 'Failed') {
	$spoligo = "FAIL\t";
} elsif ($in_spoligo_file =~ 'Unknown') {
	@split_spoligo_file = split ('\t',$in_spoligo_file);
	$spoligo = "PASS\t$split_spoligo_file[2] (Unknown type)";
} else {
	$in_spoligo_file =~ s/NY\_//g;
	@split_spoligo_file = split ('\t',$in_spoligo_file);
	$spoligo = "PASS\t$split_spoligo_file[4] ($split_spoligo_file[3])";
}
close FILE;

print "\n\nspoligo completed!\n\n";
# system ("rm $path/spoligo.out");
####### End Spolygotype ########


####### Read Mapping Step ######
#system ("bwa index -p $path/REFERENCE -a is $path/reference.fasta 2> /dev/null");
#print "\nreference generated\n";
print "\n$map_read_script $Read_R1 $Read_R2 $path $prefix $downsampling_factor\n";
system ("$map_read_script $Read_R1 $Read_R2 $path $prefix $downsampling_factor 2> /dev/null");
print "\nmapping completed\n";
## Get average depth and % reference coverage
system ("samtools depth -q 20 -Q 20 $path/tmp_fixed_sorted_marked_RG_realigned.bam |  awk \'{sum+=\$3} END { print \"Genome Coverage =\t\",NR/$gsize} END { print \"Average Depth =\t\",sum/$gsize}\' >> $path/Stats.out 2> /dev/null");
print "\ndepth coverage calculated\n";
system ("cat $path/Stats.out >> /$path/Reports/$prefix\_$date.txt");
print "\nreport generated\n";
####### End Mapping ########

## check for coverage  > 40x ##
print "\n\ncheck for coverage\n\n";
my $depth = 0;
open (FILE,"$path/Stats.out");
while (my $in_stat_file=<FILE>) {
	if ($in_stat_file =~ 'Average Depth') {
		chomp $in_stat_file;
		my @split_stat_file = split ('\t',$in_stat_file);
		$depth = $split_stat_file[1];
	}
}
close FILE;


if ($depth < 5) {
	open (LOG,">>$local_log");
	print LOG "** Warning Coverage below 5X. Pipeline stopped **\n";
	close LOG;

	#system ("rm $path/tmp*.*");
	#system ("rm $path/REFERENCE.*");
	#system ("rm $path/metrics.out");
	#system ("rm $path/flagstat.out");
	#system ("rm $path/indels.intervals");
	#system ("rm $path/reference*.*");
	#system ("rm $path/Stats.out");

} else {
	print "\n\n coverage passed\n\n";

	if (($depth >= 5) && ($depth < 40)) {
		open (LOG,">>$local_log");
		print LOG "\n** Warning Low Coverage **\n";
		close LOG;
	}

	#### lumpy_sv for long indel detections #####
	print "\n\n perl $lumpy_sv_script $path\n\n";
	system ("perl $lumpy_sv_script $path 2> /dev/null");


####### SNP calling Step ###

## create masked genome ##

open (FILE,"$reference");
my $header=<FILE>;
my $reference_sequence='';
while (my $in_reference=<FILE>) {
	chomp $in_reference;
	$reference_sequence .= $in_reference;  ##put the nucleotide sequence on a single line
}
close FILE;

open (FILE,$MTB_mask_region_file);
my $regions_to_mask=<FILE>;
while ($regions_to_mask=<FILE>) {
	chomp $regions_to_mask;
	my @split_locations = split ('\t',$regions_to_mask);
	my $start = $split_locations[0] - 1; ## Start of the block to mask
	my $size = $split_locations[1] - $split_locations[0]; ##size of the block to mask
	substr($reference_sequence,$start,$size+1) = "N" x ($size+1); ## replace block with N's
}
close FILE;

open (OUT,">$path/masked_genome.fasta");
print OUT "$header$reference_sequence";
close OUT;


print "\n\n samtools faidx $path/masked_genome.fasta \n\n";
system ("samtools faidx $path/masked_genome.fasta 2> /dev/null"); ## index reference
unless (-e "$path/masked_genome.dict") {  ## create dictionary
	system ("java -jar $picard CreateSequenceDictionary R=$path/masked_genome.fasta O=$path/masked_genome.dict 2> /dev/null");
}
my @job=();
my $pm = Parallel::ForkManager->new( 2 );
my $job_number;
$job[0]="$GATK -T UnifiedGenotyper -R $path/masked_genome.fasta -I $path/tmp_fixed_sorted_marked_RG_realigned.bam -o $path/GATK_SNP.vcf --output_mode EMIT_ALL_SITES -glm SNP -ploidy 2 -mbq 20 -et NO_ET -K /opt/galaxy/tools/nysdh_tools/mapping/pascal.lapierre_health.ny.gov.key 2> /dev/null";
$job[1]="$GATK -T UnifiedGenotyper -R $path/masked_genome.fasta -I $path/tmp_fixed_sorted_marked_RG_realigned.bam -o $path/GATK_INDELS.vcf --output_mode EMIT_VARIANTS_ONLY -glm INDEL -ploidy 2 -mbq 20 -et NO_ET -K /opt/galaxy/tools/nysdh_tools/mapping/pascal.lapierre_health.ny.gov.key 2> /dev/null";
### do SNP and INDELS in parallel
	for $job_number( 0..1 ) {
	$pm->start and next;
	system("java -jar $job[$job_number]");
	$pm->finish;
}

$pm->wait_all_children;

print "\n\n$annotator_script $path/GATK_INDELS.vcf $path/GATK_SNP.vcf $path $prefix >> $path/Reports/$prefix\_$date.txt 2> $path/log_gatk.out \n\n";
system ("$annotator_script $path/GATK_INDELS.vcf $path/GATK_SNP.vcf $path $prefix >> $path/Reports/$prefix\_$date.txt 2> $path/log_gatk.out");

system ("cp $path/GATK_INDELS.vcf $path/VCF/$prefix\_INDELS.vcf");
system ("cp $path/tmp.vcf $path/VCF/$prefix.vcf");
system ("gzip $path/VCF/$prefix\_INDELS.vcf --force");
system ("gzip $path/VCF/$prefix.vcf --force");



## get runtime ###
open (OUT,">$path/runtime.txt");

print OUT "WGS_ID\tIdentification\t\t$top_match\n\n";
print OUT "WGS_SPOLIGOTYPE\tWGS_SPOLIGO\t$spoligo\n\n";

my $t1 = Benchmark->new;
my $td = timediff($t1, $t0);
print OUT "\n\n-----------------------------------------------------\nTotal running time:" ,timestr($td),"\n";
close OUT;
system ("cat $path/runtime.txt >> $path/Reports/$prefix\_$date.txt");

#### cleanup ###

#system ("rm $path/tmp*.*");
#system ("rm $path/REFERENCE.*");
#system ("rm $path/masked*.*");
#system ("rm $path/metrics.out");
#system ("rm $path/flagstat.out");
#system ("rm $path/indels.intervals");
#system ("rm $path/reference*.*");
#system ("rm $path/runtime.txt");
#system ("rm $path/Stats.out");
#system ("rm $path/sample*.*");
#system ("rm $path/*.idx");
#system ("rm $path/*.vcf");
#
}
### output CLIMS ###

system ("cp $local_log $log_file");
system ("$clims_script $path/Reports/$prefix\_$date.txt $path");
system ("cat $path/$clims_filename $path/clims_temp.out >> $path/clims_temp2.out");
system ("mv $path/clims_temp2.out $path/$clims_filename");
system ("mv $path/$clims_filename $path/CLIMS/$clims_filename");
system ("mv $path/* $output_directory");
