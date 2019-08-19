#!/usr/bin/perl

use strict;

## changelog v0.2 ####
## remove reporting of mutations
## for pipeline status : if presence of other species than MTB, Mbovis, Status fail.  Low coverage <5x, Status fail

my $path = $ARGV[1];
my $report_to_parse = $ARGV[0];


my (
	## spoligotype
	$spoligo,
	
	## species ID
	$species,
	
	## reads name
	$read_R1,
	$read_R2,
	
	## Percent mapped reads
	$percent_mapped,
	
	## genome coverage
	$cov,
	
	## average reads depth
	$dp,
	
	## Array to store mutations
	@mutations,
	
	## final status of the pipeline, 'Passed' as default value
	$status,
	
	## Array containing the indels list if any
	@indels,
	
	## value to tell if sample is an MTB or Negative control. Default valus set at 0 for false.
	$mtb,
	
	$infile
);
$status="Passed";
$mtb=0;


open (FILE,$report_to_parse);
do {  ##skip header and other unwanted lines
	$infile=<FILE>;
}
until ($infile =~ 'Pipeline script');
chomp $infile;

## get sample name and rundate from line
my @split_report_line = split ('Sample ',$infile);
my @split_sample_name = split (', rundate ',$split_report_line[1]);
my $sample_name= "$split_sample_name[0]";
my @split_for_rundate = split ('\_',$split_sample_name[1]);
my $rundate= "$split_for_rundate[1]\/$split_for_rundate[2]\/$split_for_rundate[0]";


## initialize value of % match kraken to 0
my $kraken_percent=0;

while ($infile = <FILE>) {
chomp $infile;
## check for MTB complex
	if ($infile =~ 'Mycobacterium tuberculosis complex') {
		my @split_kraken_line = split ('\t',$infile);
		$split_kraken_line[1] =~ s/\%//g;
		$kraken_percent=$split_kraken_line[1];
	}


## Get Spoligotyping
	if (($infile =~ 'Spolygotyping') && ($infile =~ 'FAIL')) {
		$spoligo="Failed";
	} elsif (($infile =~ 'Spolygotyping') && ($infile =~ 'PASS')) {
		$infile=<FILE>;$infile=<FILE>;
		chomp $infile;
		
		## initialize $sit variable to not determined.
		my $sit="SIT# ND";
		if ($infile !~ 'Unknown type') {
			my @split_sit = split ('\t',$infile);
			$sit=$split_sit[3];
		}
		
		$infile = <FILE>;
		chomp $infile;
		my @split_spoligo_part = split ('\t',$infile);
		if ($infile =~ 'Unknown type') {
			if ($sit) {
				$spoligo = "$split_spoligo_part[3]($split_spoligo_part[2])";
			} else {
				$spoligo = "$split_spoligo_part[3]($split_spoligo_part[2])";
			}
		} else {
			$spoligo = "$split_spoligo_part[3]($split_spoligo_part[4])SIT# $sit";	
		}
	}
	
	
## Get Species ID
	if (($infile =~ 'Mycobacterium Species detected|No species identified')) {
	
		chomp $infile;
		my @split_infile = split ('\: ',$infile);
		my @split_species_id = split ('\,',$split_infile[1]);
		if (($split_species_id[0] =~ 'Mycobacterium tuberculosis') && ($kraken_percent >=90)) {
			$species="Mycobacterium tuberculosis";	
			$mtb=1;
		} elsif ($kraken_percent >=0) {
			$species="No species identified";	
			$mtb=1;
		}
	}

## Get read locations
	if (($infile =~ 'Map_reads.pl')) {
		chomp $infile;
		my @split_line= split (' ',$infile);
		$read_R1=$split_line[1];
		$read_R2=$split_line[2];
	}
## Get mapping stats
	if (($infile =~ 'Percent of mapped')) {
		chomp $infile;
		my @split_line = split ('\t',$infile);
		$percent_mapped=sprintf "%.2f", $split_line[1];
		$percent_mapped =~ s/\%//g;
	}
		
	if (($infile =~ 'Genome Coverage')) {
		chomp $infile;
		my @split_line = split ('\t',$infile);
		$cov=$split_line[1]*100;
		$cov  = sprintf "%.2f", $cov;
	}
	if (($infile =~ 'Average Depth')) {
		chomp $infile;
		my @split_line = split ('\t',$infile);
		$dp=$split_line[1];
		$dp  = sprintf "%.2f", $dp;
		if ($dp < 5) {
			$status ="Failed";
		}
	}
		
## Get run status

	if (($infile =~ 'Resistance Report')) {		
		do {
		$infile=<FILE>;
		if ($infile =~ "FAILED") {
			$status="Failed";
			}		
		}
		until ($infile =~ '-------');
	}

}	
close FILE;




if ($sample_name =~ '\_NT\_') {
	$mtb==1;
}


## print final report to file
if ($mtb==1) {
	open (OUT,">$path/clims_temp.out");
	print OUT "$sample_name\t$rundate\t$percent_mapped\t$cov\t$dp\t\t\t$status\t$spoligo\t\tWhole Genome Sequencing by Illumina MiSeq of MTB complex\tShotgun MTB complex\tMiSeq deep shotgun sequencing of cultured isolate\tWGS\tGENOMIC\tRANDOM\tPAIRED\tILLUMINA\tIllumina MiSeq\tfastq\t$read_R1\t$read_R2\t\t\t\t\n";
} 

close OUT;

