#!/usr/bin/perl
##################################################
#
#	Author: shizhongbin
#
##################################################
use strict;
use warnings;

die "Usage: <pairFile1> <pairFile2> <outputFile>\n" if(@ARGV < 3);

my $pairFile1 = shift;
my $pairFile2 = shift;
my $outputFile = shift;

my $pairName1;
my $pairSeq1;
my $iPair1;
my $pairName2;
my $pairSeq2;
my $iPair2;
my $pairReadsum = 0;
my $singleReadsum = 0;
open PAIROUT,">${outputFile}.pair" or die "$!\n";
open SINGLEOUT,">${outputFile}.single" or die "$!\n";
open READSUMOUT,">${outputFile}.readsum" or die "$!\n";
open (PAIR1IN, $pairFile1) || die "fail open $pairFile1\n";
open (PAIR2IN, $pairFile2) || die "fail open $pairFile2\n";

if (!eof(PAIR1IN) && !eof(PAIR2IN)) {

	$pairName1 = <PAIR1IN>;
	$pairSeq1 = <PAIR1IN>;
	$pairName2 = <PAIR2IN>;
	$pairSeq2 = <PAIR2IN>;
	chomp $pairName1;
	chomp $pairSeq1;
	chomp $pairName2;
	chomp $pairSeq2;

	$pairName1 =~ /\[(\d*)\]/;
	$iPair1 = $1;
	$pairName2 =~ /\[(\d*)\]/;
	$iPair2 = $1;

	while (!eof(PAIR1IN) && !eof(PAIR2IN)) {
		
		if ($iPair1 < $iPair2) {
			
			$singleReadsum++;
			print SINGLEOUT "$pairName1\n$pairSeq1\n";
			$pairName1 = <PAIR1IN>;
			$pairSeq1 = <PAIR1IN>;
			chomp $pairName1;
			chomp $pairSeq1;
			$pairName1 =~ /\[(\d*)\]/;
			$iPair1 = $1;
		}
		elsif ($iPair1 > $iPair2) {

			$singleReadsum++;
			print SINGLEOUT "$pairName2\n$pairSeq2\n";
			$pairName2 = <PAIR2IN>;
			$pairSeq2 = <PAIR2IN>;
			chomp $pairName2;
			chomp $pairSeq2;
			$pairName2 =~ /\[(\d*)\]/;
			$iPair2 = $1;
		}
		else {
			$pairReadsum += 2;
			print PAIROUT "$pairName1\n$pairSeq1\n";
			$pairName1 = <PAIR1IN>;
			$pairSeq1 = <PAIR1IN>;
			chomp $pairName1;
			chomp $pairSeq1;
			$pairName1 =~ /\[(\d*)\]/;
			$iPair1 = $1;

			print PAIROUT "$pairName2\n$pairSeq2\n";
			$pairName2 = <PAIR2IN>;
			$pairSeq2 = <PAIR2IN>;
			chomp $pairName2;
			chomp $pairSeq2;
			$pairName2 =~ /\[(\d*)\]/;
			$iPair2 = $1;
		}
	}

	if ($iPair1 == $iPair2) {
		$pairReadsum += 2;
		print PAIROUT "$pairName1\n$pairSeq1\n";
		print PAIROUT "$pairName2\n$pairSeq2\n";
	}
	else {
		$singleReadsum += 2;
		print SINGLEOUT "$pairName1\n$pairSeq1\n";
		print SINGLEOUT "$pairName2\n$pairSeq2\n";
	}
}

while (!eof(PAIR1IN)) {
	
	$pairName1 = <PAIR1IN>;
	$pairSeq1 = <PAIR1IN>;
	chomp $pairName1;
	chomp $pairSeq1;
	$singleReadsum++;
	print SINGLEOUT "$pairName1\n$pairSeq1\n";
}

while (!eof(PAIR2IN)) {

	$pairName2 = <PAIR2IN>;
	$pairSeq2 = <PAIR2IN>;
	chomp $pairName2;
	chomp $pairSeq2;
	$singleReadsum++;
	print SINGLEOUT "$pairName2\n$pairSeq2\n";
}



print READSUMOUT "$outputFile\t$pairReadsum\t$singleReadsum\n"; 

close PAIROUT;
close SINGLEOUT;
close READSUMOUT;
close PAIR1IN;
close PAIR2IN;


