#!/usr/bin/perl
##################################################
#
#	Author: shizhongbin
#
##################################################
use strict;
use warnings;

die "Usage: <fileName> \n" if(@ARGV < 1);

my $fileName = shift;
my $pairFile1;
my $pairFile2;
my $outputFile;

#my $i=0;
open (IN, $fileName) || die "fail open $fileName\n";
while (<IN>) {

	$pairFile1 = $_;
	$pairFile2 = <IN>;
	chomp $pairFile1;
	chomp $pairFile2;
	
	#$pairFile1 .= ".corr";
	#$pairFile2 .= ".corr";

	$pairFile1 =~ /(.*)_\w*1\.\w+\.corr$/;
	$outputFile = $1;
	
	print "merging $outputFile\n";
	system "perl merge_pair.pl $pairFile1 $pairFile2 $outputFile";
	
#	$i++;
#	if ($i%25==0) {
#		wait;
#	}
}
close IN;
