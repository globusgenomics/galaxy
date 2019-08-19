#!/usr/bin/perl

use strict;
use warnings;

my $usage = "coverageStats.pl [ bed file ]
\n\tthis script calculates the coverage statistics from a bed file\n\n";

die $usage unless (@ARGV==1);

my $file = $ARGV[0];
my $test_line = `head -1 $file`;

my %Coverage;
my $Nts = 0;

die "$file not in correct format" unless $test_line =~ /^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)/;

open (IN, $file) or die;

LINE:while (<IN>){
    my ($chr, $start, $end, $cov) = $_ =~/(\S+)\s+(\d+)\s+(\d+)\s+(\d+)/;

    for (my $i=$start+1; $i<=$end; $i++){
	$Coverage{$cov}++;
	$Nts++;
    }
}


print "$Nts nucleotides in $file\n";

my $Remaining = 100;
foreach my $x (sort {$a<=>$b} keys %Coverage){
    my $percent = sprintf("%.3f",$Coverage{$x}*100/$Nts);
    print "$x\t$Coverage{$x}\t$percent\t$Remaining\n";
    $Remaining = sprintf ("%.3f", $Remaining -= $percent);
}
