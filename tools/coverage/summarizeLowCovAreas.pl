#!/usr/bin/perl

use strict;
use warnings;

my %Ex_Coords;
my %Hash;

my $usage = "\n\nfindLowCovAreas.pl zero_coverage_file exome_bed_file\n\n";

die $usage unless (@ARGV==2);
my $file = $ARGV[0];
open (IN2, $file) or die "Can't open $file\n"; ### test to see that the file opens before reading the exon file.

my $exon_file = $ARGV[1];

open (IN, $exon_file) or die;

my $Old = "old";
my $Count = 0;

while (<IN>){

    if (/chr/){
	my ($chr, $start, $end) = $_ =~ /(chr\S+)\s+(\d+)\s+(\d+)/;
	$Count = 0 if ($chr ne $Old);
	my $coords = "$start..$end";
#	print "$Count\t$chr\t$coords\n";
	push @{$Ex_Coords{$chr}}, $coords;
	$Hash{$chr}{$coords} = $Count;
	$Count++;
	$Old = $chr;
    }
}

#&test;


my $Last = "0..1";
my $String = "";
my $Length = 0;

print "COVERAGE GAPS\n\n";

while (<IN2>){

    if (/chr/){
	my ($chr, $start, $end) = $_ =~ /(chr\S+)\s+(\d+)\s+(\d+)\s+\d+\s+\d+$/;
#	print "$chr $start $end\n";
	my $coords   = "$start..$end";
	my $count    = $Hash{$chr}{$coords};
	my $previous = "x";
	$previous = ${$Ex_Coords{$chr}}[$count-1] unless ($count==1);
#	print "$Length\t$chr\t$coords\t$count\t$previous\n\t$_" if (! $count);
	if ($Last ne $previous || eof ){
	    if ($Length > 4){
		print "\n   Gap = $Length exons\n$_$String\n";
	    }
	    $Last = $coords;
	    $String = "";
	    $Length = 0;
	}
	else{
	    $Length++;
	    $String .= $coords." ";
	}
	$Last = $coords
    }
}


sub test{
    foreach my $x (sort keys %Ex_Coords){
	print "$x\n";
	foreach my $y (@{$Ex_Coords{$x}}){
#	    print "\t$y $Hash{$x}{$y}\n";
	}
    }
    exit
}
