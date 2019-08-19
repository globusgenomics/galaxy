#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $INTERSECTBED = '/mnt/galaxyTools/tools/bedtools/2.16.2/intersectBed'; #*

my $usage = "\nScript annotates a list of indels in bed format using BedTool's intersect bed.  The output of intersect bed is printed to a file.  The output is also formatted and printed to a tab-separated file.

Usage: $0 -i <indels in exon> -b <bedintersect output> -a <annotated indels> -e <exons>";

my ($help, $indels_in_exon, $bedintersect_output, $indels_in_exon_annotated, $EXONS);

GetOptions(
          "h|help|?"                 => \$help,
          "i|input-indels=s"         => \$indels_in_exon,
          "b|bedintersect-output=s"  => \$bedintersect_output,
          "a|annotated-indels=s"     => \$indels_in_exon_annotated,
          "e|exons=s"                => \$EXONS,
          );

if ($help) {

    die $usage;
}

unless ((defined $bedintersect_output) &&
       (defined $indels_in_exon_annotated) &&
       (defined $indels_in_exon_annotated) &&
       (defined $EXONS)) {

    die "Incorrect number of aruments.\n\n$usage";
}

unless (-e $indels_in_exon) {
    die "Error: indels in exon file '$indels_in_exon' does not exist\n$usage";
}

unless (-s $indels_in_exon) {

    warn "Warning: indels in exon file '$indels_in_exon' is empty.  Script will assume that there are no indels\n";
}

my $command = "$INTERSECTBED -a $indels_in_exon -b $EXONS -wa -wb > $bedintersect_output";
system($command) == 0 or die "Comparison to CCDS_exon file Failed\n\n";


open (IN, $bedintersect_output) or die "Can't open semi annotated indels in exon file '$bedintersect_output': $!\n\n";

open (OUT, '>', $indels_in_exon_annotated) or die "Can't open fully annotated indels in exon file '$indels_in_exon_annotated': $!\n\n";

print OUT "Chromosome\tPosition\tIndel\tAlleles\tCon_qual\tSNP_qual\tMap_qual\tRead_depth\tGene\tGene_id\tCCDS_id\n";

my $Last_pos = 0;
my @Output;
    
while (<IN>) {

    chomp $_;

    my @columns = split (/\t/, $_);
    my $pos       = $columns[2];
    my $ccds      = $columns[15];

    if ($pos == $Last_pos){ # redundant information
	$Output[10] .= ",$ccds";
    }	      
    else {
	if (@Output){ # first line doesn't have a value	    
	    &print_out;
	}
	
	@Output = ("$columns[0]", "$columns[1]", "$columns[3]", "$columns[4]", "$columns[5]", "$columns[6]", "$columns[7]", "$columns[8]", "$columns[12]", "$columns[13]", "$columns[15]");
	
    }
    $Last_pos = $pos;
    
}

&print_out; # print the last line;

close(IN);
close(OUT);

sub print_out {
    foreach my $x (@Output){
        print OUT "$x\t";
    }
    print OUT "\n";
}


