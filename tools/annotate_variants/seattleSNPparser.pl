#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;


my $usage = "\nParses the relevant information from a SeattleSeq SNP file\n\n\tseattleSNPparser.pl --output-from=seattle-seq <SeattleSeq output file> --input-to-seattle-seq <SeattleSeq input file> --annotated-snps-output <annotated snps output file>\n\n";


my ($help, $seattle_seq_output, $seattle_seq_input, $annotated_snps);
GetOptions(
          'h|help|?'                    => \$help,
          'so|output-from-seattle-seq=s' => \$seattle_seq_output,
          'si|input-to-seattle-seq=s'  => \$seattle_seq_input,
          'a|annotated-snps-output=s'   => \$annotated_snps,
          );

if ($help) {

    die $usage;
}

unless ((defined $seattle_seq_output) &&
       (defined $seattle_seq_input) &&
       (defined $annotated_snps)) {

    die "Incorrect number of arguments.\n\n$usage";
}

my %Hash;
&get_scores;

open (IN, $seattle_seq_output) or die "Can't open $seattle_seq_output: $!\n\n";
open (OUT, '>', $annotated_snps) or die "Can't open $annotated_snps: $!\n\n";

my $Last_pos = "033";
my $Last_func= "diddly";
my @Output;

print OUT "Chromosome\tPosition\tRef_nt\tAlleles\tConsensus_qual\tSNP_qual\tMap_qual\tRead_depth\tGene\tNCBI_id\tCCDS_id\tGVS_func\tAminoAcid\tAA_position\tPolyPhen\tDatabase\n";

while (<IN>){
    
    if ($_ !~ /^#/){ #skip their first header line
	my @clmns     = split(/\t/,$_);
	my $chr       = $clmns[1];
	my $pos       = $clmns[2];
	my $acc       = $clmns[7]; #combine multiple accessions to print
	my $func      = $clmns[8];
	$func =~ s/coding-synonymous/cod-syn/;
#	print "\t### $pos $Last_pos $acc $func $Last_func\n";
	if ($pos == $Last_pos){ # possibly redundant information
	    # see if there is any new information
	    if ($Last_func !~/$func/){ # print both lines of data, because of unique function
		# append new data to the output;
		$Output[7] .= ",$func";
	    }

	    else { #don't print both lines of data, just add acc#
		if ($acc =~/NM/){
		    $Output[6] .= ",$acc";
		}
		else {
		    if ($Output[7] eq "ccds"){
			$Output[7] = $acc;
		    }
		    else {$Output[7] .= ",$acc";
		    }
		}
	    }
	}
	else {
	    if (@Output){ # first line doesn't have a value
	       
		&print_out;
	    }
# in seattleSNP131 column 19 NickLab score is deleted
#	    must change $clmns[19] => $clmns[18] to grab the gene name
	    
	    @Output = ("chr$clmns[1]", $clmns[2], $clmns[3], $clmns[5], "$Hash{$chr}{$pos}", $clmns[18], $clmns[7], "ccds", $clmns[8], $clmns[11], $clmns[12], $clmns[13], $clmns[0]);
	}
	$Last_pos = $pos;
	$Last_func = $func;
	
    }
#    elsif ($_ !~ /#\s*in/){ print $_}
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

sub get_scores{

    open (IN1,$seattle_seq_input) || die "Could not open $seattle_seq_input for getting scores: $!\n\n";

    while (<IN1>){
	if (/chr(\S+)\s+(\S+)\s+\S+\s+\S+\s+(\d+\s+\d+\s\d+\s+\d+)/){
	    my $chr    = $1;
	    my $pos    = $2;
	    my $scores = $3;

	    $Hash{$chr}{$pos}=$scores;
	}
    }
    close(IN1);
}
