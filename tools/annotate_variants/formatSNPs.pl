#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $usage = "\nThis script formats BED SNP files for submission to SeattleSEQ\n\n\tformatSNPs.pl (snpInExon.txt in dir)\n\nIt reads the snpInExome.txt file and creates files:\nsnpInExome_SSformat.txt AND indelsInExome.txt\n\n";

my ($help, $snps_in_exon, $snps_in_exon_ss_format, $indels_in_exon);

GetOptions(
          "h|help|?"         => \$help,
          "i|input-file=s"   => \$snps_in_exon,
          "os|output-snps=s"   => \$snps_in_exon_ss_format,
          "oi|output-indels=s" => \$indels_in_exon,
          ); 

if ($help) {

	die $usage;
}

unless ((defined $snps_in_exon) &&
       (defined $snps_in_exon_ss_format) &&
       (defined $indels_in_exon)) {

	die "Incorrect number of arguments given.  $usage";
}

open (IN,$snps_in_exon) or die $usage;
open (OUT1,">",$snps_in_exon_ss_format) or die "Can't open output file\n";
open (OUT2,">",$indels_in_exon) or die "Can't open output file\n";


while (<IN>){
    
    my @columns = split(/\t/,$_);

    if ($columns[3] eq "*"){
    print OUT2 "$columns[0]\t$columns[1]\t$columns[2]\t$columns[3]\t$columns[4]\t$columns[5]\t$columns[6]\t$columns[7]\t$columns[8]\n";
    }
    else{
    # chr  pos  refbase  genotype  consensus_qual  SNP_qual  map_qual  coverage
	print OUT1 "$columns[0]\t$columns[2]\t$columns[3]\t$columns[4]\t$columns[5]\t$columns[6]\t$columns[7]\t$columns[8]\n";
    } 
}

close(IN);
close(OUT1);
close(OUT2);
