#!/usr/bin/perl
use strict;
use Getopt::Long;
use Cwd;
use Cwd 'abs_path';
use Bio::SearchIO; 

my ($input, $qty, $help, $evalue, $identity, $output);
GetOptions('input|i=s' => \$input, 'help|h' => \$help, 'evalue|e=s' => \$evalue, 'percentage|p=s' => \$identity, 'output|o=s' => \$output, 'number|n=s' => \$qty
    );
my $usage = "filter_fasta [--evalue | --identity] < --output output_file > <--number qty of top blast hits for a query>  < --input filename >\n";

if ($help)
{
    die $usage;
}

if (! $qty)
{
    $qty = 1;
}

# without bioperl
open (FH, $input) || die "could not find file $input\n\nNeed input as:\n$usage\n";
open (FO, ">$output") || die "need output file name\n\n$usage\n";

my $seen={};
while (my $line = <FH>)
{
    chomp $line;
    my ($query,$hit,$hidentity,$hlength,$qgap,$hgap,$qstart,$hstart,$qend,$hend,$hevalue,$hscore) = split (/\t/, $line);
    next if ($seen->{$query} && $seen->{$query} >= $qty);
    $seen->{$query}++;
    
    if (($identity && $hidentity >= $identity) ||
        ($evalue && $hevalue <= $evalue))
    {
	print FO $line . "\n";
    }
}
close FH;
close FO;
