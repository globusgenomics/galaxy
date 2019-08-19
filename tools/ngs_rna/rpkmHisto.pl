#!/usr/bin/env perl

use strict;
use warnings;

use English;
use Getopt::Long;
use IO::File;
use Pod::Usage;


my %Hash;

my $Total   = 0;
my $Genes   = 0;
my $Over0   = 0;
my $Over1   = 0;
my $Over5   = 0;
my $Over10  = 0;
my $Over25  = 0;
my $Over50  = 0;
my $Over100 = 0;

my $input_file;

my $old_cufflinks = 0;
my $help          = 0;

GetOptions(
           "i|input-file=s"  => \$input_file,
           "o|old-cufflinks" => \$old_cufflinks,
          ) or pod2usage(2);

pod2usage(1) if $help;
pod2usage(2) unless (defined($input_file));

my $rpkm_index = 9;
my $chr_index  = 6;

if ($old_cufflinks) { 
    $rpkm_index = 10;
    $chr_index  = 6;
}

my $fh = IO::File->new($input_file, 'r') || die "failed to open $input_file for reading: $OS_ERROR";

# discard header line
<$fh>;

LINE: 
while (my $line = <$fh>){

    chomp $line;

    my @cols = split /\t/, $line;

    $Genes++;
    
    my $rpkm = $cols[$rpkm_index];
	
    next LINE if ($rpkm==0); # only count genes that have a rpkm > 0;
	
    $rpkm = int $rpkm;
	
    $Hash{$rpkm}++;
    $Total++

}

my $Remaining = 100;

print "Genes: $Genes\n===============\n";

foreach my $x (sort {$a<=>$b} keys %Hash){

    my $percent = sprintf("%.1f",$Hash{$x}*100/$Genes);
    $Remaining  = sprintf ("%.1f", $Remaining - $percent);

    print "$x\t$Hash{$x}\t$Remaining\n";

#    if ($x ==0){$Zero += $Hash{$x}}
    if ($x >= 0){$Over0 += $Hash{$x}}
    if ($x > 1){$Over1 += $Hash{$x}}
    if ($x > 4){$Over5 += $Hash{$x}}
    if ($x > 9){$Over10 += $Hash{$x}}
    if ($x > 24){$Over25 += $Hash{$x}}
    if ($x > 49){$Over50 += $Hash{$x}}
    if ($x > 99){$Over100 += $Hash{$x}}

}
print "\n$Genes Genes\n";

print "$Over0 Genes detected\n";
print "$Over1 Genes expressed at 1 RPKM or greater\n";
print "$Over5 Genes expressed at 5 RPKM or greater\n";
print "$Over10 Genes expressed at 10 RPKM or greater\n";
print "$Over25 Genes expressed at 25 RPKM or greater\n";
print "$Over50 Genes expressed at 50 RPKM or greater\n";
print "$Over100 Genes expressed at 100 RPKM or greater\n\n";


__END__

=head1 NAME

rpmkHisto.pl - A script that produces a histogram from Cufflinks *.expr files

=head1 SYNOPSIS

rpmkHisto.pl --input-file genes.fpkm_tracking 

rpmkHisto.pl -i isoforms.fpkm_tracking 

rpmkHisto.pl --input-file transcripts.expr --old-cufflinks

=head1 OPTIONS

=over 8

=item B<-h>

Alais for -help.

=item B<-help>

Prints a brief help message and exits.

=item B<-o>

Alias for -old-cufflinks.

=item B<-old-cufflinks>

Specifies that input file format is from pre-1.0 cufflinks.

=item B<-i>

Alias for -input-file.

=item B<-input-file>

Specifies input file to process.

=back

=head1 DESCRIPTION

B<rpmkHisto.pl> processes the specified input file and outputs a histogram.

=cut
