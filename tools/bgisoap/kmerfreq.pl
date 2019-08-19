#! /usr/bin/perl -w
use strict;
use Getopt::Std;

my %opt;
getopts('k:c:q:m:b:t:p:1:2:3:4:',\%opt);
$opt{"p"} ||= "kmerfreq.$$";

my ($out1,$out2,$out3,$out4)= ($opt{'1'}, $opt{'2'}, $opt{'3'}, $opt{'4'});
delete ($opt{'1'});
delete ($opt{'2'});
delete ($opt{'3'});
delete ($opt{'4'});

my @para='';
foreach my $key (sort keys %opt){
	push @para, "-$key $opt{$key}";
}

open(FQLIST,">kmerfreq.$$.list");
foreach my $file (@ARGV){
	print FQLIST "$file\n";
}

my $command= "kmerfreq ".join (" ", @para)." "."kmerfreq.$$.list".">kmerfreq.$$.cout 2>kmerfreq.$$.cerr";
print "==command==\n";
print "$command\n";
print "==running==\n";
system "$command";
print "==output==\n";
rename "kmerfreq.$$.freq.cz", $out1;
rename "kmerfreq.$$.genome_estimate", $out2;
rename "kmerfreq.$$.freq.stat", $out3;
rename "kmerfreq.$$.freq.cz.len", $out4;
print "Out:$out1,$out2,$out3,$out4\n==Done==\n";
unlink glob ("kmerfreq.$$.*");
