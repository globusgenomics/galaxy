#! /usr/bin/perl -w
use strict;
use Getopt::Std;
#use Cwd;
#my $path= cwd;
#print "$path\n";
my %opt;
getopts('a:b:o:O:Z:2:3:l:u:s:m:k:t:f:L:M:H:x:',\%opt);    #x: additional para, O: cross_connect.fq, Z:log file
$opt{'2'} ||= "fail.$$.read1.fq";
$opt{'3'} ||= "fail.$$.read2.fq";
$opt{'x'} ||= ' ';

my $xpara= $opt{'x'};
my $cross= $opt{'O'};
my $log= $opt{'Z'};
delete($opt{'x'});
delete($opt{'O'});
delete($opt{'Z'});

my @para;
foreach my $key (sort keys %opt){
	push @para, "-$key $opt{$key}";
}

my $command= "cope ".join (" " ,@para)." ".$xpara." >cope.$$.log 2>cope.$$.err";
print "$command\n";

if($opt{'m'} >0){
	die("-k -t -f are needed in k-mer assisted mode\n") unless( exists $opt{'k'} && exists $opt{'t'} && $opt{'f'});
}
mkdir "cope$$", 0755;
chdir "cope$$";
print "==running==\n";
system "$command";
rename "cope.$$.log", "$log";
rename "cross_connect.fq", "$cross" if(-e "cross_connect.fq");
unlink glob "*";
chdir "../";
rmdir "cope$$";
print "==out==\n";
print "$opt{'o'}\n$log\n";
print "$cross\n" if(defined $cross);
print "==done==\n";
