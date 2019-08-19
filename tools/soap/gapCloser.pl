#! /usr/bin/perl -w
use strict;

die("parameter number wrong\n") unless (@ARGV == 7);

my ($scaffold, $library, $readlength, $overlap, $thread, $output, $stat)= @ARGV;

#my $cmd= "/usr/local/bgisoap/GapCloser -a $scaffold -b $library -o $output -l $readlength -p $overlap -t $thread";
my $cmd= "GapCloser -a $scaffold -b $library -o $output";

system $cmd;
if(! rename "$output.fill", "$stat"){
	warn "file exists, replace with the new file";
	unlink "$stat";
	rename "$output.fill", "$stat";
}

