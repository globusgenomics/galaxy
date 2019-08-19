#! /usr/bin/perl -w
use strict;

die("Usage: perl $0 <mode> <target sequence> <in format> <out format> <coverage> <out detail> <map file1> .. <mapfileN>\n") if( @ARGV<7);

my $mode= shift @ARGV;
my $target= shift @ARGV;
my $format= shift @ARGV;
my $outtype= shift @ARGV;
my $coverage= shift @ARGV;
my $outdetail= shift @ARGV;
my $maplist= join " ", @ARGV;


my $cmd;
if($outtype eq "fasta"){
	$cmd= "soap.coverage -$mode -refsingle $target -i $maplist -o $coverage -depthsingle $outdetail";
	$cmd= "$cmd -sam" if($format eq "sam");
	system "$cmd";
}elsif($outtype eq "wig"){
	$cmd= "soap.coverage -$mode -refsingle $target -i $maplist -o $coverage -depthsingle $outdetail.tmp";
	$cmd= "$cmd -sam" if($format eq "sam");
	system "$cmd";
	convert($outdetail);
	unlink "$outdetail.tmp";
}

#####################################
sub convert{
	my $cor=1;
	open(WIG,">$_[0]");
	open (TMP, "$_[0].tmp") || die "fail to open the file";
	while(<TMP>){
		chomp;
		if(/^>/){
			/^>(\S+)/;
			print WIG "variableStep chrom=$1\n";
			$cor=1;
		}else{
			my @seq= split / /;
			foreach my $base (@seq){
				if(defined $base){
					print WIG "$cor $base\n";
				}
				$cor++;
			}
		}
	}
	close TMP;
	close WIG;
}
