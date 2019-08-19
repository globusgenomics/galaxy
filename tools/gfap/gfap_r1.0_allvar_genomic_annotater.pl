#!/usr/bin/perl

use strict;
use warnings FATAL => qw[ numeric uninitialized ];
use File::Basename;
use Getopt::Long;

sub sepind{
		$_=shift @_ foreach my($str, $sep);
		my($pos, @pos);
		$pos=0;
		while(1){
			$pos=index($str, $sep, $pos);
			last if($pos<0);
			push @pos, $pos++;
		}
		return \@pos;
	}

my($varfile, $buildver, $refseq_dir, $cosmic_dir, $refseq_release, $cosmic_release, $annovar_release, $outdir, $noncoding, $coding, $cos, $ogs, $mid, $pid, $cno, $pno);
my(@buffer, @header, @ogs, @mid, @cno, @pno, @sep, %buffer, %mid, %ogs, %opts);

GetOptions(\%opts, "varfile=s", "buildver=s", "refseq_dir=s", "refseq_release=s", "cosmic_dir=s", "cosmic_release=s", "annovar_release=s", "outdir=s", "noncoding=s", "coding=s");
$varfile         = $opts{varfile};
$buildver        = $opts{buildver};
$refseq_dir      = $opts{refseq_dir};
$refseq_release  = $opts{refseq_release};
$cosmic_dir      = $opts{cosmic_dir};
$cosmic_release  = $opts{cosmic_release};
$annovar_release = $opts{annovar_release};
$outdir          = $opts{outdir};
$noncoding       = $opts{noncoding};
$coding          = $opts{coding};

my %legend=(
	'unk' => 'undefined column',
	'chr' => 'chromosome identifier',
	'start' => "${buildver} 1-based start position",
	'end' => "${buildver} 1-based end position",
	'ref' => 'reference allele',
	'alt' => 'alternate allele',
	'annot' => 'ig:intergenic; pp:1kb-upstream; 5|3u:UTR; in:intronic; ss:splice; nc:ncRNA',
	'ogs' => 'official gene symbol(s)',
	'cos' => "gene listed in cosmic ${cosmic_release} release",
	'mid' => "RefSeq mRNA identifier(s) from human.protein.gpff ${refseq_release} release",
	'pid' => "RefSeq protein identifier(s) from human.protein.gpff ${refseq_release} release",
	'c.x' => 'ATG-based variant descriptor in mRNA',
	'p.x' => 'ATG-based variant descriptor in protein'
);

my $annovar_src_dir = 'inc/annovar';
my $annovar_db_dir = "db/annovar/${annovar_release}";

my $fname = readlink($varfile) || $varfile;
$fname = basename($fname);

`${annovar_src_dir}/annotate_variation.pl -buildver $buildver ${outdir}/${fname} $annovar_db_dir 2> /dev/null` and die $!;

open IN, "<${refseq_dir}/mid2pid_${refseq_release}.txt" or die $!;
while(<IN>){
		next if /^#/;
		/^(\S+)\s+(\S+)/;
		$mid{$1}=$2;
	}
close IN;

open IN, "<${cosmic_dir}/${buildver}_cosmic_ogs_${cosmic_release}.txt" or die $!;
chomp and $ogs{$_}++ while(<IN>);
close IN;

open IN, "<${outdir}/${fname}" or die $!;
while(<IN>){
		last if $_!~/^#/;
		last if $_!~/=/;
		chomp;
		/^#(\S+)\s{1}=\s{1}(.+)/;
		push @header, $1;
		$legend{$1}=$2;
	}
if(!scalar(@header)){
		@header=('chr', 'start', 'end', 'ref', 'alt');
		$_=readline(IN);
		@_=split /\t/, $_;
		$_=$#_-4;
		push @header, ('unk')x$_ if($_!=0);
	}
close IN;
push @header, ('annot', 'ogs', 'cos');
open OUT, ">${outdir}/${fname}.nc" or die $!;
print OUT "#", join(' = ', $_, $legend{$_}), "\n" foreach @header;
print OUT "#", join("\t", @header), "\n";
open IN, "<${outdir}/${fname}.variant_function" or die $!;
while(<IN>){
		next if /exonic/;
		s/^downstream/ig/;
		s/;downstream//;
		s/,/:/g;
		s/(UTR(3|5))|(upstream)|(intronic)|(splicing)|(ncRNA)|(intergenic)/$1?"${2}u":$3?'pp':$4?'in':$5?'ss':$6?'nc':'ig'/eg;
		chomp;
		@buffer=split /\t/, $_;
		$buffer[1]='na' if $buffer[0] eq 'ig';
		$buffer[1]=~s/([^;]+);(?:\S+)$/$1/ if $buffer[0]!~/;/;
		print OUT join("\t", @buffer[2..$#buffer, 0..1], ($buffer[1] eq 'na')?'na':(exists $ogs{$buffer[1]})?'TRUE':'FALSE'), "\n";
	}
close IN;
close OUT;

$legend{annot}='fd:frameshift deletion; fi:frameshift insertion; nd:nonframeshift deletion; ni:nonframeshift insertion; bs:block substitution; ss:synonymous SNV; ns:nonsynonymous SNV; sg:stopgain; sl:stoploss; na:unknown';
push @header, ('mid', 'pid', 'c.x', 'p.x');

open IN, "${outdir}/${fname}.exonic_variant_function" or die $!;
open OUT, ">${outdir}/${fname}.cds" or die $!;
print OUT "#", join(' = ', $_, $legend{$_}), "\n" foreach @header;
print OUT "#", join("\t", @header), "\n";
while(<IN>){
		next if /unknown/;
		s/^\S+\s+//;
		chomp;
		%buffer=();
		@{$_}=() foreach (\@ogs, \@mid, \@cno, \@pno, \@sep);
		@buffer=split /\t/, $_;
		$buffer[0]=~s/(nonf\w+\s{1}(d|i|s)\w+)|(\w+\s{1}(d|i)\w+)|(stop(\w){1}.+)|(^(n|s).+)|(.+)/$1?(($2 eq 's')?'b':'n').$2:$3?"f$4":$5?"s$6":$7?"${8}s":'na'/eg;
		foreach (split /,/, $buffer[1]){
				@_=split /:/, $_;
				splice(@_, 2, 1);
				$_=shift(@_) || 'na' foreach ($ogs, $mid, $cno, $pno);
				$buffer{ogs}->{$ogs}->{$cno}->{$mid}++;
				$buffer{ono}->{$cno}=$pno;
			}
		$cos=0;
		foreach $ogs (@ogs=keys %{$buffer{ogs}}){
				push @cno, join('|', (@_=keys %{$buffer{ogs}->{$ogs}}));
				unshift @pno, $buffer{ono}->{$_} foreach reverse @_;
				$pno=join('|', @pno[0..$#_]);
				splice(@pno, 0, $#_+1);
				push @pno, $pno;
				foreach $cno (@_){
						push @mid, join(':', keys %{$buffer{ogs}->{$ogs}->{$cno}});
					}
				$cos++ if exists $ogs{$ogs};
			}
		$mid=join('|', @mid);
		$cno=join(';', @cno);
		if($#ogs!=0){
				(my $sep=$cno)=~s/[^;\|:]+//g;
				@_=@{sepind($mid, '|')}[@{sepind($sep, ';')}];
				substr($mid, $_, 1)=';' foreach @_;
			}
		($pid=$mid)=~s/([^;\|:]+)/$mid{$1} || 'na'/eg;
		push @buffer, shift @buffer, join(';', @ogs), ($cos!=0)?'TRUE':'FALSE', $mid, $pid, $cno, join(';', @pno);
		shift @buffer;
		print OUT join("\t", @buffer), "\n";
	}
close IN;
close OUT;
system "rm $noncoding $coding ${outdir}/${fname}*variant_function ${outdir}/${fname}*invalid*; ln -s ${outdir}/${fname}.nc $noncoding; ln -s ${outdir}/${fname}.cds $coding" and die $!;