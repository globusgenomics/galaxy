#!/usr/bin/perl

use strict;
#use lib 'inc/perlmod';
#use ngsutil qw[ :DEFAULT &varscan ];
use warnings FATAL => qw[ numeric uninitialized ];
use File::Basename;
use Getopt::Long;

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#	TEMP include ngsutil.pm
sub explode_varcall{
		my $N=0;
		$_=shift @_ foreach my($POS, $REF, $ALT);
		$_=$POS foreach my($START, $END);
		my(@length, @range, @idx, @VAR, @POS);
		@{$_}=() foreach (\@length, \@range, \@idx, \@VAR, \@POS);
		push @length, length($_) foreach ($REF, $ALT);
		@range=sort{ $a<=>$b } @length;
		if($range[0]==1){
			if($range[1]!=1){
				foreach ($REF, $ALT){
						$_=substr($_, 1);
						$_=~s/^$/-/;
					}
				if($length[0]!=1){
						$END+=$length[0]-1;
						$START++;
					}
			}
			push @POS, $START, $END;
			push @VAR, $REF, $ALT;
		}else{
			my @N=();
			undef $_ foreach my ($i, $VAR);
			$_-=2 foreach (@length, @range);
			$_++ foreach ($START, $END);
			$_=substr($_, 1) foreach ($REF, $ALT);
			my $indel='-' x ($range[1]-$range[0]);
			$VAR.=($_>$range[0])?
				('-'):((substr($REF, $_, 1) ne substr($ALT, $_, 1))?
					0:1) for 0 .. $range[1];
			$N++ while $VAR =~ /0/g;
			if($length[0]<$length[1]){
				@VAR=($VAR);
				@N=($N);
				$N=0;
				undef($VAR);
				$VAR.=($_>$range[0])?
					('-'):((substr($REF, $length[0]-$_, 1) ne substr($ALT, $length[1]-$_, 1))?
						0:1) for reverse 0 .. $range[1];
				$N++ while $VAR =~ /0/g;
				if($N>=$N[0]){ $N=shift(@N); $VAR=shift(@VAR); }
				else{ $REF=$indel . $REF; }
			}else{ $ALT.=$indel; }
			foreach (qw[ 0 \- ]){
					push @idx, [ $-[0], $+[0]-$-[0] ] while ($VAR =~ /$_+/g);
				}
			@{$_}=() foreach (\@VAR, \@POS);
			foreach my $k (@idx){
					push @VAR, substr($_, ${$k}[0], ${$k}[1]) || '-' foreach ($REF, $ALT);
					push @POS, ${$k}[0], sum(@{$k})-1;
				}
			$_+=$START foreach @POS;
			$_=~s/\-+/\-/ foreach @VAR;
			for($i=0; $i<$#POS; $i+=2){ $POS[$i+1]=$POS[$i] if $VAR[$i] eq '-'; }
		}
		return(\@POS, \@VAR);
	}

sub varscan{
		$_=shift @_ foreach my($kname, $fpath, $href);
		my($k, @buffer);
		open IN, "<$fpath" or die $!;
		while(<IN>){
				next if /^#/;
				chomp;
				@buffer=split /\s+/, $_;
				next if !exists $$href{($k=join(':', @buffer[0..2]))};
				next if $$href{$k}->{ref} !~ $buffer[3];
				next if $$href{$k}->{alt} !~ $buffer[4];
				splice(@buffer, 0, 5);
				$$href{$k}->{$kname}=join(':', @buffer);
			}
		close IN;
	}
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

my($varfile, $buildver, $outdir, $dir_1000g, $dir_dbsnp, $dir_cosmic, $release_1000g, $release_dbsnp, $release_cosmic, $outfile, $k, @buffer, @varlist, %opts, %varlist);

GetOptions(\%opts, "varfile=s", "buildver=s", "outdir=s", "dir_1000g=s", "dir_dbsnp=s", "dir_cosmic=s", "release_1000g=s", "release_dbsnp=s", "release_cosmic=s", "outfile=s");
$varfile        = $opts{varfile};
$buildver       = $opts{buildver};
$outdir         = $opts{outdir};
$dir_1000g      = $opts{dir_1000g};
$dir_dbsnp      = $opts{dir_dbsnp};
$dir_cosmic     = $opts{dir_cosmic};
$release_1000g  = $opts{release_1000g};
$release_dbsnp  = $opts{release_dbsnp};
$release_cosmic = $opts{release_cosmic};
$outfile        = $opts{outfile};

my $fname = readlink($varfile) || $varfile;
$fname = basename($fname);

my %k=(
	'1000g' => {
		'dir' => $dir_1000g, 'release' => $release_1000g, 'value' => join(':', ('0.00000')x5), 'header' => join(':', 'AF_ALL', 'AF_AFR', 'AF_AMR', 'AF_ASN', 'AF_EUR')
	}, 'dbsnp' => {
		'dir' => $dir_dbsnp, 'release' => $release_dbsnp, 'value' => join(':', ('na')x2), 'header' => join(':', 'rs', 'dbsnp')
	}, 'cosmic_var' => {
		'dir' => $dir_cosmic, 'release' => $release_cosmic, 'value' => join(':', '0.00000', 'na'), 'header' => join(':', 'AF_COS', 'cid')
	}
);

my %legend=(
	'chr' => 'chromosome identifier',
	'start' => "${buildver} 1-based start position",
	'end' => "${buildver} 1-based end position",
	'ref' => 'reference allele',
	'alt' => 'alternate allele',
	'QC' => 'Phred-scaled call quality',
	'NRF' => '#reads consistent w/ the reference allele on the F-strand',
	'NRR' => '#reads consistent w/ the reference allele on the R-strand',
	'NAF' => '#reads consistent w/ the alternate allele on the F-strand',
	'NAR' => '#reads consistent w/ the alternate allele on the R-strand',
	'DP' => 'total #reads in call ie. NRF+NRR+NAF+NAR',
	'AD' => 'total #reads consistent w/ the alternate allele ie. NAF+NAR',
	'AF' => 'alternate allele ratio ie. AD/DP',
	'VCF.FILTER' => 'FILTER field from the input vcf file',
	'DPT.FILTER' => 'check for heterogeneous depth in substituted blocks',
	'VAR.FILTER' => 'GFAP default FILTER to discriminate between TP and FP variants',
	'P.str' => 'NRF+NAF vs. NRR+NAR binomial test P-value ie. total strand bias',
	'P.ref' => 'NRF vs. NRR binomial test P-value ie. reference allele strand bias',
	'P.alt' => 'NAF vs. NAR binomial test P-value ie. alternate allele strand bias',
	'AF_ALL' => "global AF in ${release_1000g} 1000g data",
	'AF_AFR' => "AF in AFR ${release_1000g} 1000g data",
	'AF_AMR' => "AF in AMR ${release_1000g} 1000g data",
	'AF_ASN' => "AF in ASN ${release_1000g} 1000g data",
	'AF_EUR' => "AF in EUR ${release_1000g} 1000g data",
	'AF_COS' => "AF in ${release_cosmic} cosmic data",
	'rs' => "dbsnp rs identifier(s) from ${release_dbsnp} release",
	'dbsnp' => "dbsnp build version(s) from ${release_dbsnp} release",
	'cid' => "cosmic mutation identifier from ${release_cosmic} release"
);
my @header=('chr', 'start', 'end', 'ref', 'alt', 'DPT.FILTER', 'QC', 'NRF', 'NRR', 'NAF', 'NAR', 'VCF.FILTER', 'P.str', 'P.ref', 'P.alt', 'DP', 'AD', 'AF', 'VAR.FILTER');
my @k=qw[ 1000g dbsnp cosmic_var ];

open IN, "<$varfile" or die $!;
while(<IN>){
		chomp;
		@buffer=split /\s+/, $_;
		$buffer[0]=~s/^chr(.+)$/$1/;
		push @varlist, ($k=join(':', @buffer[0..2]));
		shift(@buffer) for 0..2;
		$varlist{$k}->{$_}=shift(@buffer) foreach qw[ ref alt ];
		$varlist{$k}->{cov}=join(':', (($buffer[0] eq 'unk')?'SKIP':'PASS'), @buffer[1..$#buffer]);
	}
close IN;

foreach $k (@k){
		push @header, split(/:/, $k{$k}->{header});
		varscan($k, $k{$k}->{file}, \%varlist);
	}

my @idx=(0..4,7..10,15..17,6,12..14,11,5,18..23,26..27,24..25);
open OUT, ">${outdir}/${fname}.dbi" or die $!;
print OUT '#', join(' = ', $_, $legend{$_}), "\n" foreach @header[@idx];
print OUT '#', join("\t", @header[@idx]), "\n";
foreach $k (@varlist){
		@buffer=(split(/:/, 'chr'.$k), $varlist{$k}->{ref}, $varlist{$k}->{alt});
		push @buffer, split(/:/, ($varlist{$k}->{$_} || $k{$_}->{value})) foreach ('cov', @k);
		print OUT join("\t", @buffer[@idx]), "\n";
	}
close OUT;

system "rm $outfile; ln -s ${outdir}/${fname}.dbi $outfile" and die $!;