#!/usr/bin/perl

use strict;
use warnings FATAL => qw[ numeric uninitialized ];
use List::Util qw[ sum min max ];
use List::MoreUtils qw[ first_index ];
use File::Basename;
use Getopt::Long;

my($varfile, $buildver, $dbdir, $release, $outdir, $outfile, $max, $i, $k);
my(@buffer, @legend, @header, @k, @score, @tools, @Temp, %buffer, %AAS, %dbscore, %dbtools, %opts);

GetOptions(\%opts, "varfile=s", "buildver=s", "dbdir=s", "release=s", "outdir=s", "outfile=s");
$varfile  = $opts{varfile};
$buildver = $opts{buildver};
$dbdir    = $opts{dbdir};
$release  = $opts{release};
$outdir   = $opts{outdir};
$outfile  = $opts{outfile};

my $fname = readlink($varfile) || $varfile;
my $dbfile="${dbdir}/${buildver}_dbnsfp_${release}.txt";
$fname = basename($fname);

open IN, "<$varfile" or die $!;
open OUT, ">${outdir}/${fname}.Temp" or die $!;
while(<IN>){
		push @legend, $1 and next if /^#(.+=.+)/;
		next if $_!~/\b(?:s(?:g|l))|ns\b/;
		next if /\s+-\s+/;
		/^(?:chr)*(\S+)\s+(\S+)/;
		@{$buffer{($k=join('_', $1, $2))}->{dbnsfp}}=();
		push @k, $k;
		print OUT $_;
	}
close IN;
close OUT;

$i=first_index{ /^annot/ } @legend;
@_=$legend[$i]=~/((?:s(?:g|l)|ns):[^;|\s]+)/g;
$legend[$i]=join(' = ', 'annot', join('; ', @_));
push @legend, (
	'AAS = Amino Acid Substitution(s)',
	"FIS = Functional Impact Score(s) from dbnsfp ${release} release",
	"OCC = number of tools from which FIS was/were calculated",
	"FIS.max = highest score among FIS",
	"OCC.max = number of tools from which FIS.max was calculated",
	"PRED = qualitative ternary classifier ie. [L]ow; [M]edium; [H]igh"
);
foreach (@legend){
		/^(\S+)/;
		push @header, $1;
	}

open IN, "<$dbfile" or die $!;
while(<IN>){
		next if /^#/;
		/^(\S+)\s+(\S+)(?:\s+\S+){2}\s+(.+)/;
		next if !exists $buffer{($k=join('_', $1, $2))};
		push @{$buffer{$k}->{dbnsfp}}, join(':', split /\t/, $3);
	}
close IN;
open IN, "<${outdir}/${fname}.Temp" or die $!;
open OUT, ">${outdir}/${fname}.dbnsfp" or die $!;
print OUT "#", $_, "\n" foreach @legend;
print OUT "#", join("\t", @header), "\n";
foreach $k (@k){
		$i=0;
		$_=readline(IN);
		chomp;
		@buffer=split /\s+/, $_;
		%{$_}=() foreach (\%AAS, \%dbscore, \%dbtools);
		foreach (split(/[;\|]/, $buffer[-1])){
				$AAS{$1.$2}++ if /^p\.(\w{1})\d+(\w{1})$/;
			}
		if($#{$buffer{$k}->{dbnsfp}}<0){
			unshift @buffer, (%AAS)?(join(':', keys %AAS)):('na'), (join(':', ('na') x max(scalar(keys %AAS), 1))) x 2;
		}elsif(%AAS){
			foreach (@{$buffer{$k}->{dbnsfp}}){
					@Temp=split /:/, $_;
					$k=shift @Temp;
					@{$_}=split(/;/, pop @Temp) foreach (\@tools, \@score);
					foreach (split /;/, shift @Temp){
							$dbscore{$k.$_}=shift @score;
							$dbtools{$k.$_}=shift @tools;
						}
				}
			foreach (keys %AAS){
				push @score, $dbscore{$_} || 'na';
				push @tools, $dbtools{$_} || 'na';
			}
			unshift @buffer, join(':', keys %AAS), join(':', @score), join(':', @tools);
		}
		push @buffer, shift @buffer for 1..3;
		@{$_}=grep{ !/na/ } split(/:/, $buffer[--$i]) foreach (\@tools, \@score);
		$max=max(@score) || 'na';
		push @buffer, (($max ne 'na')?($score[($i=first_index{ $max } @score)], $tools[$i], ($max<.3)?'L':($max<.7)?'M':'H'):(('na')x3));
		print OUT join("\t", @buffer), "\n";
	}
close IN;
close OUT;
system "rm ${outdir}/${fname}*Temp $outfile; ln -s ${outdir}/${fname}.dbnsfp $outfile" and die $!;