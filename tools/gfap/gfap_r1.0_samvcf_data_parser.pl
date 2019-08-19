#!/usr/bin/perl

use strict;
# use lib 'inc/perlmod';
# use ngsutil qw[ :DEFAULT &explode_varcall ];
use warnings FATAL => qw[ numeric uninitialized ];
use List::Util qw[ sum min max ];
use File::Basename;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#	PATH TO YOUR R-bin DIRECTORY
my $rbin = 'R';
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

my $annovar_dir = '/nfs/software/galaxy/tools/gfap/inc/annovar';
my $rdep = '/nfs/software/galaxy/tools/gfap/inc/R';

my($varfile, $outdir, $outfile, $i, @DP4, @buffer, @Temp, @previous, @fnames, %opts, %chr);

GetOptions(\%opts, "varfile=s", "outdir=s", "outfile=s");
$varfile = $opts{varfile};
$outdir  = $opts{outdir};
$outfile = $opts{outfile};

my $fname = readlink($varfile) || $varfile;
$fname = basename($fname);

my %fh=(
	'chr1' => *chr1,	'chr2' => *chr2,	'chr3' => *chr3,	'chr4' => *chr4,	'chr5' => *chr5,
	'chr6' => *chr6,	'chr7' => *chr7,	'chr8' => *chr8,	'chr9' => *chr9,	'chr10' => *chr10,
	'chr11' => *chr11,	'chr12' => *chr12,	'chr13' => *chr13,	'chr14' => *chr14,	'chr15' => *chr15,
	'chr16' => *chr16,	'chr17' => *chr17,	'chr18' => *chr18,	'chr19' => *chr19,	'chr20' => *chr20,
	'chr21' => *chr21,	'chr22' => *chr22,	'chrX' => *chrX,	'chrY' => *chrY,	'chrM' => *chrM
);

`${annovar_dir}/convert2annovar.pl -format vcf4 $varfile -includeinfo > ${outdir}/${fname}_Temp-00 2> /dev/null` and die $!;

open($fh{$_}, ">${outdir}/${fname}_${_}.Temp-00") or die $! foreach keys %fh;
open IN, "<${outdir}/${fname}_Temp-00" or die $!;
while(<IN>){
		/^(\S+)\s+(?:\S+\s+){2}(\S+)\s+(\S+)/;
		next if !exists $fh{$1};
		if(min(length($2), length($3))!=1){
				chomp;
				@buffer=split /\s+/, $_;
				@Temp=explode_varcall(@buffer[1,3..4]);
				for($i=0; $i<$#{$Temp[0]}; $i+=2){
						print{ $fh{$buffer[0]} } join("\t", $buffer[0], @{$Temp[0]}[$i..$i+1], @{$Temp[1]}[$i..$i+1], @buffer[6..$#buffer]), "\n";
					}
				next;
			}
		print{ $fh{$1} } $_;
		$chr{$1}++;
	}
close IN;
foreach (keys %fh){
		close($fh{$_});
		next if !exists $chr{$_};
		`sort -k2,2n -k3,3n ${outdir}/${fname}_${_}.Temp-00 > ${outdir}/${fname}_${_}.Temp-01` and die $!;
		open IN, "<${outdir}/${fname}_${_}.Temp-01" or die $!;
		open OUT, ">${outdir}/${fname}_${_}.Temp-02" or die $!;
		$_=readline(IN);
		/^((?:\S+\s+){7})(?:\S+\s+){8}(\S+\s+\S+)/;
		@buffer=split /\s+/, $1.$2;
		($_=pop(@buffer))=~s/.+DP4=([^;]+).+/$1/;
		@DP4=split /,/, $_;
		push @buffer, @DP4;
		@previous=@buffer;
		MAINLOOP: while(<IN>){
				/^((?:\S+\s+){7})(?:\S+\s+){8}(\S+\s+\S+)/;
				@buffer=split /\s+/, $1.$2;
				($_=pop(@buffer))=~s/.+DP4=([^;]+).+/$1/;
				@DP4=split /,/, $_;
				push @buffer, @DP4;
				while(($previous[0] eq $buffer[0]) && ($buffer[2]==$previous[2]+1) && (join('', @previous[3..4]) !~ /-/) && (join('', @buffer[3..4]) !~ /-/)){
						$previous[2]=$buffer[2];
						$previous[$_].=$buffer[$_] for 3..4;
						$previous[5]='unk' if $previous[5] ne $buffer[5];
						$previous[7]='SKIP' if $previous[7] ne $buffer[7];
						for (6){
								$previous[$_]+=$buffer[$_];
								$previous[$_]/=2;
							}
						next MAINLOOP;
					}
				$previous[7]='NONE' if $previous[7] eq '.';
                                #print join ("\t", @buffer) . "\n"; 
                                for (my $i=6;$i<=11; $i++)
                                {
                                    next if ($i eq '7');
                                #    print $i . "\t" . $previous[$i] . "\n";
                                #    if (looks_like_number($previous[$i]))
                                #    {
                                     $previous[$i]=$previous[$i];
				#        $previous[$i]=sprintf("%.0f", $previous[$i]);
                                #    }
                                }
				print OUT join("\t", @previous[0..6,8,7]), "\n";
				@Temp=@previous if eof;
				@previous=@buffer;
			}
		$previous[7]='NONE' if $previous[7] eq '.';
		#$previous[$_]=sprintf("%.0f", $previous[$_]) for (6,8);
		print OUT join("\t", @previous[0..6,8,7]), "\n" if(join('_', @Temp[1..2]) ne join('_', @previous[1..2]));
		close IN;
		close OUT;
	}
foreach (1..22, 'X', 'Y', 'M'){
		push @fnames, "${outdir}/${fname}_chr${_}.Temp-02" if exists $chr{"chr$_"};
	}
system join(' ', 'cat', @fnames, '>', "${outdir}/${fname}.Temp.2R") and die $!;
`${rbin} --vanilla --slave --args ${outdir}/${fname}.Temp.2R < ${rdep}/samvcf_data_parser.R` and die $!;
system "rm ${outdir}/${fname}*Temp* $outfile; ln -s ${outdir}/${fname}.var $outfile" and die $!;
