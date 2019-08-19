package ngsutil;
use Exporter;
our @ISA = qw[ Exporter ];
our @EXPORT = qw[ &explode_varcall &varscan ];
use strict;
use warnings FATAL => qw[ numeric uninitialized ];
use List::Util qw[ sum ];

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
1;