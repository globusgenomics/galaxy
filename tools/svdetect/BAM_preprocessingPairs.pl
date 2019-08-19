#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;
my $version = '0.5_galaxy';

my $SAMTOOLS_BIN_DIR="/bioinfo/local/samtools";

my %opts = ( t=>1, p=>1, n=>1000000, f=>3, s=>0, S=>10000, o=>"." );

getopts('dt:p:n:f:s:S:o:b:l:x:N:', \%opts); #GALAXY 

my $working_dir=($opts{o} ne ".")? $opts{o}:"working directory";

my $pt_bad_mates_file=$opts{b};  #GALAXY 
my $pt_log_file=$opts{l}; #GALAXY 
my $pt_good_mates_file=$opts{x} if($opts{d}); #GALAXY 


die(qq/
    
Description:
    
    Preprocessing of mates to get anomalously mapped mate-pair\/paired-end reads as input
    for SVDetect.

    From all pairs mapped onto the reference genome, this script outputs abnormal pairs:
        - mapped on two different chromosomes
        - with an incorrect strand orientation and\/or pair order
        - with an insert size distance +- sigma threshold
    into a file <prefix.ab.bam\/sam>
    
    -BAM\/SAM File input format only, sorted by read names

    Version : $version
    SAMtools required for BAM files
    
    
Usage:   BAM_preprocessingPairs.pl [options] <all_mate_file.sorted.bam\/sam>

Options: -t BOOLEAN   read type: =1 (Illumina), =0 (SOLiD) [$opts{t}]
         -p BOOLEAN   pair type: =1 (paired-end), =0 (mate-pair)  [$opts{p}]
         -n INTEGER   number of pairs for calculating mu and sigma lengths [$opts{n}]
	 -s INTEGER   minimum value of ISIZE for calculating mu and sigma lengths [$opts{s}]
	 -S INTEGER   maximum value of ISIZE for calculating mu and sigma lengths [$opts{S}]
         -f REAL      minimal number of sigma fold for filtering pairs [$opts{f}]
         -d           dump normal pairs into a file [<prefix.norm.bam\/sam>] (optional)
	 -o STRING    output directory [$working_dir]

\n/) if (@ARGV == 0 && -t STDIN);

unless (-d $opts{o}){
	mkdir $opts{o} or die;
}
$opts{o}.="/" if($opts{o}!~/\/$/);

my $mates_file=shift(@ARGV);

$mates_file=readlink($mates_file);

my $bad_mates_file=(split(/\//,$mates_file))[$#_];

if($bad_mates_file=~/.(s|b)am$/){
    $bad_mates_file=~s/.(b|s)am$/.ab.sam/;
    $bad_mates_file=$opts{o}.$bad_mates_file;
}

else{
    die "Error: mate_file with the extension <.bam> or <.sam> needed !\n";
}

my $good_mates_file;
if($opts{d}){
    $good_mates_file=(split(/\//,$mates_file))[$#_];
    $good_mates_file=~s/.(b|s)am$/.norm.sam/;
    $good_mates_file=$opts{o}.$good_mates_file;
}

my $log_file=$opts{o}.$opts{N}.".svdetect_preprocessing.log"; #GALAXY 

#------------------------------------------------------------------------------#
#Calculate mu and sigma

open LOG,">$log_file" or die "$0: can't open ".$log_file.":$!\n";

print LOG "\# Calculating mu and sigma lengths...\n";
print LOG "-- file=$mates_file\n";
print LOG "-- n=$opts{n}\n";
print LOG "-- ISIZE min=$opts{s}, max=$opts{S}\n";

my ($record, $sumX,$sumX2) = (0,0,0);
my $warn=$opts{n}/10;
my $prev_pair="FIRST";

my $bam=($mates_file =~ /.bam$/)? 1:0;

if($bam){
    open(MATES, "${SAMTOOLS_BIN_DIR}/samtools view $mates_file |") or die "$0: can't open ".$mates_file.":$!\n";
}else{
    open MATES, "<".$mates_file or die "$0: can't open ".$mates_file.":$!\n";
}

while(<MATES>){
    
    my @t=split;
    
    next if ($t[0]=~/^@/);
    
    my $current_pair=$t[0];
    next if($current_pair eq $prev_pair);
    $prev_pair=$current_pair;                                                   
    
    my ($chr1,$chr2,$length)=($t[2],$t[6],abs($t[8]));
    
    next if (($t[1]&0x0004) || ($t[1]&0x0008));
    next if ($length<$opts{s} || $length>$opts{S}) ;
    
    if($chr2 eq "="){

        $sumX += $length;							#add to sum and sum^2 for mean and variance calculation
	$sumX2 += $length*$length;
        $record++;
    }

    if($record>$warn){
	print LOG "-- $warn pairs analysed\n";
        $warn+=$warn;
    }
    
    last if ($record>$opts{n});
    
}
close (MATES);

$record--;
my $mu = $sumX/$record;
my $sigma = sqrt($sumX2/$record - $mu*$mu); 

print LOG "-- Total : $record pairs analysed\n";
print LOG "-- mu length = ".decimal($mu,1).", sigma length = ".decimal($sigma,1)."\n";

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#Preprocessing pairs

$warn=100000;

$record=0;
my %count=( ab=>0, norm=>0, chr=>0, sense=>0, dist=>0, unmap=>0);

my $read_type=($opts{t})? "Illumina":"SOLiD";
my $pair_type=($opts{p})? "paired-end":"mate-paired";

print LOG "\# Preprocessing pairs...\n";
print LOG "-- file= $mates_file\n";
print LOG "-- type= $read_type $pair_type reads\n";
print LOG "-- sigma threshold= $opts{f}\n";
print LOG "-- using ".decimal($mu-$opts{f}*$sigma,4)."-".decimal($mu+$opts{f}*$sigma,4)." as normal range of insert size\n";

my @header;

if($bam){
    open(HEADER, "${SAMTOOLS_BIN_DIR}/samtools view -H $mates_file |") or die "$0: can't open ".$mates_file.":$!\n";
    @header=<HEADER>;
    close HEADER;
    open(MATES, "${SAMTOOLS_BIN_DIR}/samtools view $mates_file |") or die "$0: can't open ".$mates_file.":$!\n";
}else{
    open MATES, "<".$mates_file or die "$0: can't open ".$mates_file.":$!\n";
}

open AB, ">$bad_mates_file" or die "$0: can't write in the output: $bad_mates_file :$!\n";
print AB @header if($bam);

if($opts{d}){
    open NORM, ">$good_mates_file" or die "$0: can't write in the output: $good_mates_file :$!\n";
    print NORM @header if($bam);
}

$prev_pair="FIRST";
my $prev_bad;

while(<MATES>){
    
    my @t=split;
    my $bad=0;
    
    if ($t[0]=~/^@/){
        print AB;
        print NORM if ($opts{d});
        next;
    }
    
    my $current_pair=$t[0];
    if($current_pair eq $prev_pair){
        next if($prev_bad==-1);
        if($prev_bad){
            print AB;
        }elsif(!$prev_bad){
            print NORM if($opts{d});
        }
        next;
    }
    
    $prev_pair=$current_pair;
    
    my ($chr1,$chr2,$pos1,$pos2,$length)=($t[2],$t[6],$t[3],$t[7], abs($t[8]));
    
    if (($t[1]&0x0004) || ($t[1]&0x0008)){
        $prev_bad=-1;
        $count{unmap}++;
        $record++;
        next;
        
    }
    
    my $strand1 = (($t[1]&0x0010))? 'R':'F';
    my $strand2 = (($t[1]&0x0020))? 'R':'F';
    my $order1  = (($t[1]&0x0040))? '1':'2';
    my $order2  = (($t[1]&0x0080))? '1':'2';
    
    if($order1 == 2){
        ($strand1,$strand2)=($strand2,$strand1);
        ($chr1,$chr2)=($chr2,$chr1);
        ($pos1,$pos2)=($pos2,$pos1);
        ($order1,$order2)=($order2,$order1);
    }
    
    my $sense=$strand1.$strand2;
    
    if($chr1 ne "=" && $chr2 ne "="){
        $bad=1;
        $count{chr}++;
    }
    
    if($opts{p}){ #paired-end
        if(!(($sense eq "FR" && $pos1<$pos2) || ($sense eq "RF" && $pos2<$pos1))){
            $bad=1;
            $count{sense}++;
        }
    }else{ #mate-pair
        if($opts{t}){ #Illumina
            if(!(($sense eq "FR" && $pos2<$pos1) || ($sense eq "RF" && $pos1<$pos2))){
            $bad=1;
            $count{sense}++;
            }
        }else{ #SOLiD
            if(!(($sense eq "FF" && $pos2<$pos1) || ($sense eq "RR" && $pos1<$pos2))){
            $bad=1;
            $count{sense}++;
            }
        }
    }
    
    if(($chr1 eq "=" || $chr2 eq "=") && ($length <$mu - $opts{f}*$sigma || $length>$mu + $opts{f}*$sigma)){
        $bad=1;
        $count{dist}++;
    }
    
    if($bad){
        print AB;
        $count{ab}++;
        $prev_bad=$bad;
    }else{
        print NORM if ($opts{d});
        $count{norm}++;
        $prev_bad=$bad;
    }
    
    $record++;
    
    if($record>$warn){
        print LOG "-- $warn pairs analysed\n";
        $warn+=100000;
    }
}

close AB;
close NORM if($opts{d});

print LOG "-- Total : $record pairs analysed\n";
print LOG "-- $count{unmap} pairs whose one or both reads are unmapped\n";
print LOG "-- ".($count{ab}+$count{norm})." mapped pairs\n";
print LOG "---- $count{ab} abnormal mapped pairs\n";
print LOG "------ $count{chr} pairs mapped on two different chromosomes\n";
print LOG "------ $count{sense} pairs with incorrect strand orientation and\/or pair order\n";
print LOG "------ $count{dist} pairs with incorrect insert size distance\n";
print LOG "--- $count{norm} correct mapped pairs\n";

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#OUTPUT

if($bam){
    
    my $bam_file=$bad_mates_file;
    $bam_file=~s/.sam$/.bam/;
    print LOG "\# Converting sam to bam for abnormal mapped pairs\n";
    system("${SAMTOOLS_BIN_DIR}/samtools view -bS $bad_mates_file > $bam_file 2>".$opts{o}."samtools.log");
    unlink($bad_mates_file);
    print LOG "-- output created: $bam_file\n";

    system "rm $pt_bad_mates_file ; ln -s $bam_file $pt_bad_mates_file"; #GALAXY
    
    if($opts{d}){
        $bam_file=$good_mates_file;
        $bam_file=~s/.sam$/.bam/;
        print LOG "\# Converting sam to bam for correct mapped pairs\n";
        system("${SAMTOOLS_BIN_DIR}/samtools view -bS $good_mates_file > $bam_file 2>".$opts{o}."samtools.log");
        unlink($good_mates_file);
        print LOG "-- output created: $bam_file\n";

	system "rm $pt_good_mates_file ; ln -s $bam_file $pt_good_mates_file"; #GALAXY

    }

}

else{
    print LOG "-- output created: $bad_mates_file\n";
    print LOG "-- output created: $good_mates_file\n" if($opts{d});
}

close LOG;

system "rm $pt_log_file ; ln -s $log_file $pt_log_file"; #GALAXY


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub decimal{
    
  my $num=shift;
  my $digs_to_cut=shift;

  $num=sprintf("%.".($digs_to_cut-1)."f", $num) if ($num=~/\d+\.(\d){$digs_to_cut,}/);

  return $num;
}
#------------------------------------------------------------------------------#
