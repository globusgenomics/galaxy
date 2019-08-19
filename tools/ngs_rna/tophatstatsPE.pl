#!/usr/bin/perl -w

############################################
# tophatstatsPE
# John Garbe
# May 29, 2012
#
# Calculate statistics about the data going into and out of tophat:
# -Count the number of raw reads, the number of aligned reads, the number 
# of properly paired reads, etc
#
###########################################

die "USAGE: tophatstats.pl accepted_hits.bam fastqfile.fastq\n" unless ($#ARGV == 1);

$bamfile = $ARGV[0];
$fastqfile = $ARGV[1];
die "cannot open $fastqfile" unless (-e $fastqfile);
die "cannot open $bamfile" unless (-e $bamfile);

my $junk;
$divide = 4; # four lines per read in a fastq file
$insertsum = 0;
$insertcount = 0;
# count total number of raw reads
if (1) {
    print "Input files: $bamfile\t $fastqfile\n";
    $results = `wc -l $fastqfile`;
    ($results, $junk) = split ' ', $results;
    $total = $results / $divide;
    print "$total total read pairs in fastq file\n";
}

# read in bam file as a pipe from samtools
open IFILE, "samtools view $bamfile |" or die "cannot open $bamfile\n";

# initialize array counting multiple alignments
for $i (0..200) {
    $count[$i] = 0;
    for $j (0..40) {
	$unique[$i][$j] = 0;
    }
}

$insertsum = 0;
$insertcount = 0;
# read through bam file, counting alignments
while ($line = <IFILE>) {
    @line = split /\t/, $line, 3;
    $count[$line[1]]++;

    # pull out number of alignments for this read
    @line = split /\t/, $line, 12;
    ($junk, $int) = split /NH:i:/, $line[11];
    ($int, $junk) = split /\W/, $int;
    if ($int == 1) { # count unique alignments
	$unique[$line[1]][1]++;
    } elsif ($int > 1) { # count non-unique alignments
	$unique[$line[1]][2]+= 1 / $int;
    }
    if ($line[1] == 99 || $line[1] == 83) { # if proper insert size use for insert size average
	@line = split /\t/, $line, 10;
	$insertsum += $line[8];
	$insertcount++;
    }
}

# singleton
$singleton_unique = 0;
$singleton_multiple = 0;
for $i (73, 133, 89, 121, 165, 181, 101, 117, 153, 185, 69, 137) {
    $singleton_unique += $unique[$i][1];
    $singleton_multiple += int($unique[$i][2]);
}
$singleton = $singleton_unique + $singleton_multiple;
$singleton_percent = $singleton / $total * 100;
printf "$singleton (%.2f%%) read pairs with only one read in the pair mapped ($singleton_unique with unique alignments)\n", $singleton_percent;

# correct correct
$correctcorrect_unique = 0;
$correctcorrect_multiple = 0;
for $i (99, 83) {
    $correctcorrect_unique += $unique[$i][1];
    $correctcorrect_multiple += int($unique[$i][2]);
}
$correctcorrect = $correctcorrect_unique + $correctcorrect_multiple;
$correctcorrect_percent = $correctcorrect / $total * 100;
printf "$correctcorrect (%.2f%%) read pairs mapped with correct orientation and insert size ($correctcorrect_unique with unique alignments)\n", $correctcorrect_percent;

# correct wrong insert size
$correctwrong_unique = 0;
$correctwrong_multiple = 0;
for $i (81, 97, 65, 113) {
    $correctwrong_unique += $unique[$i][1];
    $correctwrong_multiple += int($unique[$i][2]);
}
$correctwrong = $correctwrong_unique + $correctwrong_multiple;
$correctwrong_percent = $correctwrong / $total * 100;
printf "$correctwrong (%.2f%%) read pairs mapped with correct orientation but wrong insert size ($correctwrong_unique with unique alignments)\n", $correctwrong_percent;

# correct wrong orientation
$correctswitched_unique = 0;
$correctswitched_multiple = 0;
for $i (67, 115) {
    $correctswitched_unique += $unique[$i][1];
    $correctswitched_multiple += $unique[$i][2];
}
$correctswitched = $correctswitched_unique + $correctswitched_multiple;
$correctswitched_percent = $correctswitched / $total * 100;
printf "$correctswitched (%.2f%%) read pairs mapped with wrong orientation but correct insert size ($correctswitched_unique with unique alignments)\n", $correctswitched_percent;

# no mapping
$nomapping = $total - ($correctcorrect + $correctswitched + $correctwrong + $singleton);
$nomapping_percent = $nomapping / $total * 100;
printf "$nomapping (%.2f%%) read pairs with no mapping\n", $nomapping_percent;

# insert size
$insertavg = 0;
$insertavg = $insertsum / $insertcount unless ($insertcount == 0);
#printf "%.2fbp average inner distance between read pairs (of read pairs with correct insert size)\n", $insertavg; # Not ready for primetime yet - needs verification

exit;

############################################################

