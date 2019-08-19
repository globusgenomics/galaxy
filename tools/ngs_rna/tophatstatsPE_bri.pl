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

$samfile = $ARGV[0];
$fastqfile = $ARGV[1];
die "cannot open $fastqfile" unless (-e $fastqfile);
die "cannot open $samfile" unless (-e $samfile);

my $junk;
$divide = 4; # four lines per read in a fastq file
# count total number of raw reads
$results = `wc -l $fastqfile`;
($results, $junk) = split ' ', $results;
$total = $results / $divide;
print "$total\ttotal reads in fastq file\n";

# read through sam file, counting alignments
open IFILE, "$samfile" or die "cannot open $samfile\n";

%count = ();
while ($line = <IFILE>) 
{
    @line = split('\t', $line);
    $tempKey = $line[0];
    if(exists($count{$tempKey})){$count{$tempKey}++;}else{$count{$tempKey} = 1;}
}
@segs = keys %count;
$numAligned = @segs;
print "$numAligned\treads aligned in sam file\n";

$percAligned = $numAligned / $total *100;
printf("%.1f%%\taligned\n", $percAligned);

$multAligned = 0;
for (keys %count)
{
    if( $count{$_} > 1){$multAligned++;}
}

print "$multAligned\treads with multiple alignments\n";
$multPercAligned = $multAligned / $numAligned *100;
printf("%.1f%%\tof aligned segments had multiple alignments\n", $multPercAligned);


close IFILE;

exit;

############################################################


