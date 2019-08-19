#!/usr/bin/perl 


use Data::Dumper;
use Getopt::Long;
use Pod::Usage;



#pod2usage(-verbose => 1) if ($help == 1);
#if (@ARGV == 0) {
#	pod2usage(-msg => "Invalid number of arguments!", -exitval => 2, -verbose => 2);
#}

#my $rsem_version = "/opt/rsem-1.1.17";
my $minL = 1;
my $maxL = 1000;
my $NMB = 1024;

# Extra file output #beta
# --isoformfile $isoforms 
# --thetafile $theta 
# --cntfile $cnt 
# --modelfile $model 
# --bamfile $bam_res

GetOptions(
    "log=s"                  => \$log,
    "bam_genome=s"           => \$bam_genome,
    "bamtype=s"              => \$bamtype,
    "isoformfile=s"          => \$isoforms,
    "reference=s"            => \$dbref,
    "sampling-for-bam=s"     => \$samplingbam,
    "thetafile=s"            => \$theta,
    "cntfile=s"              => \$cnt,
    "modelfile=s"            => \$model,
    "bamfile=s"              => \$bamfile,
    "output=s"               => \$output,
    "single_fasta=s"         => \$single_fasta,
    "fasta1=s"               => \$fasta1,
    "fasta2=s"               => \$fasta2,
    "single_fastq=s"         => \$single_fastq,
    "fastq1=s"               => \$fastq1,
    "fastq2=s"               => \$fastq2,
    "no-qualities"           => \$no_qual,
    "paired-end"             => \$paired_end,
    "sam"                    => \$is_sam,
    "bam"                    => \$is_bam,
    "sam-header-info=s"      => \$fn_list,
    "tag=s"                  => \$tagName,
    "seed-length=i"          => \$L,
    "bowtie-path=s"          => \$bowtie_path,
    "bowtie-n=i"             => \$C,
    "bowtie-e=i"             => \$E,
    "bowtie-m=i"             => \$maxHits,
    "phred33-quals"          => \$phred33,
    "phred64-quals"          => \$phred64,
    "solexa-quals"           => \$solexa,
    "forward-prob=f"         => \$probF,
    "fragment-length-min=i"  => \$minL,
    "fragment-length-max=i"  => \$maxL,
    "fragment-length-mean=f" => \$mean,
    "fragment-length-sd=f"   => \$sd,
    "estimate-rspd=s"        => \$estRSPD,
    "num-rspd-bins=i"        => \$B,
    "p|num-threads=i"        => \$nThreads,
    "output-genome-bam"      => \$genBamF,
    "calc-ci=s"              => \$calcCI,
    "ci-memory=i"            => \$NMB,
    "time"                   => \$mTime,
    "q|quiet"                => \$quiet,
    "keep_unmapped1=s"       => \$keep_unmapped_1,
    "keep_unmapped2=s"       => \$keep_unmapped_2,
) or pod2usage( -exitval => 2, -verbose => 2 );

#check parameters and options

if ($is_sam || $is_bam) {
    pod2usage(-msg => "from rsem-wrapper->Invalid number of arguments!", -exitval => 2, -verbose => 2) if (scalar(@ARGV) != 4);
    pod2usage(-msg => "--sam and --bam cannot be active at the same time!", -exitval => 2, -verbose => 2) if ($is_sam == 1&& $is_bam == 1);
    pod2usage(-msg => "--bowtie-path, --bowtie-n, --bowtie-e, --bowtie-m, --phred33-quals, --phred64-quals or --solexa-quals cannot be set if input is SAM/BAM format!", -exitval => 2, -verbose => 2) if ($bowtie_path ne "" || $C != 2 || $E != 99999999 || $maxHits != 200 || $phred33 || $phred64 || $solexa);
}
#else {
#    pod2usage(-msg => "from rsem-wraper->Invalid number of arguments!", -exitval => 2, -verbose => 2) 
#	if (!$paired_end && scalar(@ARGV) != 1 || $paired_end && scalar(@ARGV) != 1);    
#    pod2usage(-msg => "Only one of --phred33-quals --phred64-quals/--solexa1.3-quals --solexa-suqls can be active!", -exitval => 2, -verbose => 2) if ($phred33 + $phred64 + $solexa > 1);    
#    podwusage(-msg => "--sam , --bam or --sam-header-info cannot be set if use bowtie aligner to produce alignments!", -exitval => 2, -verbose => 2) if ($is_sam | $is_bam || $fn_list ne "");
#}

pod2usage(-msg => "Forward probability should be in [0, 1]!", -exitval => 2, -verbose => 2) if ($probF < 0 || $probF > 1);
pod2usage(-msg => "Min fragment length should be at least 1!", -exitval => 2, -verbose => 2) if ($minL < 1);
pod2usage(-msg => "Min fragment length should be smaller or equal to max fragment length!", -exitval => 2, -verbose => 2) if ($minL > $maxL);
pod2usage(-msg => "The memory allocated for calculating credibility intervals should be at least 1 MB!\n", -exitval => 2, -verbose => 2) if ($NMB < 1);
pod2usage(-msg => "Number of threads should be at least 1!\n", -exitval => 2, -verbose => 2) if ($nThreads < 1);

# IO Redirection to log file
use IO::Handle;
open OUTPUT, '>', $log or die "cant open file $log $!\n";;
open ERROR,  '>>', $log  or die "cant open file $log $!\n";
STDOUT->fdopen( \*OUTPUT, 'w' ) or die "cant open file $!\n";
STDERR->fdopen( \*ERROR,  'w' ) or die "cant open file $!\n";
# 

# print timestamp
my $start = `date`;
print "Start time: $start\n";

my @options;

# generates new output called sample_name.genome.bam 
# with alignments
# mapped to genomic coordinates and annotated with their posterior
# probabilities. In addition, RSEM will call samtools (included in
# RSEM package) to sort and index the bam file.
# 'sample_name.genome.sorted.bam' and
# 'sample_name.genome.sorted.bam.bai' will be generated. (Default: off)

if ($bamtype eq "yes") {
    my $bam_genome_par = "--output-genome-bam";
    push @options, $bam_genome_par;
}
if ($samplingbam eq "yes") {
    my $samplingbam = "--sampling-for-bam";
    push @options, $samplingbam;
}
if ($estRSPD eq "yes") {
    my $rspd = "--estimate-rspd";
    push @options, $rspd;
}
$probF = "--forward-prob $probF";
push @options, $probF;

if ($calcCI eq "yes") {
    my $calcCI = "--calc-ci";
    push @options, $calcCI;
    my $cimem = "--ci-memory $NMB";
    push @options, $cimem;
}
if ($tagName) {
    my $tagName = "--tag $tagName";
    push @options, $tagName;
}
if ($L) {
	my $L = "--seed-length $L";
	push @options, $L;
}
if ($C) {
    my $C = "--bowtie-n $C";
    push @options, $C;
}
if ($E) {
    my $E = "--bowtie-e $E";
    push @options, $E;
}
if ($maxHits) {
    my $maxHits = "--bowtie-m $maxHits";
    push @options, $maxHits;
}
if ($minL != 1) {
    my $minL = "--fragment-length-min $minL";
    push @options, $minL;
}
if ($maxL != 1000) {
    my $maxL = "--fragment-length-max $maxL";
    push @options, $maxL;
}
if ($mean) {
    my $mean = "--fragment-length-mean $mean";
    push @options, $mean;
}
if ($sd) {
    my $sd = "--fragment-length-sd $sd";
    push @options, $sd;
}
if ($keep_unmapped_1) {
    my $keep_unmap = "--keep-intermediate-files";
    push @options, $keep_unmap;
}
my $options=  join(" ", @options);

# make temporary directory where output will be written to. Must be deleted at end of run
my $wd_tmpdir = $output."_rsemtmpdir"; ## make tmpdir
unless(-e $wd_tmpdir or mkdir $wd_tmpdir) {
    die "Unable to create $wd_tmpdir";
}
my $output_tmp = "$wd_tmpdir/rsemout.dat";

#BUILD COMMAND BASED ON PARSED OPTIONS
if ($no_qual) { 
        #reads are in fasta file format
	if ($paired_end) { # reads are in paired end
	    my $cmd = "rsem-calculate-expression --quiet --no-qualities --paired-end -p $nThreads $options $fasta1 $fasta2 $dbref $output_tmp";
 	    print "RSEM Parameters used by Galaxy:\n$cmd\n";
	    system($cmd);
	}
	#run single end with one fasta file
	else {
	my $cmd = "rsem-calculate-expression --quiet --no-qualities -p $nThreads $options $single_fasta $dbref $output_tmp";
 	print "RSEM Parameters used by Galaxy:\n$cmd\n";
	system($cmd);
	}
}
else {
    # reads are in fastq file format
    # type of fastq file?
	my $fastqtype;
	if ($phred33) {
		$fastqtype = "--phred33-quals";
	}
	elsif ($phred64) {
		$fastqtype = "--phred64-quals";
	}
	elsif ($solexa) {
		$fastqtype = "--solexa-quals";
	}
	if ($paired_end) { 
	#reads in paired end 
        #run paired end with two fasq files
	    my $cmd = "rsem-calculate-expression --quiet --paired-end -p $nThreads $options $fastqtype $fastq1 $fastq2 $dbref $output_tmp";
 	    print "RSEM Parameters used by Galaxy:\n$cmd\n";
	    system($cmd);
	}
	else { 
	my $cmd = "rsem-calculate-expression --quiet -p $nThreads $options $fastqtype $single_fastq $dbref $output_tmp";
 	print "RSEM Parameters used by Galaxy:\n$cmd\n";
	system($cmd);
	}
}

 #Rename files for galaxy
my $mv_genes = "mv $output_tmp.genes.results $output";
my $mv_isoforms = "mv $output_tmp.isoforms.results $isoforms";

#print "bamtype-parameter=$bamtype\n";
my $mv_bam_transcript;
my $mv_bam_genome;
if ($bamtype eq "yes") {
#  $mv_bam_genome = "mv $output.genome.sorted.bam $bam_genome";
  $mv_bam_genome = "mv $output_tmp.genome.bam $bam_genome";
  system($mv_bam_genome);
}

#$mv_bam_transcript = "mv $output.transcript.sorted.bam $bamfile";
$mv_bam_transcript = "mv $output_tmp.transcript.bam $bamfile";

my @rsem_dir = split(/\//, $output_tmp);
my $short_output = $rsem_dir[-1];
my $mv_theta = "mv $output_tmp.stat/$short_output.theta $theta";
my $mv_cnt = "mv $output_tmp.stat/$short_output.cnt $cnt";
my $mv_model = "mv $output_tmp.stat/$short_output.model $model";

my ($mv_unmapped_1, $mv_unmapped_2);
if ($keep_unmapped_1)
{
   $mv_unmapped1 = "mv $output_tmp.temp/$short_output" . "_un_1.fq $keep_unmapped_1";
   system($mv_unmapped1);
}
if ($keep_unmapped_2)
{
   $mv_unmapped2 = "mv $output_tmp.temp/$short_output" . "_un_2.fq $keep_unmapped_2";
   system($mv_unmapped2);
}


system($mv_genes);
system($mv_isoforms);
system($mv_bam_transcript);

# print timestamp
my $end = `date`;
print "End time: $end\n";

system($mv_theta);
system($mv_cnt);
system($mv_model);
print "LOG $mv\n";

## delete the temporary directory
my $rm_tmp = "rm -rf $wd_tmpdir";
system($rm_tmp)

