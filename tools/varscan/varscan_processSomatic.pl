#!/usr/bin/perl


use strict;
use Cwd;

die qq(
Bad numbr of inputs

) if(!@ARGV);

my $options ="";
my $command="";
my $input="";
#my $output="";
my $working_dir = cwd();
my $log = '';
my %outfiles;
foreach my $arg (@ARGV) 
{
	my @tmp = split "::", $arg;
	if($tmp[0] eq "COMMAND") 
	{
		$command = $tmp[1];
	} 
	elsif($tmp[0] eq "INPUT")
	{
		$input = $tmp[1];
	}
	elsif($tmp[0] eq "OPTION") 
	{
		my @p = split(/\s+/,$tmp[1]);
		$options = "$options ${tmp[1]}";
	}
	elsif($tmp[0] eq "OUTPUT") 
	{
		my @p = split(/\s+/,$tmp[1]);
		$p[0] = substr($p[0],2);
		$outfiles{$p[0]} = $p[1];
	}
	
	elsif($tmp[0] eq "LOG")
        {
                $log = $tmp[1];
        }
	else 
	{
		die("Unknown Input: $input\n");
	}
}
system("ln -s $input $working_dir/in.dat");



## RUN
$command = "$command $working_dir/in.dat $options > $log 2>&1";
system($command);

## if tabular files are kept, write them to galaxy-datafile
if ($outfiles{'loh'} ne 'None') {
	## 
	system("cp '$working_dir/in.dat.LOH' '".$outfiles{'loh'}."'");
	system("cp '$working_dir/in.dat.LOH.hc' '".$outfiles{'loh'}."'");
	system("cp '$working_dir/in.dat.Germline' '".$outfiles{'germ'}."'");
	system("cp '$working_dir/in.dat.Germline.hc' '".$outfiles{'germ_hc'}."'");
	system("cp '$working_dir/in.dat.Somatic' '".$outfiles{'som'}."'");
	system("cp '$working_dir/in.dat.Somatic.hc' '".$outfiles{'som_hc'}."'");
}
## if vcf files are kept, generate them, and write to galaxy-datafile
if ($outfiles{'loh_hc_vcf'} ne 'None') {
	tab2vcf($working_dir,'in.dat.LOH.hc',$outfiles{'loh_hc_vcf'});
	tab2vcf($working_dir,'in.dat.Germline.hc',$outfiles{'germ_hc_vcf'});
	tab2vcf($working_dir,'in.dat.Somatic.hc',$outfiles{'som_hc_vcf'});
}
exit;

sub tab2vcf
{
	my $wd = shift;
	my $in = shift;
	my $out = shift;
	open OUT, ">$out";
	print OUT "##fileformat=VCFv4.1\n";
	print OUT "##source=processSomatic\n";
	print OUT '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth of quality bases">'."\n";
	print OUT '##INFO=<ID=SS,Number=1,Type=String,Description="Somatic status of variant (Germline,Somatic,LOH)">'."\n";
	print OUT '##INFO=<ID=VPV,Number=1,Type=Float,Description="Variant P-value">'."\n";
	print OUT '##INFO=<ID=SPV,Number=1,Type=Float,Description="Somatic Variant P-value">'."\n";
	print OUT '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'."\n";
	print OUT '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">'."\n";
	print OUT '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">'."\n";
	print OUT '##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">'."\n";
	print OUT '##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">'."\n";
	print OUT '##FORMAT=<ID=FREQ,Number=1,Type=String,Description="Variant allele frequency">'."\n";
	print OUT '##FORMAT=<ID=DP4,Number=1,Type=String,Description="Strand read counts: ref/fwd, ref/rev, var/fwd, var/rev">'."\n";
	print OUT '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NORMAL	TUMOR'."\n";
	open IN, "$wd/$in";
	my $head = <IN>;
	while (<IN>) {
		chomp;
		my @p = split(/\t/,$_);
		my $td = $p[4] + $p[5] + $p[8] + $p[9];
		my $vpv = sprintf('%.3E',$p[13]);
		my $spv = sprintf('%.3E',$p[14]);
		my $info = "DP=$td;SS=$p[12];VPV=$vpv;SPV=$spv";
		## gt normal string
		my $gtnormal = '';
		if ("$p[7]" eq "$p[2]") {
			$gtnormal = "0/0";
		}
		elsif ("$p[7]" eq "$p[3]") {
			$gtnormal = "1/1";
		}
		else {
			$gtnormal = "0/1";
		}
		my $nsum = $p[4] + $p[5];
		$gtnormal .= ":.:$nsum:$p[4]:$p[5]:$p[6]:$p[19],$p[20],$p[21],$p[22]";
		## gt tumor string
		my $gttumor = '';
		if ("$p[11]" eq "$p[2]") {
			$gttumor = "0/0";
		}
		elsif ("$p[11]" eq "$p[3]") {
			$gttumor = "1/1";
		}
		else {
			$gttumor = "0/1";
		}
		my $tsum = $p[8] + $p[9];
		$gttumor .= ":.:$tsum:$p[8]:$p[9]:$p[10]:$p[15],$p[16],$p[17],$p[18]";

		## outline 
		my $line = "$p[0]\t$p[1]\t.\t$p[2]\t$p[3]\t.\tPASS\t$info\tGT:GQ:DP:RD:AD:FREQ:DP4\t$gtnormal\t$gttumor\n";
		print OUT $line;
	}
	close IN;
	close OUT;
}
sub vs2vcf 
{

	#
	# G l o b a l     v a r i a b l e s 
	#
	my $version = '0.1';

	#
	# Read in file
	#
	my $input = shift;
	my $output = shift;
	my $chr_ord = shift;
	open(IN, $input) or die "Can't open $input': $!\n";
	open(OUT, ">$output") or die "Can't create $output': $!\n";
	my %output;

	while ( <IN> )
	{
		if ( /^#/ )
		{
			print OUT;
			next;
		}
		chomp;
		my $line = $_;

		my @flds = split ( "\t", $line );
		my $ref = $flds[3];
		my $alt = $flds[4];
		#
		# Deletion of bases
		#
		if ( $alt =~ /^\-/ )
		{
			($flds[3], $flds[4]) = ($ref.substr($alt,1), $ref);
		}

		#
		# Insertion of bases
		#
		if ( $alt =~ /^\+/ )
		{
			$flds[4] = $ref.substr($alt,1);
		}
		## insert dot for reference positions.
		if ($flds[4] eq '') {
			$flds[4] = '.';
		}
		print OUT join( "\t", @flds),"\n" unless defined $chr_ord;
		$output{$flds[0]}{$flds[1]} = join( "\t", @flds)."\n" if defined $chr_ord;
	}
	close(IN);
	# if chromosome order given return in sorted order
	if(defined $chr_ord) 
	{
		for my $chrom (@{ $chr_ord }) 
		{
			for my $pos (sort {$a<=>$b} keys %{ $output{$chrom} }) 
			{
				print OUT $output{$chrom}{$pos};
			}
		}
	}
	close(OUT);
}


sub chromosome_order 
{
	my $input = shift;
	# calculate flagstats
	my $COMM = "samtools view -H $input | grep '^\@SQ'";
	my @SQ = `$COMM`;
	chomp @SQ;
	for(my $i = 0; $i <= $#SQ; $i++) 
	{
		$SQ[$i] =~ s/^\@SQ\tSN:(.*?)\tLN:\d+$/$1/;
	} 
	return(@SQ);
}


