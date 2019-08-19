use strict;
use warnings;
use File::Basename; 
use Cwd;
use File::Path qw(make_path remove_tree);
die qq(
Bad numbr of inputs

) if(!@ARGV);

my $options ="";
my $normal="";
my $tumor="";
my $output="";


foreach my $input (@ARGV) 
{
	my @tmp = split "::", $input;
	if($tmp[0] eq "NORMAL") 
	{
		$normal = $tmp[1];
	} 
	elsif($tmp[0] eq "TUMOR") 
	{
		$tumor = $tmp[1];
	}
	elsif($tmp[0] eq "OUTPUT") 
	{
		$output = $tmp[1];
	}
	elsif($tmp[0] eq "OPTION") 
	{
		$options = "$options ${tmp[1]}";
	}
	
	else 
	{
		die("Unknown Input: $input\n");
	}
}


my $working_dir = cwd();

system ("ln -s $normal $working_dir/normal.bam");
system ("samtools index $working_dir/normal.bam");
 
system ("ln -s $tumor $working_dir/tumor.bam");
system ("samtools index $working_dir/tumor.bam");

system ("bam-somaticsniper $options $working_dir/tumor.bam $working_dir/normal.bam $output");


