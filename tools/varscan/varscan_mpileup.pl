#!/usr/bin/perl

use strict;
use Cwd;

die qq(
Bad numbr of inputs

) if(!@ARGV);

my $options ="";
my $file="";
my $command="";
my $output="";
my $working_dir = cwd();
my $temp_vcf = "$working_dir/temp";
my $log="";

foreach my $input (@ARGV) 
{
	my @tmp = split "::", $input;
	if($tmp[0] eq "COMMAND") 
	{
		$command = $tmp[1];
	} 
	elsif($tmp[0] eq "INPUT") 
	{
		$file = $tmp[1];
	}
	elsif($tmp[0] eq "OPTION") 
	{
		my @p = split(/\s+/,$tmp[1]);
		if ($p[0] =~ m/variants|strand-filter|output-vcf/ && $p[1] == 0) {
			next;
		}

		$options = "$options ${tmp[1]}";
	}
	elsif($tmp[0] eq "OUTPUT") 
	{
		$output = $tmp[1];
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

system ("$command $file $options 1>$temp_vcf 2>$log");

vs2vcf($temp_vcf, $output);


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

