use strict;
use warnings;

use File::Find;
use File::Path qw(make_path);
use File::Spec;


die "Bad number of inputs" if(!@ARGV);

my $player_options = "";
my ($contra_output,$contra_dir,$contra_vcf,$contra_txt);


foreach my $input (@ARGV) 
{
	my @tmp = split "::", $input;
	
	if($tmp[0] eq "PLAYEROPTION") 
	{
		my $variable = $tmp[1];
		$variable =~ s/=/ /g;
		
		$player_options = "$player_options $variable";
	}
	elsif($tmp[0] eq "CONTRAOUTPUT") 
	{
		$contra_output = $tmp[1];
	}  
	elsif($tmp[0] eq "CONTRADIR") 
	{
		$contra_dir = $tmp[1];
	}  
        elsif($tmp[0] eq "CONTRAVCF")
        {
                $contra_vcf = $tmp[1];
        }
        elsif($tmp[0] eq "CONTRATXT")
        {
                $contra_txt = $tmp[1];
        }
	else 
	{
		die("Unknown input: $input\n");
	}
}



my $working_dir = File::Spec->catfile($contra_dir, 'CONTRA_OUTPUT');
make_path($contra_dir);
#remove extension

my $cmd = "contra.py -o $working_dir $player_options";
print "Command to be executed: $cmd\n";
system($cmd);


# set the output to the correct locations
if (-f $File::Find::name)

#set html
#print "$contra_output - $working_dir\n";
open(HTML, ">$contra_output");
print HTML "<html><head><title>Contra: Copy Number Analysis for Targeted Resequencing</title></head><body><h3>Contra Output Files:</h3><p><ul>\n";
find({wanted => \&add_file, preprocess => sub {sort @_}}, $working_dir);
print HTML "</ul></p>\n";
close(HTML);

sub add_file {
	if (-f $File::Find::name) {
		my $rel_path = File::Spec->abs2rel($File::Find::name, $working_dir);
		print ("adding $rel_path\n");
		print HTML "<li><a href=\"CONTRA_OUTPUT/$rel_path\">$rel_path</a></li>\n";
 	}
}

