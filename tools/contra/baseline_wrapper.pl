#sshpass -p pmac1512 ssh -o StrictHostkeyChecking=no galaxy@pmc-bioinf03 "gmt music play $*"
#echo "gmt music play $*"


use strict;
use warnings;
use File::Basename; 
use Cwd;
die qq(
Bad numbr of inputs

) if(!@ARGV);


my $player_options = "";
my $baseline_output;

my $dir = getcwd;
my $variable = "";
my $files = "";
foreach my $input (@ARGV) 
{
	my @tmp = split "::", $input;
	
	if($tmp[0] eq "PLAYEROPTION") 
	{
		$variable = $tmp[1];
		$variable =~ s/=/ /g;
		$player_options = "$player_options $variable";
	}
	elsif($tmp[0] eq "BASELINEOUTPUT") 
	{
		$baseline_output = $tmp[1];
	}  
	elsif($tmp[0] eq "BAMLISTENTRY") 
	{
		$files = "$files ${tmp[1]}";
	}
        elsif($tmp[0] eq "BAMLISTFILE")
        {
                open(FH, '<:encoding(UTF-8)', $tmp[1]);
                while (my $line = <FH>)
                {
                    chomp($line);
                    $files .= $line . " ";
                }
                close(FH)
        }
	else 
	{
		die("Unknown Input: $input\n");
	}
}


my $working_dir = "BASELINE_OUTPUT";
#remove extension

#Create Contra Output dir 
if (! -d $working_dir)
{
    system ("mkdir $working_dir");
}

#run baseline

system ("baseline.py --files $files --output $working_dir $player_options > /dev/null");

#Search control file in output dir
opendir(DIR, $working_dir);
my @FILES= readdir(DIR); 
foreach my $file (@FILES) 
{
	my ($filename,$directory,$extension) = fileparse($file, qr/\.[^.]*/);
	if ($extension eq ".txt")
	{
		system ("mv $working_dir/$file $baseline_output");
	}
}
closedir(DIR);

