#!/usr/bin/perl -w


#####################################################
# On January 8th,2018, fixed an issue with the agrep counting portion.
# Originally, the pipeline was 'catting' the kmer file with the option '-n'
# to report the line on which spacers sequences were found.  In one instances,
# agrep somehow ignore the line # in the output file, causing the split function to
# be unphased when getting the count number.  I correct the issue by removing the '-n'
# option in the catting portion.
#####################################################

#### Changelog #########

### spoligo_kmers_parralle_v1.0.pl ##
## Agrep has an bug that sometime falsely giving 0 counts, resulting in an error in
## the spoligotype prediction.  In order to make the prediction a bit more robust,
## the script was changed to 1st perform a strait grep search.  If that search comes back
## empty, agrep is then performed

use LWP::Simple;
use Parallel::ForkManager;
use strict;

my $Read_R1 = $ARGV[0];
my $Read_R2 = $ARGV[1];
my $run_path = $ARGV[2];

#my $jellyfish = "/mnt/galaxyTools/tools/jellyfish/2.2.10/bin/jellyfish";
my $spoligo_data = "/opt/galaxy/tools/nysdh_tools/spoligotyper/spoligo_data.txt";
my $NY_spoligotypes = "/opt/galaxy/tools/nysdh_tools/spoligotyper/NY_spolygo_types.txt";
### run jellyfish version 2.2.4
print "\nzcat $Read_R1 $Read_R2 > $run_path/tmp_kmer.fastq\n\n";
print "\njellyfish count -C -m 25 -s 10M -t 12 $run_path/tmp_kmer.fastq -o $run_path/mer_counts.jf \n\n";
print "\n jellyfish dump $run_path/mer_counts.jf -c -t -L 4 > $run_path/mers_count_25.out\n\n";
system ("zcat $Read_R1 $Read_R2 > $run_path/tmp_kmer.fastq");
system ("jellyfish count -C -m 25 -s 10M -t 12 $run_path/tmp_kmer.fastq -o $run_path/mer_counts.jf");
system ("jellyfish dump $run_path/mer_counts.jf -c -t -L 4 > $run_path/mers_count_25.out");


# hash for bin to octa
my %octa;
$octa{"000"} = 0;
$octa{"001"} = 1;
$octa{"010"} = 2;
$octa{"011"} = 3;
$octa{"100"} = 4;
$octa{"101"} = 5;
$octa{"110"} = 6;
$octa{"111"} = 7;



# Spacer sequences
# full spacer 3 = CCGTGCTTCCAGTGATCGCCTTCTA
# example agrep with 2 missmatches -> system ("zcat $Read_R1 $Read_R2 | agrep -2 -c  $query > count.out"); 
my @seq_spacers=('ATAGAGGGTCGCCGGCTCTGGATCA','CCTCATGCTTGGGCGACAGCTTTTG','CCGTGCTTCCAGTGATC','ACGTCATACGCCGACCAATCATCAG','TTTTCTGACCACTTGTGCGGGATTA','CGTCGTCATTTCCGGCTTCAATTTC','GAGGAGAGCGAGTACTCGGGGCTGC','CGTGAAACCGCCCCCAGCCTCGCCG','ACTCGGAATCCCATGTGCTGACAGC','TCGACACCCGCTCTAGTTGACTTCC','GTGAGCAACGGCGGCGGCAACCTGG','ATATCTGCTGCCCGCCCGGGGAGAT','GACCATCATTGCCATTCCCTCTCCC','GGTGTGATGCGGATGGTCGGCTCGG','CTTGAATAACGCGCAGTGAATTTCG','CGAGTTCCCGTCAGCGTCGTAAATC','GCGCCGGCCCGCGCGGATGACTCCG','CATGGACCCGGGCGAGCTGCAGATG','TAACTGGCTTGGCGCTGATCCTGGT','TTGACCTCGCCAGGAGAGAAGATCA','TCGATGTCGATGTCCCAATCGTCGA','ACCGCAGACGGCACGATTGAGACAA','AGCATCGCTGATGCGGTCCAGCTCG','CCGCCTGCTGGGTGAGACGTGCTCG','GATCAGCGACCACCGCACCCTGTCA','CTTCAGCACCACCATCATCCGGCGC','GGATTCGTGATCTCTTCCCGCGGAT','TGCCCCGGCGTTTAGCGATCACAAC','AAATACAGGCTCCACGACACGACCA','GGTTGCCCCGCGCCCTTTTCCAGCC','TCAGACAGGTTCGCGTCGATCAAGT','GACCAAATAGGTATCGGCGTGTTCA','GACATGACGGCGGTGCCGCACTTGA','AAGTCACCTCGCCCACACCGTCGAA','TCCGTACGCTCGAAACGCTTCCAAC','CGAAATCCAGCACCACATCCGCAGC','CGCGAACTCGTCCACAGTCCCCCTT','CGTGGATGGCGGATGCGTTGTGCGC','GACGATGGCCAGTAAATCGGCGTGG','CGCCATCTGTGCCTCATACAGGTCC','GGAGCTTTCCGGCTTCTATCAGGTA','ATGGTGGGACATGGACGAGCGCGAC','CGCAGAATCGCACCGGGTGCGGGAG');

## import Spoligo SIT # and Clades
my %spoligo; ## hash to store CDC spolygotype
my %ny_spoligo; ## hash to store CDC spolygotype
open (FILE,$spoligo_data);  ## CDC spoligotype file
while (my $in_spoligo_data=<FILE>) {
	chomp $in_spoligo_data;
	my @split_spoligo_data = split (' ',$in_spoligo_data);
	my @split_line = split ('\) ',$in_spoligo_data);
	if ($in_spoligo_data =~ '\!|\"') {
		$spoligo{$split_spoligo_data[2]} = "$split_spoligo_data[0]\t$split_line[-1]";
	}
}
close FILE;


open (FILE,$NY_spoligotypes);  ## NYS spoligotype file
my $in_ny_file=<FILE>;
while ($in_ny_file=<FILE>) {
	chomp $in_ny_file;
	my @split_line = split ('\t',$in_ny_file);
	$ny_spoligo{$split_line[3]} = "$in_ny_file";
}
close FILE;





### Read mapping, SNP calling and Consensus sequence ###
open (OUT,">$run_path/spoligo.out");
&spoligo;
#system ("rm $run_path/count.out");

my $counter=0;
do {
	#system ("rm $run_path/$counter\_count.out");
	#system ("rm $run_path/$counter\_rcount.out");
	++$counter;
}
until $counter==43;
#system ("rm $run_path/tmp_kmer.fastq");
#system ("rm $run_path/mers_count_25.out");
#system ("rm $run_path/mer_counts.jf");

exit;




sub spoligo {
	my $triplet="";
	my $octa_num="";
	my $stringx="";

	my $control="TCTCGGGGTTTTGGGTCTGACGACC"; 

	my $reverse_control = reverse $control;
	$reverse_control =~ tr/ATCG/TAGC/;

        print "cat $run_path/mers_count_25.out | grep -E '$control|$reverse_control' > $run_path/count.out\n"; 
	system ("cat $run_path/mers_count_25.out | grep -E '$control|$reverse_control' > $run_path/count.out"); 
	open (FILE_KMER,"$run_path/count.out");
	my $count=0;
	while (my $in_kmer_file=<FILE_KMER>) {
		my @split_kmer_file = split ('\t',$in_kmer_file);
		$count += $split_kmer_file[1];
	}
	close FILE_KMER;

	print OUT "---------------------------------------------------\n";
	print OUT "Spolygotyping\t";
	if ($count >=100) { ## threshold for pass/fail of the control
		print OUT "PASS\t";
	} else {
		print OUT "FAILED\t";
	}
	print OUT "Control Count = $count\n\n";
	my $count_spacers = 0;
	if ($count >=100) {
		## changed the lines 'cat -n $run_path/mers_count_25.out  ..' to 'cat $run_path/mers_count_25.out  ..' to correct for a glitch with agrep.

		## set the subroutine to submit all 43 kmer search in parallel..it hurts the system a bit but just for a few seconds.
		## if not ran in parrallel, this step takes 20+ minutes
		
		my $pm = Parallel::ForkManager->new(32);
                print "forked for 32";
                print "$pm";
		## run forward sequence spacers
		for my $job( 0..42 ) {
                    $pm->start and next;

			
		    if ($job == 2) {  ## because spacer 3 is shorter, agrep with the -w option was causing miss counting in some cases
                        system("cat $run_path/mers_count_25.out | grep $seq_spacers[$job] > $run_path/$job\_count.out");
                        if (-s "$run_path/$job\_count.out") {
                            ## skip file not empty
                        } else {
      			    #system("cat $run_path/mers_count_25.out | agrep -1 $seq_spacers[$job] > $run_path/$job\_count.out");
                            system("cat $run_path/mers_count_25.out > $run_path/x.txt; agrep -1 $seq_spacers[$job] $run_path/x.txt  > $run_path/$job\_count.out");
                        }
		    } else {
                        system("cat $run_path/mers_count_25.out | grep -w $seq_spacers[$job] > $run_path/$job\_count.out");
                        if (-s "$run_path/$job\_count.out") {
                            ## skip file not empty
                        } else {
                            #system("cat $run_path/mers_count_25.out | agrep -1 -w $seq_spacers[$job] > $run_path/$job\_count.out");
                            system("cat $run_path/mers_count_25.out > $run_path/x.txt; agrep -1 -w $seq_spacers[$job] $run_path/x.txt > $run_path/$job\_count.out");
                        }
		    }
		    $pm->finish;
		}

		$pm->wait_all_children;

		## run reverse-complement sequence spacers
		for my $job( 0..42) {
		    $pm->start and next;
		    my $reverse_spacer = reverse $seq_spacers[$job];
		    $reverse_spacer =~ tr/ATCG/TAGC/;
		    if ($job ==2) {
                        system("cat $run_path/mers_count_25.out | grep $reverse_spacer > $run_path/$job\_rcount.out");
                        if (-s "$run_path/$job\_rcount.out") {
                            ## skip file not empty
                        } else {
		            #system("cat $run_path/mers_count_25.out | agrep -1 $reverse_spacer > $run_path/$job\_rcount.out");
                            system("cat $run_path/mers_count_25.out >$run_path/x.txt; agrep -1 $reverse_spacer  $run_path/x.txt > $run_path/$job\_rcount.out");
                        }
	            } else {
                        system("cat $run_path/mers_count_25.out | grep -w $reverse_spacer > $run_path/$job\_rcount.out");
                        if (-s "$run_path/$job\_rcount.out") {
                            ## skip file not empty
                        } else {
		            #system("cat $run_path/mers_count_25.out | agrep -1 -w $reverse_spacer > $run_path/$job\_rcount.out");
                            system("cat $run_path/mers_count_25.out >$run_path/x.txt; agrep -1 -w $reverse_spacer $run_path/x.txt > $run_path/$job\_rcount.out");
                        }
		    }
		    $pm->finish;
		}
		$pm->wait_all_children;






		## count number of each spacers and generate spoligotype
		my $spacer=0;
		
		do {
			$count_spacers=0;
			open (FILE_COUNT,"$run_path/$spacer\_count.out");
			while (my $in_count = <FILE_COUNT>) {
				chomp $in_count;
				my @split_count_line = split ('\t',$in_count);
				$count_spacers += $split_count_line[1];
			}
			close FILE_COUNT;


			open (FILE_COUNT,"$run_path/$spacer\_rcount.out");
			while (my $in_count=<FILE_COUNT>) {
				chomp $in_count;
				my @split_count_line = split ('\t',$in_count);
				$count_spacers += $split_count_line[1];
			}
			close FILE_COUNT;


			## this step was to normalize the counts for the overlapping kmers
			$count_spacers = $count_spacers / (25 - length($seq_spacers[$spacer]) +1);

			if ($count_spacers >=10) {
				$stringx .= "1";
				$triplet .= "1";
			} else {
				$triplet .= "0";
				$stringx .= "0";
			}
			if (length($triplet) ==3) {
				$octa_num .= $octa{$triplet};
				$triplet="";
			}
			++$spacer;
		}
		until $spacer == 43;
				
		$octa_num .= $triplet;
                print "\n\n$octa_num\n\n";

		if ($spoligo{$octa_num}) {
			print OUT "CDC code :\t$stringx\t$octa_num\t$spoligo{$octa_num}\n";
		} else {
			print OUT "CDC code :\t$stringx\t$octa_num\tUnknown type\n";
		}

		if ($ny_spoligo{$octa_num}) {
			print OUT "NY code :\t$ny_spoligo{$octa_num}\n\n";
		} else {
			print OUT "NY :\t$stringx\t$octa_num\tUnknown type\n\n";
		} 

	} else {  ##failed spolygotype
		print OUT "\nFailed control In-silico spoligotyping (n = $count_spacers)\n\n"; 
	}


}

