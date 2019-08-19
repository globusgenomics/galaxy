#!/usr/bin/perl
use Data::Dumper;
use strict;
use Parallel::ForkManager;  #Forkmanager module
use IPC::Run 'run'; ##run [ "command", "arguments", "here" ], ">", \my $stdout;
use Benchmark;
use Data::Printer;

=for comments
This script is used to quickly screened all the consensus sequences for matches bellow
a certain SNP threshold. It does so by using the 'hamming' subrouting in parallel and 
print any matches to stdout
=cut


my $njobs=100; 
my $file1 = $ARGV[0]; #seed sequence to compare
my $run_directory = $ARGV[1];
my $results;


#load in seed consensus sequence
open (FILE,$file1);
my @in=<FILE>;
close FILE;
my $joined = join ('',@in);
my @tmp = split ('\n',$joined);
my $prefix1 = $tmp[0];
$prefix1 =~ s/\>//g;
shift @tmp;
my $query1 = join ('',@tmp);


# block to colect the data coming out of the parallel loop
my @collect_processed_data;
my $pm = Parallel::ForkManager->new( $njobs ); 
$pm->run_on_finish( sub {
    my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
	push @collect_processed_data, $data_structure_reference
      if defined($data_structure_reference);
    #$overall->{$pid} = $data_structure_reference;
});

#print "xx $run_directory/Consensus/*.fa\n";
#my @glob = glob ("$run_directory/Consensus/*.fa"); 
my @glob = glob ("/mnt/galaxyIndices/genomes/Mtuberculosis/Consensus/*.fa");
foreach my $file2(@glob) {	## iterate through all file in the query directory

	$pm->start and next; ## start jobs and go the next job	
	#run ["$script","$file1","$file2"],">", \my $results;
	#print $results;
	#$list .= $results;
	
	open (FILE,$file2);
	@in=<FILE>;
	close FILE;
	$joined = join ('',@in);
	@tmp = split ('\n',$joined);
	my $prefix2 = $tmp[0];
	$prefix2 =~ s/\>//g;
	shift @tmp;
	my $query2 = join ('',@tmp);
	
	my $different_positions=hd($query1,$query2);
	my $shared_N_positions=hdn($query1,$query2);
	my $n1=n($query1);
	my $n2=n($query2);
	my $SNP_differences = $different_positions - (($n1-$shared_N_positions) + ($n2-$shared_N_positions));
	
	my $results;
	if ($n2 < 1800000) {
		$results = "$prefix2\t$SNP_differences\n";	
	}
	$pm->finish (0,\$results );	## signal indicating that a job has finished

}
$pm->wait_all_children;   ## wait until all jobs are done to proceed furthe
#$Data::Dumper::Indent = 1;
#$Data::Dumper::Sortkeys = 1;
#print Dumper($overall);
# print @collect_processed_data;

#print "\@collect_processed_data has ", scalar @collect_processed_data, "elements\n";

## put results into hash
my %hash_results;
foreach my $ref (@collect_processed_data) {
	chomp $$ref;
	my @split_results = split ('\t',$$ref);
	if ($prefix1 !~ $split_results[0]) {
    	$hash_results{$split_results[0]} = $split_results[1];
    }
}
foreach my $dataset_name (sort { $hash_results{$a} <=> $hash_results{$b} } keys %hash_results) {
	print "$dataset_name\t$hash_results{$dataset_name}\n";
}
#my $t1 = Benchmark->new;
#my $td = timediff($t1, $t0);
#print  "\n\nruntime:",timestr($td),"\n";
exit;

sub hd {
    return ($_[0] ^ $_[1]) =~ tr/\001-\255//;
}
sub hdn {
    return ($_[0] & $_[1]) =~ tr/N//;
}

sub n {
    return ($_[0]) =~ tr/N//;
}
sub valid {
    return ($_[0] & $_[1]) =~ tr/ATCG\-//;
}
		






