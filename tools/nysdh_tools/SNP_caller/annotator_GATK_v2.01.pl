#!/usr/bin/perl

use strict;
use Parallel::ForkManager;
use IPC::Run 'run'; ##run [ "command", "arguments", "here" ], ">", \my $stdout;
###Changelog Version 0.2###
## Modify script to list out the resistance associated mutations for each high confidence loci individually.  Previous script was listed mutation based on the type of antimicrobial resistance, bypassing any loci without mutations
##  If a loci has no resistance associated mutation, listed as "no mutation detected"
##  also added a status and failed position if any for each assessed loci.

###Changelog Version 0.3###
## Change all final section of report to be CLIMS friendly

###Changelog Version 0.4###
## Final Annotation now only say INSERTION OR DELETION

###Changelog Version 1.0### 
## change Kanamycin to Kanamycin/Amikacin
## add ethA as list of HC target to detect any frameshift mutations
## add Stop codon list as resistant for specified loci
## add annotation for neutral (not associated with resistance) mutation from list Neutral.lst


###Changelog Version 2.0###
## mutation in the final reports will now be separated with ';' instead of ',' at @split_final_annotation_line
## add hashes for HC mutation and neutral mutation for annotation purposed
## added annotation bloc for variant to add HC mutation, Neutral or Unknown TAG
## added a procedure so if a loci has an unknown mutation (and no HC mutation),
## resistance will be tagged as unknown
## added SNP comparisons between consensus sequences 'find_match_fork' script
## remove subrouting where flag Isoniazid as passed even in presence of failed position at non-HC position 

###Changelog Version 2.01###
### Fixed various little bugs that were uncovered during testing
### Change script so only 1st line drugs are reported as Unknown for the final reports

use List::MoreUtils qw(uniq);
my ($day, $month, $year) = (localtime)[3,4,5];
my $date = sprintf("%02d/%02d/%04d",$month+1, $day, $year+1900);
my $path_to_GATK_INDELS = $ARGV[0];
my $path_to_GATK_SNP = $ARGV[1];
my $path=$ARGV[2];
my $prefix=$ARGV[3];
my $run_dir = $ARGV[4];
my $commandline = $0 . " ". (join " ", @ARGV);
my $loci = "/opt/galaxy/tools/nysdh_tools/SNP_caller/Resistance_loci_v3.txt";  ## resistance loci list	
my $hq = "/opt/galaxy/tools/nysdh_tools/SNP_caller/High_confidence_list_v3.txt"; ## High confidence mutation list 
my $neutral_mutations_list = "/opt/galaxy/tools/nysdh_tools/SNP_caller/Neutral_v2.lst";
my $hamming_fork = "/opt/galaxy/tools/nysdh_tools/SNP_caller/find_match_fork_V1.pl";

## list of high confidence loci for resistance
my @loci_list=('rpoB','katG','oxyR-ahpC promoter region','inhA','mabA-inhA promoter region','mabA','pncA','pncA promoter region','embB','embC-embA promoter region','rpsL','rrs','eis promoter region','gyrA','gyrB','ethA');
my %hash_of_loci = map { $_ => 1 } @loci_list;

print "------------------------------------------------------------------\n";
print "SNP discovery step:\n$commandline\nLoci = $loci, High confidence mutation list = $hq\n\n";

##Antibio with specific position to screen for
my @antimicrobial = ('Ethambutol','Fluoroquinolones','Isoniazid','Pyrazinamide','Rifampin','Streptomycin','Kanamycin/Amikacin','Kanamycin','Ethionamid');  

open (FILE,$loci);  ## list of all loci of interrest, put in array
my @tmp=<FILE>;
close FILE;
my $join_loci_file = join ('',@tmp);
my @loci = split ('\n',$join_loci_file); ## array with loci list

## list of the high confidence mutation, put in an array
open (FILE,$hq);  
@tmp=<FILE>;
close FILE;
my $join_hq_list = join ('',@tmp);
my @high_confidence_mutations = split ('\n',$join_hq_list);
shift @high_confidence_mutations; ## remove header
my %high_confidence_mutations_list=();
my %hash_of_HC;
my %hash_of_all_positions_changes;


foreach my $in_HC_list(@high_confidence_mutations) {
	my @split_HC_line = split ('\t',$in_HC_list);
	$high_confidence_mutations_list{$split_HC_line[0]} .= $in_HC_list;  # consolidate by antibiotic in the HC list in a hash
 
	if ($split_HC_line[8] =~ 'Non\-coding') {
		my @split_mutant_NT = split ('\,',$split_HC_line[5]);
		foreach my $individual_mutation(@split_mutant_NT) {	
			#$split_HC_line[1] =~ s/rrs1/rrs/g;
			#$split_HC_line[1] =~ s/rrs2/rrs/g;
			$hash_of_HC{$split_HC_line[1].$split_HC_line[3]."$split_HC_line[4] -> $individual_mutation"} = 1;
		}
	} else {
		my @split_mutant_AA = split ('\,',$split_HC_line[7]);
		foreach my $individual_mutation(@split_mutant_AA) {	
			$hash_of_HC{$split_HC_line[1].$split_HC_line[3]."$split_HC_line[6] -> $individual_mutation"} = 1;
		}
	}
	
}

### get indels position and stats from VCF file ####
### INDELS are called separatly for better accuracy #####

my %position_indel=();
my %heterogenecity;
my @split_line;
open (FILE,$path_to_GATK_INDELS);
while (my $in_VCF_indels=<FILE>) {  ## go through indels VCF and filters
	if ($in_VCF_indels !~ '\#') {	## if non header line	
		@split_line = split ('\t',$in_VCF_indels);
		
		### filter low qual indels ###
		my @split_at_DP = split ('DP\=',$in_VCF_indels); ## depth field
		my @split_downstream_of_DP = split ('\;',$split_at_DP[1]);
		my $dp=$split_downstream_of_DP[0];
		my @split_at_MQ = split ('MQ\=',$in_VCF_indels); ## mapping quality field
		my @split_downstream_of_MQ = split ('\;',$split_at_MQ[1]);
		my $mq=$split_downstream_of_MQ[0];
		my @split_at_QD = split ('QD\=',$in_VCF_indels); ## quality per depth field
		my @split_downstream_of_MQ = split ('\;',$split_at_QD[1]);
		my $qd=$split_downstream_of_MQ[0];
		my @split_at_FS = split ('FS\=',$in_VCF_indels); ## Fisher test field
		my @split_downstream_of_FS = split ('\;',$split_at_FS[1]);
		my $fs=$split_downstream_of_FS[0];
		my $ReadPosRankSum=0;
		if ($in_VCF_indels =~ 'ReadPosRankSum') {
			my @split_at_ReadPosRankSum = split ('ReadPosRankSum\=',$in_VCF_indels);
			my @split_downstream_of_ReadPosRankSum = split ('\;',$split_at_ReadPosRankSum[1]);
			$ReadPosRankSum=$split_downstream_of_ReadPosRankSum[0];
		}
		@split_line = split ('\t',$in_VCF_indels);
		
		if ($in_VCF_indels =~ '1\/1') {
			$heterogenecity{$split_line[1]} = "1/1 (Pure)";
		}
		if ($in_VCF_indels =~ '0\/1') {
			$heterogenecity{$split_line[1]} = "0/1 (Mixte)";

		}	
		if (($in_VCF_indels!~ 'LowQual') && ($dp >= 10) && ($mq >= 40) && ($qd >=2) && ($fs <=200) && ($ReadPosRankSum >= -20) && ($split_line[4] !~ '\,')) {			
			$position_indel{$split_line[1]} = $in_VCF_indels;  ##  If indel pass filter, stored as valid indels		
		}		
	}
}	
close FILE;

###########################################
## Parsing and combining both VCF files  ##
###########################################

my %hash_aa=(); # store codon position + Amino acids
my %hash_nt=(); # store nucleotide position + base
my %hash_indel=(); # store indels positions + types
my @frameshift=(); # store frameshift position + types
my %status=(); #
my @deleted_orf=(); ## keep list of any deleted orfs
my %antibio=(); #


open (OUT,">$path/tmp.vcf"); ## will generate a temporary VCF file with both combined VCF (indels + SNPs)
open (FILE,$path_to_GATK_SNP);

my $in_VCF_file=<FILE>;
print OUT $in_VCF_file;

if ($in_VCF_file !~ 'fileformat\=VCF') {
	die "Malformed or Missing VCF file\n";
}


## initialize values used in the VCF parsing

## previous genomic location
my $prev=0;

## number of sites assessed
my $sites=0;

##current genomic location
my $loc=1; 

## number of missing genomic locations
my $missing=0;

## number of location that failed QC thresholds
my $fail=0;

## number of location that passed QC threshold
my $pass=0;

## keep array of all genomic depth for stats calculation
my @depth;

## keep hash of wildtype positions
my %wt;

##hash of vcf lines per genome positions
my %position;
my $genome_size;


### get into hash list of neutral position to avoid flagging them as resistant
my %neutral_pos; ## hash to store neutral (not associated with resistance) position
open (TMP,$neutral_mutations_list);
while (my $length_WT_call=<TMP>) {
	chomp $length_WT_call;
	@split_line = split ('\t',$length_WT_call);
	$neutral_pos{$split_line[0].$split_line[1].$split_line[2]} = 1; ## neutral position stored in this hash
}
close TMP;


while (my $in_VCF_file=<FILE>) {  ## go through VCF, put in hash position{nt} and Status of the position
	if ($in_VCF_file =~ '\#\#contig\=') {		##contig=<ID=H37Rv,length=4411532>
		chomp $in_VCF_file;
		@split_line = split ('length\=',$in_VCF_file);
		$split_line[1] =~ s/\>//g;
		$genome_size = $split_line[1]; ## reference genome size		
	}				
	if ($in_VCF_file !~ '\#') {	## if not part of the header
		@split_line  = split ('\t',$in_VCF_file);	
		if (($in_VCF_file =~ '1\/1') && (!exists $heterogenecity{$split_line[1]})) { ## set homozygote of not already set by the indels call
			$heterogenecity{$split_line[1]} = "1/1 (Pure)";
		}
		if (($in_VCF_file =~ '0\/1') && (!exists $heterogenecity{$split_line[1]})) { ## set heterozygote of not already set by the indels call
			$heterogenecity{$split_line[1]} = "0/1 (Mixte)";
		}
		
		if ($position_indel{$split_line[1]}) { ## if the position has an indel, replace current position with the indel line 
			@split_line = split ('\t',$position_indel{$split_line[1]});
			$in_VCF_file=$position_indel{$split_line[1]};
		}
	
		## since GATK spits out two lines for indel non-indel position, substract 1 to all counters
		## to avoid overcounting
		if ($prev == $split_line[1]) {  
			--$loc;
			--$sites;
			--$missing;
			--$fail;
			--$pass;
		}
		if ($loc > $split_line[1]) {  ## if change in chromosome, reset genome position
			$loc=1;
		}
		

		if ($split_line[1] > $loc) { ##Fill up missing positions with Ns if not reported in the VCF file
			do {
				print OUT "$split_line[0]\t$loc\t$split_line[2]\tN\tN\t.\t0\t.\tAN=0;DP=0;MQ=0;MQ0=0\tNo mapping Information or Masked position\n";
				++$missing;
				push (@depth,"0");
				$status{$split_line[1]} = $split_line[9];
				$wt{$split_line[1]} = "N";
				$position{$split_line[1]} = "N";
				++$loc;
				++$sites
			}
			until $loc==$split_line[1];
		}


		if (($in_VCF_file =~ '\.\/\.') || ($split_line[3] =~ 'N')){
			$in_VCF_file = "$split_line[0]\t$split_line[1]\t$split_line[2]\tN\tN\t.\t0\t.\tAN=0;DP=0;MQ=0;MQ0=0\n";
			++$missing;
		}
	
	
		##retrieve VCF stats##
		my @split_at_DP = split ('DP\=',$in_VCF_file);
		my @split_downstream_DP = split ('\;',$split_at_DP[1]);
		my $dp=$split_downstream_DP[0];
		push (@depth,$dp);
		my @split_at_MQ = split ('MQ\=',$in_VCF_file);
		my @split_downstream_MQ = split ('\;',$split_at_MQ[1]);
		my $mq=$split_downstream_MQ[0];
		my $fs=60;
		if ($in_VCF_file =~ 'FS') {
			my @split_at_FS = split ('FS\=',$in_VCF_file);
			my @split_downstream_FS = split ('\;',$split_at_FS[1]);
			$fs=$split_downstream_FS[0];
		}
		my $ReadPosRankSum=0;
		if ($in_VCF_file =~ 'ReadPosRankSum') {
			my @split_at_ReadPosRankSum = split ('ReadPosRankSum\=',$in_VCF_file);
			my @split_downstream_ReadPosRankSum = split ('\;',$split_at_ReadPosRankSum[1]);
			my $ReadPosRankSum=$split_downstream_ReadPosRankSum[0];
		}
		my $qd=2;
		if ($in_VCF_file =~ 'QD') {
			my @split_at_QD = split ('QD\=',$in_VCF_file);
			my @split_downstream_QD = split ('\;',$split_at_QD[1]);
			$qd=$split_downstream_QD[0];
		}
		my $MQRankSum=0;
		if ($in_VCF_file =~ 'MQRankSum') {
			my @split_at_MQRankSum = split ('MQRankSum\=',$in_VCF_file);
			my @split_downstream_split_at_MQRankSum = split ('\;',$split_at_MQRankSum[1]);
			$MQRankSum=$split_downstream_split_at_MQRankSum[0];
		}	
	

		
		if (($in_VCF_file=~ 'LowQual') || ($dp < 10) || ($mq < 40) || ($fs > 60) || ($MQRankSum < -12.5) || ($ReadPosRankSum < -8) || ($qd < 2) && ((length($split_line[3]) == 1) || (length($split_line[4]) == 1))) {
			++$fail;
			$split_line[4] = "N";
			$in_VCF_file = join "\t",@split_line;
			$wt{$split_line[1]} = $split_line[3];  #get WT NT
			$position{$split_line[1]} = "N";
			if ($prev != $split_line[1]) {
				## skip
			}
		} else {
			++$pass;	
			$wt{$split_line[1]} = $split_line[3];  #get WT NT
			if ($split_line[4] =~ '\.') {		
				$position{$split_line[1]} = $split_line[3];  ##wild type position
			} else {
				$position{$split_line[1]} = $split_line[4];  ##variant position
			}		
		}	
		print OUT $in_VCF_file;
		++$loc;
		++$sites;
		$prev = $split_line[1];
	} else {
		print OUT $in_VCF_file;
	}
}
close FILE;
close OUT;



#### add indels



#print "Sites with no mapping information or masked : $missing\nSites below filters: $fail\nSites passing filters : $pass\n";
if ($sites < $genome_size) {  ## if VCF file too short compared to the reference genome
	die "VCF file too short\n";
}
###### end VCF parsing ######
print "\nAll Mutations in screened loci:\n\n";
##############################################################################################################################################
### SNP annotation part using all position stored in hash.  Annotation is done loci per loci and then in each loci, position per position. ###
### If loci is coding for a protein, annotation is done codon per codon (triplet of nucleotides)                                           ###
##############################################################################################################################################

my $string="";  ## placeholder for all loci mutations
my %hash_codon; ## hash of hash_codon codon returned by the codon sub-routine
my $amino_acid_wt;
my $amino_acid;
my %hash_aa_wt;
my %hash_aa;
my %hash_loci;
my $codon;
my $codon_wt;
my $length_variant_call;
my $length_WT_call;
my %hash_nt_wt;
my $start_position_3prime;		
my $stop_position_5prime;

foreach my $process_loci(@loci) { # go through Loci list				
	my @mutations=();
	my @failed=();
	$pass="PASS";
	chomp $process_loci;
	my @split_loci_line = split ('\t',$process_loci);	
	
	
	
	
	
	#### if positive strand and not intergenic loci ####
	
	
	if (($split_loci_line[2] =~ '\+') && ($split_loci_line[3] !~ 'intergenic|rrs')) {  
		my $start = $split_loci_line[0]; ## loci start position
		my $stop = $split_loci_line[1]; ## loci stop position
		my $codon_position=1; ##initialize codon position
		my $size = $stop - $start +1 ; ## length of loci
		## look for Indels ##
		my $pos=0;
		do {  # go position per position
			if ($position{$start} =~ 'N') {
				$pass="FAIL";
				push (@failed,$start);
			}
			my $real_loci_start = $start+1;
			
			### look for Insertions or Deletions
			$length_variant_call = length($position{$start});
			$length_WT_call = length($wt{$start});
			
			if ($length_variant_call > 1) { ## mean insertion in variant
				my @split_variant_based = split ('',$position{$start});
				shift @split_variant_based;
				my $size_insertion = @split_variant_based;
				my $join = join ('',@split_variant_based);
				if (($size_insertion % 3)) {		## check for frameshift mutation using modulo			
					$string .= "$real_loci_start\t$pos\tFrameshift Insertion\t+$join\t$split_loci_line[3]\t$split_loci_line[4]\t$split_loci_line[5]\t$heterogenecity{$start}\n";
					
					## if insertion in rpoB, keep as potential resistant
					if ($split_loci_line[4] =~ 'rpoB') {
						if (!exists $neutral_pos{$split_loci_line[4].$pos."+$join"}) { ##ignore if part of neutral mutation list
							push @frameshift,"$split_loci_line[5]\t$split_loci_line[4]\t$pos\t+$join\tINSERTION\t$heterogenecity{$start}";
							$hash_loci{$split_loci_line[4]} .= "separatorINSERTION";
							#$hash_loci{$split_loci_line[4]} .= "separator$pos+$join";								
							$antibio{$split_loci_line[5]} = 1;
						}
					}
				
				
				} else {			## inframe insertion					
					my $frag="";
					my $frag_start=0;
					do {
						my $codon = $split_variant_based[$frag_start].$split_variant_based[$frag_start+1].$split_variant_based[$frag_start+2];				
						my $amino_acid = &codon2aa($codon);
						$frag .= $amino_acid;
						$frag_start += 3;	
					}
					until $frag_start == $size_insertion;
					$string .= "$real_loci_start\t$pos\tIn-Frame Insertion\t+$join\t+$frag\t$split_loci_line[3]\t$split_loci_line[4]\t$split_loci_line[5]\t$heterogenecity{$start}\n";
					
					if ($split_loci_line[4] =~ 'rpoB') {
						if (!exists $neutral_pos{$split_loci_line[4].$pos."+$join"}) { ##ignore if part of neutral mutation list
							push @frameshift,"$split_loci_line[5]\t$split_loci_line[4]\t$pos\t+$join\t+$frag\tINSERTION\t$heterogenecity{$start}";	
							$hash_loci{$split_loci_line[4]} .= "separatorINSERTION";
							#$hash_loci{$split_loci_line[4]} .= "separator$pos+$join";								
							$antibio{$split_loci_line[5]} = 1;
						}
					}
				}	
			} elsif ($length_WT_call > 1) { ## mean deletion in variant
				my @split_WT_based = split ('',$wt{$start});
				shift @split_WT_based;								
				$size = @split_WT_based;
				my $join = join ('',@split_WT_based);
				
				if (($size % 3)) {						
					$string .= "$real_loci_start\t$pos\tFrameshift Deletion\t-$join\t\t$split_loci_line[3]\t$split_loci_line[4]\t$split_loci_line[5]\t$heterogenecity{$start}\n";	
					if ($split_loci_line[4] =~ 'rpoB') {
						if (!exists $neutral_pos{$split_loci_line[4].$pos."-$join"}) { ##ignore if part of neutral mutation list
							push @frameshift,"$split_loci_line[5]\t$split_loci_line[4]\t$pos\t-$join\tDELETION\t$heterogenecity{$start}";
							$hash_loci{$split_loci_line[4]} .= "separatorDELETION";	
							#$hash_loci{$split_loci_line[4]} .= "separator$pos-$join";
							$antibio{$split_loci_line[5]} = 1;	
						}	
					}	
				} else {
					my $frag="";
					my $frag_start=0;
					do {
						$codon = $split_WT_based[$frag_start].$split_WT_based[$frag_start+1].$split_WT_based[$frag_start+2];	
						$hash_codon{$codon} = "XXX"; ## to avoid codon translation error if empty
						$amino_acid = &codon2aa($codon);
						$frag .= $amino_acid;						
						$frag_start += 3;	
					}
					until $frag_start == $size;					
					$string .= "$real_loci_start\t$pos\tIn-Frame Deletion\t-$join\t-$frag\t$split_loci_line[3]\t$split_loci_line[4]\t$split_loci_line[5]\t$heterogenecity{$start}\n";
					
					if ($split_loci_line[4] =~ 'rpoB') {
						if (!exists $neutral_pos{$split_loci_line[4].$pos."-$join"}) { ##ignore if part of neutral mutation list
							push @frameshift,"$split_loci_line[5]\t$split_loci_line[4]\t$pos\t-$join\t-$frag\tDELETION\t$heterogenecity{$start}";
							$hash_loci{$split_loci_line[4]} .= "separatorDELETION";
							#$hash_loci{$split_loci_line[4]} .= "separator$pos-$join";
							$antibio{$split_loci_line[5]} = 1;
						}
					}
				}
						
			}
			

			++$start;
			++$pos;

		}
		until $start >= $stop + 1;

### look for Codon based mutations
		$start = $split_loci_line[0];
		do {	
			$codon = $position{$start}.$position{$start+1}.$position{$start+2};
			$codon_wt = $wt{$start}.$wt{$start+1}.$wt{$start+2};
			
			if ($position{$start+1} !~ $wt{$start+1}) {
				$heterogenecity{$start} = $heterogenecity{$start+1};
				} elsif ($position{$start+2} !~ $wt{$start+2}) {
				$heterogenecity{$start} = $heterogenecity{$start+2};				
			}
			
			
			
			$length_variant_call = length($codon);
			$length_WT_call = length($codon_wt);

			if (($length_WT_call == 3 )&& ($length_variant_call == 3)) {			
				$hash_codon{$codon} = "XXX";
				$amino_acid_wt = &codon2aa($codon_wt);
				$amino_acid = &codon2aa($codon);
				$hash_aa_wt{$split_loci_line[4].$codon_position} = $codon_wt;	
				$hash_aa{$split_loci_line[4].$codon_position} = $codon;
						
				if (($codon !~ $codon_wt) && ($codon !~ 'N')){
				
				if ($amino_acid_wt =~ $amino_acid) {
					$string .= "$start\t$codon_position\t$codon_wt -> $codon\t$amino_acid_wt -> $amino_acid\t$split_loci_line[3]\t$split_loci_line[4]\t$split_loci_line[5] (Silent)\t$heterogenecity{$start}\n";
				
				} else {			
					$string .= "$start\t$codon_position\t$codon_wt -> $codon\t$amino_acid_wt -> $amino_acid\t$split_loci_line[3]\t$split_loci_line[4]\t$split_loci_line[5]\t$heterogenecity{$start}\n";
				}
					$antibio{$split_loci_line[5]} = 1;
					
				}
			}
			++$codon_position;
			$start += 3;
		}
		until $start >= $stop + 1;
		

	####  End In-ORFs SNPs positive strand ############


	### look for mutation on positive strand in intergenic positions ######

	} elsif (($split_loci_line[2] =~ '\+') && ($split_loci_line[3] =~ 'intergenic')) {	
		my $start = $split_loci_line[0];
		my $stop = $split_loci_line[1];
		$loc=$start-$stop-1;
		my $real_loci_start = $start + 1;
		do {	
			if ($position{$start} =~ 'N') {
				$pass="FAIL";
				push (@failed,$start);
			}
			$hash_nt{$split_loci_line[4].$loc} = $position{$start};
			$length_variant_call = length($position{$start});
			$length_WT_call = length($wt{$start});
			$hash_nt_wt{$split_loci_line[4].$loc} = $wt{$start};
			if (($position{$start} !~ $wt{$start}) && ($position{$start} !~ 'N')) {			
				if ($length_variant_call > 1) {			## look for Indels ##				
					@tmp = split ('',$position{$start});
					shift @tmp;
					my $join = join ('',@tmp);
					$string .= "$real_loci_start\t$loc\tIntergenic Insertion\t+$join\t$split_loci_line[3]\t$split_loci_line[4]\t$split_loci_line[5]\t$heterogenecity{$start}\n";
					push (@mutations,"$real_loci_start\t$loc\tItergenic INSERTION\t+$join\t$split_loci_line[3]\t$split_loci_line[4]\t$split_loci_line[5]\t$heterogenecity{$start}\n");
					$antibio{$split_loci_line[5]} = 1;			
				} elsif ($length_WT_call > 1) {
					@tmp = split ('',$wt{$start});
					shift @tmp;
					my $join = join ('',@tmp);
				
				
				
					$string .= "$real_loci_start\t$loc\tIntergenic deletion\t-$join\t$split_loci_line[3]\t$split_loci_line[4]\t$split_loci_line[5]\t$heterogenecity{$start}\n";
					push (@mutations,"$real_loci_start\t$loc\tItergenic DELETION\t-$join\t$split_loci_line[3]\t$split_loci_line[4]\t$split_loci_line[5]\t$heterogenecity{$start}\n");
					$antibio{$split_loci_line[5]} = 1;
				} else {	

					$string .= "$start\t$loc\t$wt{$start} -> $position{$start}\t$split_loci_line[3]\t$split_loci_line[4]\t$split_loci_line[5]\t$heterogenecity{$start}\n";		
					push (@mutations,"$start\t$loc\t$wt{$start} -> $position{$start}\t$split_loci_line[3]\t$split_loci_line[4]\t$split_loci_line[5]\t$heterogenecity{$start}\n");
					$antibio{$split_loci_line[5]} = 1;
				}		
			}				

			++$start;
			++$loc;		
		}
		until $start >= $stop +1;
	
	
	
	
	### look for mutation on positive strand in ribosomal RNA positions ######
		
	} elsif (($split_loci_line[2] =~ '\+') && ($split_loci_line[3] =~ 'rrs')) {
		my $start = $split_loci_line[0];
		my $real_loci_start = $start + 1;
		my $stop = $split_loci_line[1];
		$loc=1;
		do {	
			if ($position{$start} =~ 'N') {
				$pass="FAIL";
				push (@failed,$start);
			}
			$hash_nt{$split_loci_line[4].$loc} = $position{$start};
		if (($position{$start} !~ $wt{$start}) && ($position{$start} !~ 'N')) {
			$length_variant_call = length($position{$start});
			$length_WT_call = length($wt{$start});
			if ($length_variant_call>1) {			## look for Indels ##			
				@tmp = split ('',$position{$start});
				shift @tmp;shift @tmp;
				my $join = join ('',@tmp);
				$string .= "$real_loci_start\t$loc\trrs Insertion\t+$join\t$split_loci_line[3]\t$split_loci_line[4]\t$split_loci_line[5]\t$heterogenecity{$start}\n";
				push (@mutations,"$real_loci_start\t$loc\trrs INSERTION\t+$join\t$split_loci_line[3]\t$split_loci_line[4]\t$split_loci_line[5]\t$heterogenecity{$start}\n");	
				$antibio{$split_loci_line[5]} = 1;			
			} elsif ($length_WT_call>1) {
				$string .= "$real_loci_start\t$loc\trrs deletion\t$position{$start}\t$split_loci_line[3]\t$split_loci_line[4]\t$split_loci_line[5]\t$heterogenecity{$start}\n";	
				push (@mutations,"$real_loci_start\t$loc\trrs DELETION\t$position{$start}\t$split_loci_line[3]\t$split_loci_line[4]\t$split_loci_line[5]\t$heterogenecity{$start}\n");
				$antibio{$split_loci_line[5]} = 1;
			} else {	
				if (($loc == 1400) && ($split_loci_line[3] =~ 'rrs')) {
					$split_loci_line[5] = "Kanamycin/Amikacin";
				}
				$string .= "$start\t$loc\t$wt{$start} -> $position{$start}\t$split_loci_line[3]\t$split_loci_line[4]\t$split_loci_line[5]\t$heterogenecity{$start}\n";	
				$antibio{$split_loci_line[5]} = 1;	
			}	
		}	

			++$start;
			++$loc;
		#print "$loc $wt{$start} -> $position{$start}\n";
	}
	until $start >= $stop +1;	

	
	
	
	####################################################################
	##### mutation in Negative strand  codong regions ##################
	####################################################################
	
		
	} elsif (($split_loci_line[2] =~ '\-') && ($split_loci_line[3] !~ 'intergenic|rrs')) {
		my $start = $split_loci_line[1];
		my $stop = $split_loci_line[0];
		my $codon_position=1;
		my $size = $start - $stop +1;
		$loc=1;
			do {
				if ($position{$start} =~ 'N') {
				$pass="FAIL";
				push (@failed,$start);
				}
			
			$length_variant_call = length($position{$start});
			$length_WT_call = length($wt{$start});
			
			
			### look for Insertions or Deletions
			if ($length_variant_call > 1) {
				@tmp = split ('',$position{$start});
				shift @tmp;
				$size = @tmp;
				my $join = join ('',@tmp);
				$join = reverse $join;
				$join =~ tr/ATCG/TAGC/;
				@tmp = split ('',$join);
				$start_position_3prime = $start + $size;		
				$stop_position_5prime = $start_position_3prime - 1;									
				if (($size % 3)) {					
					$string .= "$start_position_3prime\t$loc\tFrameshift Insertion\t+$join\t$split_loci_line[3]\t$split_loci_line[4]\t$split_loci_line[5]\t$heterogenecity{$stop_position_5prime}\n";					                                                                         
					if ($split_loci_line[4] !~ 'gidB|iniB|embR') {
						if (!exists $neutral_pos{$split_loci_line[4].$loc."+$join"}) { ##ignore if part of neutral mutation list
							push @frameshift,"$split_loci_line[5]\t$split_loci_line[4]\t$loc\t+$join\tINSERTION\t$heterogenecity{$stop_position_5prime}";
							$hash_loci{$split_loci_line[4]} .= "separatorINSERTION";
							#$hash_loci{$split_loci_line[4]} .= "separator$loc+$join";
							$antibio{$split_loci_line[5]} = 1;
						}
					}
				} else {								
					my $frag="";
					my $frag_start=0;
					do {
						$codon = $tmp[$frag_start].$tmp[$frag_start+1].$tmp[$frag_start+2];	
						$hash_codon{$codon} = "XXX";					
						$amino_acid = &codon2aa($codon);
						$frag .= $amino_acid;
						$frag_start += 3;	
					}
					until $frag_start == $size;
					$string .= "$start_position_3prime\t$loc\tIn-Frame Insertion\t+$join\t+$frag\t$split_loci_line[3]\t$split_loci_line[4]\t$split_loci_line[5]\t$heterogenecity{$stop_position_5prime}\n";
					if ($split_loci_line[4] =~ 'katG|pncA') {
						if (!exists $neutral_pos{$split_loci_line[4].$loc."+$join"}) { ##ignore if part of neutral mutation list
							push @frameshift,"$split_loci_line[5]\t$split_loci_line[4]\t$loc\t+$join\t+$frag\tINSERTION\t$heterogenecity{$stop_position_5prime}";
							$hash_loci{$split_loci_line[4]} .= "separatorINSERTION";
							#$hash_loci{$split_loci_line[4]} .= "separator$loc+$join";
							$antibio{$split_loci_line[5]} = 1;
						}
					}
					
				}	
			} elsif ($length_WT_call > 1) {				
				@tmp = split ('',$wt{$start});

				shift @tmp;				
				my $size = @tmp;
				my $join = join ('',@tmp);
				$join = reverse $join;
				$join =~ tr/ATCG/TAGC/;
				@tmp = split ('',$join);
				$start_position_3prime = $start + $size;				
				if ($size > 0) {
					if (($size % 3)) {
						$string .= "$start_position_3prime\t$loc\tFrameshift Deletion\t-$join\t\t$split_loci_line[3]\t$split_loci_line[4]\t$split_loci_line[5]\t$heterogenecity{$stop_position_5prime}\n";							
						if ($split_loci_line[4] =~ 'katG|pncA|ethA') {
							if (!exists $neutral_pos{$split_loci_line[4].$loc."-$join"}) { ##ignore if part of neutral mutation list
								push @frameshift,"$split_loci_line[5]\t$split_loci_line[4]\t$loc\t-$join\tDELETION\t$heterogenecity{$stop_position_5prime}";
								$hash_loci{$split_loci_line[4]} .= "separatorDELETION";	
								#$hash_loci{$split_loci_line[4]} .= "separator$loc-$join";
								$antibio{$split_loci_line[5]} = 1;		
							}
						}			
					} else {
						my $frag="";
						my $frag_start=0;
						do {
							$codon = $tmp[$frag_start].$tmp[$frag_start+1].$tmp[$frag_start+2];	
							$amino_acid = &codon2aa($codon);
							$frag .= $amino_acid;						
							$frag_start += 3;	
						}
						until $frag_start == $size;
						$start_position_3prime	= $start + $size;
						$string .= "$start_position_3prime\t$loc\tIn-Frame Deletion\t-$join\t-$frag\t$split_loci_line[3]\t$split_loci_line[4]\t$split_loci_line[5]\t$heterogenecity{$stop_position_5prime}\n";
						if ($split_loci_line[4] =~ 'katG|pncA') {
							if (!exists $neutral_pos{$split_loci_line[4].$loc."-$join"}) { ##ignore if part of neutral mutation list
								$hash_indel{$split_loci_line[4].$loc} = "-$join";	
								push @frameshift,"$split_loci_line[5]\t$split_loci_line[4]\t$loc\t-$join\t$frag\tDELETION\t$heterogenecity{$stop_position_5prime}";	
								$hash_loci{$split_loci_line[4]} .= "separatorDELETION";	
								#$hash_loci{$split_loci_line[4]} .= "separator$loc-$join";
							}
						}
					}
				}		
			}
			

			--$start;
			++$loc;
			

		}
		until $start <= $stop - 1;
				
		#### In ORFs mutations #####
		my $start = $split_loci_line[1];
		my $stop = $split_loci_line[0];
		my $codon_position=1;
		$length_variant_call = split ('',$position{$start});
		$length_WT_call = split ('',$wt{$start});
		do {
			$codon_wt = $wt{$start}.$wt{$start-1}.$wt{$start-2};
			$codon = $position{$start}.$position{$start-1}.$position{$start-2};
			
			if ($position{$start-1} !~ $wt{$start-1}) {
				$heterogenecity{$start} = $heterogenecity{$start-1};
			} elsif ($position{$start-2} !~ $wt{$start-2}) {
				$heterogenecity{$start} = $heterogenecity{$start-2};				
			}			
			$length_variant_call = length($codon);
			$length_WT_call = length($codon_wt);
			if (($length_variant_call==3)&& ($length_WT_call==3)) {
			
			$codon =~ tr/ATCG/TAGC/;
			$codon_wt =~ tr/ATCG/TAGC/;
			$hash_codon{$codon} = "XXX";
			$amino_acid_wt = &codon2aa($codon_wt);
			$amino_acid = &codon2aa($codon);
			$hash_aa{$split_loci_line[4].$codon_position} = $codon;
			$hash_aa_wt{$split_loci_line[4].$codon_position} = $codon_wt;	
			if (($codon !~ $codon_wt) && ($codon !~ 'N')){
				if ($amino_acid_wt =~ $amino_acid) {
					$string .= "$start\t$codon_position\t$codon_wt -> $codon\t$amino_acid_wt -> $amino_acid\t$split_loci_line[3]\t$split_loci_line[4]\t$split_loci_line[5] (Silent)\t$heterogenecity{$start}\n";
				} else {
					$string .= "$start\t$codon_position\t$codon_wt -> $codon\t$amino_acid_wt -> $amino_acid\t$split_loci_line[3]\t$split_loci_line[4]\t$split_loci_line[5]\t$heterogenecity{$start}\n";
				}
				$antibio{$split_loci_line[5]} = 1;
				}
			}
			++$codon_position;
			$start -= 3;
		}
		until $start <= $stop - 1;

    ### annotate snp in intergenic region on negative strand ###
	} elsif (($split_loci_line[2] =~ '\-') && ($split_loci_line[3] =~ 'intergenic')) {
	
		my $start = $split_loci_line[1];
		my $stop = $split_loci_line[0];
		$loc = $stop-$start-1;
		do {
			if ($position{$start} =~ 'N') {
				$pass="FAIL";
				push (@failed,$start);
			}
			my $base = $position{$start};
			$base =~ tr/ATCG/TAGC/;
			$hash_nt{$split_loci_line[4].$loc} = $base;		
			$wt{$start} =~ tr/ATCG/TAGC/;
			$hash_nt_wt{$split_loci_line[4].$loc} = $wt{$start};	
			$wt{$start} =~ tr/ATCG/TAGC/;
			if (($position{$start} !~ $wt{$start}) && ($position{$start} !~ 'N')) {	#variant position	
				$position{$start} =~ tr/ATCG/TAGC/;
				$wt{$start} =~ tr/ATCG/TAGC/;
				$length_variant_call = split ('',$position{$start});
				$length_WT_call = split ('',$wt{$start});
				$hash_nt{$split_loci_line[4].$loc} = $position{$start};
				$hash_nt_wt{$split_loci_line[4].$loc} = $wt{$start};
			if ($length_variant_call>1) {			## look for Indels ##
				@tmp = split ('',$position{$start});
				shift @tmp;
				my $join = join ('',@tmp);
				$join = reverse $join;
				$join =~ tr/ATCG/TAGC/;
				my $size = @tmp;
				$start_position_3prime= $start + $size;
				$string .= "$start_position_3prime\t$loc\tIntergenic Insertion\t+$join\t$split_loci_line[3]\t$split_loci_line[4]\t$split_loci_line[5]\t$heterogenecity{$start}\n";
				$antibio{$split_loci_line[5]} = 1;
			} elsif ($length_WT_call>1) { ## look for deletion
				@tmp = split ('',$wt{$start});
				shift @tmp;
				my $join = join ('',@tmp);
				$join = reverse $join;
				$join =~ tr/ATCG/TAGC/;
				my $size = @tmp;
				$start_position_3prime= $start + $size;			
				$string .= "$start_position_3prime\t$loc\tIntergenic deletion\t-$join\t$split_loci_line[3]\t$split_loci_line[4]\t$split_loci_line[5]\t$heterogenecity{$start}\n";	
				$antibio{$split_loci_line[5]} = 1;
			} else {			
				$string .= "$start\t$loc\t$wt{$start} -> $position{$start}\t$split_loci_line[3]\t$split_loci_line[4]\t$split_loci_line[5]\t$heterogenecity{$start}\n";
				$antibio{$split_loci_line[5]} = 1;


			}
		}	
		
		
		++$loc;
		--$start;
		
	}
	until $start <= $stop - 1;
	
	} 

} 



##### process $string to annotate phenotypically neutral mutations
=to comment
this bloc go though the list of variant previously annotated and check present in the list
of neutral positions or list of HC mutation and add appropriate tag.
If not neutral or HC, annotated as unknown
=cut


my %unknown_per_antibiotics;
my @list = split ('\n',$string);
my %hash_for_variant_ID;

foreach my $list_of_variant(@list) {

	if ($list_of_variant !~ 'Silent') {
		@split_line = split ('\t',$list_of_variant);
		if ($split_line[3] =~ 'intergenic|rrs') {
			$split_line[6] = $split_line[5];
			$split_line[5] = $split_line[4];
			$split_line[3] = $split_line[2];
		}

		
		if ($split_line[2] =~ 'Insertion|Deletion') {
			$split_line[5] = $split_line[6];			
			$split_line[6] = $split_line[7];

		}
		
		
		if (exists $neutral_pos{$split_line[5].$split_line[1].$split_line[3]}) {
			print "$list_of_variant (Neutral)\n";
		} elsif (exists $hash_of_HC{$split_line[5].$split_line[1].$split_line[3]}) {
			print "$list_of_variant (HC mutation)\n";
	
		} else {
			print "$list_of_variant (Unknown)\n";
			$split_line[3] =~ s/ \-> /$split_line[1]/g;			
			if (exists $hash_of_loci{$split_line[5]}) {
				$unknown_per_antibiotics{$split_line[6]} = "Unknown";
			}

		}
	} else {
		print "$list_of_variant\n";
	}
} 

print "\n\n------------------------------------------------------------------\nSNP-Base identification:\n\n";
######################################################
## 	SNP based identification 	###

### lineage ID  ####
my @lineage;
print "Lineage:\n";
if ($hash_aa{"gidB110"} =~ 'GTT') {
	print "Lineage 1 (Indo-Oceanic)\n";
	push (@lineage,"Lineage 1 (Indo-Oceanic)");
}
if ($hash_aa{"rpsA212"} =~ 'CGC') {
	print "Lineage 2 (Beijing)\n";
	push (@lineage,"Lineage 2 (Beijing)");
}

if ($position{2726105} =~ 'A') {
	print "Lineage 3 (Central-Asian)\n";
	push (@lineage,"Lineage 3 (Central-Asian)");
}
if ($hash_aa{"katG463"} =~ 'CGG') {
	print "Lineage 4 (Euro-American)\n";
	push (@lineage,"Lineage 4 (Euro-American)");
}
if ($hash_aa{"ethA124"} =~ 'GAC') {
	print "Lineage 5 (West African 1)\n";
	push (@lineage,"Lineage 5 (West African 1)");
}
if ($hash_aa{"inhA78"} =~ 'GCG') {
	print "Lineage 6 (West African 2)\n";
	push (@lineage,"Lineage 6 (West African 2)");
}
print "\n";
print "Species:\n";





if (($hash_aa{"gyrB403"} =~ 'GCG') && ($hash_aa{"katG203"} =~ 'ACC')) {
	print "Mycobacterium tuberculosis\n";
}
if (($position{1673338} =~ 'A') && ($hash_aa{"ethA124"} =~ 'GAC')) {
	print "Mycobacterium africanum\n";
}
if (($hash_aa{"inhA78"} =~ 'GCG') && ($hash_aa{"atpE69"} =~ 'GCT')) {
	print "Mycobacterium africanum\n";
}
if (($hash_aa{"inhA107"} =~ 'TCG') && ($position{1473094} =~ 'C')) {
	print "Mycobacterium pinnipedii\n";
}
if (($hash_aa{"gyrB144"} =~ 'TAT') && ($position{1473079} =~ 'A')) {
	print "Mycobacterium microti\n";
}
if (($hash_aa{"gyrB171"} =~ 'GTA') && ($hash_aa{"gyrB356"} =~ 'GCG')) {
	print "Mycobacterium caprae\n";
}
if (($hash_aa{"pncA57"} =~ 'GAC') && ($hash_aa{"furA43"} =~ 'GTC')) {
	print "Mycobacterium bovis-BCG\n";
}
if (($hash_aa{"pncA57"} =~ 'GAC') && ($hash_aa{"furA43"} =~ 'GCC')) {
	print "Mycobacterium bovis\n";
}
########## load lumpy_sv output and evaluate #######
########## To have a genuine Deletion, lumpy and read coverage needs to agree 
########## (interval from lumpy must have < 1 avg dp)
########## We are only looking for large deletions. Large insertion are unlikely 
########## to happens (never been described associated with resistance)
########## and will be very difficult to pick up with confidence using lumpy
open (LUMPY,"$path/sample_svn.out");
my @sub; ## array to store depth for loci location
my $avg_dp; ## variable for average read depth at loci
my $del_length; ## Length of deletion
my $percent_deletion; ## percent of loci deleted
while (my $lumpy_line = <LUMPY>) {
	if ($lumpy_line =~ 'DELETION') {
		my @split_lumpy_line = split ('\t',$lumpy_line);
		my @split_for_location = split ('\:',$split_lumpy_line[13]);
		my $deletion_end = $split_for_location[3];
		my @split_location_start = split ('\;',$split_for_location[2]);
		my $deletion_start = $split_location_start[0];
		foreach my $process_loci(@loci) { # go through Loci list
			chomp $process_loci;
			my @split_loci_line = split ('\t',$process_loci);
			
			if (($split_loci_line[0] >= $deletion_start) && ($split_loci_line[1] <= $deletion_end)) { ##full loci deletion
				my $length_loci = $split_loci_line[1] - $split_loci_line[0] + 1;
				
				@sub  = @depth[($split_loci_line[0]-1)..($split_loci_line[1]-1)];
				$avg_dp = average(@sub);
				if ($avg_dp < 1) {
					push @deleted_orf,"$split_loci_line[4] Missing from the genome ($length_loci nt [100%] of the locus missing)\t-->\t$split_loci_line[5]";
					push @frameshift,"$split_loci_line[5]\t$split_loci_line[4]";
					if ($split_loci_line[4] =~ 'katG|rpoB|pncA|ethA') {
						$hash_loci{$split_loci_line[4]} .= "separatorDELETION";			
					}
				}
				
				$antibio{$split_loci_line[5]} = 1;
			} elsif (($split_loci_line[0] >= $deletion_start) && ($split_loci_line[1] > $deletion_end) && ($split_loci_line[0] < $deletion_end)) { ##partial upstream deletion
				my $length_loci = $split_loci_line[1] - $split_loci_line[0] + 1;
				$del_length = $deletion_end - $split_loci_line[0] + 1;
				$percent_deletion = sprintf ("%.2f",100*($del_length/$length_loci));
				
				@sub  = @depth[($split_loci_line[0]-1)..($deletion_end-1)];
				$avg_dp = average(@sub);
				if ($avg_dp < 1) {
					push @deleted_orf,"$split_loci_line[4] presence of a major deletion ($del_length nt [$percent_deletion%] of the locus missing)\t-->\t$split_loci_line[5]";
					push @frameshift,"$split_loci_line[5]\t$split_loci_line[4]";
					if ($split_loci_line[4] =~ 'katG|rpoB|pncA|ethA') {
						$hash_loci{$split_loci_line[4]} .= "separatorDELETION";			
					}
				}
				$antibio{$split_loci_line[5]} = 1;
			} elsif (($split_loci_line[1] <= $deletion_end) && ($split_loci_line[0] < $deletion_start) && ($split_loci_line[1] > $deletion_start)) { ##partial downstream deletion
				my $length_loci = $split_loci_line[1] - $split_loci_line[0] + 1;
				$del_length = $split_loci_line[1] - $deletion_start + 1;
				$percent_deletion = sprintf ("%.2f",100*($del_length/$length_loci));
				@sub  = @depth[(($deletion_start-1)..($split_loci_line[1]-1))];
				$avg_dp = average(@sub);				
				if ($avg_dp < 1) {
					push @deleted_orf,"$split_loci_line[4] presence of a major deletion ($del_length nt [$percent_deletion%] of the locus missing)\t-->\t$split_loci_line[5]";
					push @frameshift,"$split_loci_line[5]\t$split_loci_line[4]";
					if ($split_loci_line[4] =~ 'katG|rpoB|pncA|ethA') {
						$hash_loci{$split_loci_line[4]} .= "separatorDELETION";		
					}
				}
				$antibio{$split_loci_line[5]} = 1;
			} elsif (($split_loci_line[1] > $deletion_end) && ($split_loci_line[0] < $deletion_start)) { ##in-loci deletion
				my $length_loci = $split_loci_line[1] - $split_loci_line[0] + 1;
				$del_length = $deletion_end - $deletion_start + 1;
				$percent_deletion = sprintf ("%.2f",100*($del_length/$length_loci));
				@sub  = @depth[(($deletion_start-1)..($deletion_end-1))];
				$avg_dp = average(@sub);				
				if ($avg_dp < 1) {								
					push @deleted_orf,"$split_loci_line[4] presence of a major deletion ($del_length nt [$percent_deletion%] of the locus missing)\t-->\t$split_loci_line[5]";
					push @frameshift,"$split_loci_line[5]\t$split_loci_line[4]";
					if ($split_loci_line[4] =~ 'katG|rpoB|pncA|ethA') {
						$hash_loci{$split_loci_line[4]} .= "separatorDELETION";
					}
				}
				$antibio{$split_loci_line[5]} = 1;
			}
		}
	}
}
close LUMPY;			
####################################################

print "\n";
print join "\n",@deleted_orf;

########### screen for High confidence mutations #############
my $amino_acid_var; ## AA at variant site
my $amino_acid_wt;  ## Wild type AA
my @HC_aa;  ## high_confidence AA
my $screen; ## antimicrobial to screen for mutations
my $status;
my %hash_failed; ## hash of failed positions
my %resis;





foreach $screen(@antimicrobial) {
	#print "\n$screen:\n";
	$pass=0; ## set to default status 0 = pass
	$status="ND"; ## Set status default to ND.
	my @split_mutation_list_by_antibiotics = split ($screen,$high_confidence_mutations_list{$screen});
	shift @split_mutation_list_by_antibiotics;

	
	foreach my $position_mutation(@split_mutation_list_by_antibiotics) { 
		my @split_HC_mutation_line = split ('\t',$position_mutation);
		
		#### look for any stop codon in certain loci flagged 'Any_C'
				
		if ($split_HC_mutation_line[8] =~ 'Stop_codon') {
			foreach my $key (keys %hash_aa) {		 
				if ($key =~ $split_HC_mutation_line[1]) {
					$amino_acid_var = &codon2aa($hash_aa{$key});		 		
					$amino_acid_wt = &codon2aa($hash_aa_wt{$key});
					if (($amino_acid_wt =~ $amino_acid_var) && ($status=~'ND'))  {
						$status="Susceptible";
					} elsif (($amino_acid_wt !~ $amino_acid_var) && ($amino_acid_var =~ 'Stop')) {
							@HC_aa = split ($split_HC_mutation_line[1],$key);
							$hash_loci{$split_HC_mutation_line[1]} .= "separator$amino_acid_wt$HC_aa[1]$amino_acid_var";
							$status="Resistant";
							$pass=0;
					} elsif (($amino_acid_var =~ 'X') && ($status !~ 'Resistant')) {
						$pass=1;
					}	
					if ($amino_acid_var =~ 'X') {
						@HC_aa = split ($split_HC_mutation_line[1],$key);
						$hash_failed{$split_HC_mutation_line[1]} .= "separator$HC_aa[1]";				
					}
				}		 
			}
		
			#if ($split_HC_mutation_line[1] =~ 'katG|ethA') {  ## flag Isoniazid as passed even in presence of failed position at non-HC position 
			#	$pass=0;
			#}		
		#### look for any non-silent mutation in certain loci flagged 'Any'
		} elsif ($split_HC_mutation_line[8] =~ 'Any_C') {
		 	foreach my $key (keys %hash_aa) {		 
			 	if ($key =~ $split_HC_mutation_line[1]) {
		 			$amino_acid_var = &codon2aa($hash_aa{$key});		 		
		 			$amino_acid_wt = &codon2aa($hash_aa_wt{$key});
		 			@HC_aa = split ($split_HC_mutation_line[1],$key);
					if (($amino_acid_wt =~ $amino_acid_var) && ($status=~'ND'))  {
						$status="Susceptible";
					} elsif (($amino_acid_wt !~ $amino_acid_var) && ($amino_acid_var !~ 'X') && (!exists $neutral_pos{$split_HC_mutation_line[1].$HC_aa[1]."$amino_acid_wt -> $amino_acid_var"})) {
					
					#print "x $split_HC_mutation_line[1].$HC_aa[1]$amino_acid_wt -> $amino_acid_var\n";
					
							#@HC_aa = split ($split_HC_mutation_line[1],$key);
							$hash_loci{$split_HC_mutation_line[1]} .= "separator$amino_acid_wt$HC_aa[1]$amino_acid_var";
							$status="Resistant";
							$pass=0;
					} elsif (($amino_acid_var =~ 'X') && ($status !~ 'Resistant')) {
						$pass=1;
					}	
					if ($amino_acid_var =~ 'X') {
						@HC_aa = split ($split_HC_mutation_line[1],$key);
						$hash_failed{$split_HC_mutation_line[1]} .= "separator$HC_aa[1]";
						
					}
				}		 
			}		
		} elsif ($split_HC_mutation_line[8] =~ 'Any_NC') {
			 foreach my $key (keys %hash_nt) {
				if ($key =~ $split_HC_mutation_line[1]) {
					#my $key_rrs = $split_HC_mutation_line[1];
					
					if (($hash_nt_wt{$key} =~ $hash_nt{$key}) && ($status=~'ND'))  {
						$status="Susceptible";
					
					
					} elsif (($hash_nt_wt{$key} !~ $hash_nt{$key}) && ($hash_nt{$key} !~ 'N') && (!exists $neutral_pos{$key.$HC_aa[1]."$hash_nt_wt{$key} -> $hash_nt{$key}"})) {
						@HC_aa = split ($split_HC_mutation_line[1],$key);

						
						$hash_loci{$split_HC_mutation_line[1]} .= "separator$hash_nt_wt{$key}$HC_aa[1]$hash_nt{$key}";
						$status="Resistant";
						$pass=0;
							
					}	elsif (($hash_nt{$key} =~ 'N') && ($status !~ 'Resistant')) {
						$pass=1;
					}
					if ($hash_nt{$key} =~ 'N') {
						@HC_aa = split ($split_HC_mutation_line[1],$key);
						$hash_failed{$split_HC_mutation_line[1]} .= "separator$HC_aa[1]";
					}
					 
			 	}
			}
	
##### end Any mutation ######
	
		} else {		
		
			if (($split_HC_mutation_line[8] =~ 'Coding')) {
				if ($hash_aa{$split_HC_mutation_line[1].$split_HC_mutation_line[3]}) {
					$amino_acid = &codon2aa($hash_aa{$split_HC_mutation_line[1].$split_HC_mutation_line[3]});
					if ($split_HC_mutation_line[4] =~ $hash_aa{$split_HC_mutation_line[1].$split_HC_mutation_line[3]}) {
						if ($status=~'ND') {
							$status="Susceptible";
						}

					} elsif ($split_HC_mutation_line[5] =~ $hash_aa{$split_HC_mutation_line[1].$split_HC_mutation_line[3]}) {
						if ($split_HC_mutation_line[1] =~ 'rpoB') {
							$split_HC_mutation_line[3] += 81; ## add 81 to position to comply with Ecoli notation
							$hash_loci{$split_HC_mutation_line[1]} .= "separator$split_HC_mutation_line[6]$split_HC_mutation_line[3]$amino_acid";
						} else {
							$hash_loci{$split_HC_mutation_line[1]} .= "separator$split_HC_mutation_line[6]$split_HC_mutation_line[3]$amino_acid";
						}			
						$status="Resistant";
						$pass=0;
					} elsif (($hash_aa{$split_HC_mutation_line[1].$split_HC_mutation_line[3]} =~ 'N') && ($status !~ 'Resistant')) {
						$pass=1;
					}
				} else {
					$pass=1;
					$hash_failed{$split_HC_mutation_line[1]} .= "separator$split_HC_mutation_line[3]";
				}
				if ($hash_aa{$split_HC_mutation_line[1].$split_HC_mutation_line[3]} =~ 'N') {
					$hash_failed{$split_HC_mutation_line[1]} .= "separator$split_HC_mutation_line[3]";
				}
			
			} elsif (($split_HC_mutation_line[8] =~ 'Non-coding')) {  ### look for Non-coding positions ###

				 
				if ($hash_nt{$split_HC_mutation_line[1].$split_HC_mutation_line[2]}) {
		
					if ($split_HC_mutation_line[4] =~ $hash_nt{$split_HC_mutation_line[1].$split_HC_mutation_line[2]}) {
						if ($status=~'ND') {
							$status="Susceptible";
						}
					} elsif ($split_HC_mutation_line[5] =~ $hash_nt{$split_HC_mutation_line[1].$split_HC_mutation_line[2]}) {
							
						$hash_loci{$split_HC_mutation_line[1]} .= "separator$split_HC_mutation_line[4]\($split_HC_mutation_line[3]\)$hash_nt{$split_HC_mutation_line[1].$split_HC_mutation_line[2]}";
						$status="Resistant";
						$pass=0;
					} elsif (($hash_nt{$split_HC_mutation_line[1].$split_HC_mutation_line[2]} =~ 'N') && ($status !~ 'Resistant')){
						$pass=1;
					}
				} else {
					$pass=1;
					$hash_failed{$split_HC_mutation_line[1]} .= "separator$split_HC_mutation_line[2]";
				}		
				if ($hash_nt{$split_HC_mutation_line[1].$split_HC_mutation_line[2]} =~ 'N') {
					$hash_failed{$split_HC_mutation_line[1]} .= "separator$split_HC_mutation_line[2]";
				}			
			} elsif (($split_HC_mutation_line[8] =~ 'INDEL')) {
				if (($hash_indel{$split_HC_mutation_line[1].$split_HC_mutation_line[2]}) && ($split_HC_mutation_line[5] eq $hash_indel{$split_HC_mutation_line[1].$split_HC_mutation_line[2]})) {
					$hash_loci{$split_HC_mutation_line[1]} .= "separator$split_HC_mutation_line[3]In-Frame Indel$split_HC_mutation_line[5]";
					$status="Resistant";
					$pass=0;
				} 			
			}
		}
	}

	if ($pass == 1) {
		$pass = "FAIL";
		$status="ND";

	} elsif ($pass==0) {
		$pass = "PASS";
	} 
	$resis{$screen} = "$screen\t$pass\t$status";
	#print "xx $resis{$screen}\n";
}



### create consensus file section ###
open (CONSENSUS,">$path/Consensus/$prefix.fa");
print CONSENSUS ">$prefix\n";
open (FILE,"$path/tmp.vcf");
my @genome_position=();
my @split_vcf_file;
while (my $in_VCF=<FILE>) {  ## go through VCF, put in hash position{nt} and Status of the position
	if ($in_VCF !~ '\#') {		
		@split_vcf_file = split ('\t',$in_VCF);
		$genome_position[$split_vcf_file[1]] = $in_VCF;
	}
}
close FILE;



my $last_position = $split_vcf_file[1];
my $counter_position=1;
do {
	my @split_vcf_line = split ('\t',$genome_position[$counter_position]);
	if (($split_vcf_line[4] =~ '\.') || (length($split_vcf_line[4]) > 1)) {
		print  CONSENSUS $split_vcf_line[3];
	} elsif (length($split_vcf_line[3]) > 1) {
		++$counter_position;
		my $length = length($split_vcf_line[3]) + $counter_position -2;
		if ($genome_position[$counter_position] =~ '\;AF\=1\.00\;') {
			print CONSENSUS "$split_vcf_line[4]";
		} else {
			print CONSENSUS "N";
		}
		do {
			print CONSENSUS "N";
			++$counter_position;
		}
		until $counter_position > $length ;
		--$counter_position;
	} else {
		if ($genome_position[$counter_position] =~ '\;AF\=1\.00\;') {
			print CONSENSUS "$split_vcf_line[4]";
		} else {
			print CONSENSUS "N";
		}
	}
	++$counter_position;
}
until $counter_position > $last_position;
close CONSENSUS;
print "\n--------------------------------------------------------------------------\nScreening for Previous MTB matches:\n\n";



#$hamming_fork $rundir
#print "$hamming_fork $path/Consensus/$prefix.fa $run_dir\n";
my $snp_cutoff=20;
run [ "$hamming_fork", "$path/Consensus/$prefix.fa", "$run_dir"], ">", \my $stdout;
my @split_results_newline = split ('\n',$stdout);
shift @split_results_newline;
print "Closest Match in Database:\n$split_results_newline[0]\n";
my @split_line_results = split ('\t',$split_results_newline[0]);
my $snps = $split_line_results[1];
my $counter=1;
do {
@split_line_results = split ('\t',$split_results_newline[$counter]);
	if ($split_line_results[1] == $snps) {
		print "$split_results_newline[$counter]\n";
	}

++$counter;
}
until ($split_line_results[1] > $snps);
print "\n";


print "Matches <= $snp_cutoff SNPs in database:\n\n";
my @split_line_results = split ('\t',$split_results_newline[0]);
if ($split_line_results[1] > $snp_cutoff) {
	print "No close matches found\n";
} else {
	print "$split_results_newline[0]\n";
}
shift @split_results_newline;
foreach my $new_line(@split_results_newline) {
	@split_line_results = split ('\t',$new_line);
	if ($split_line_results[1] <= $snp_cutoff) {
		print "$new_line\n";
	}
}


print "\n--------------------------------------------------------------------------\nQC report for all HC loci:\n\n";

my @split_final_annotation_line;
my %loci_pass;
foreach my $next_loci(@loci_list) {
	my $next_loci_renamed = $next_loci;
	#$next_loci_renamed =~ s/rrs1/rrs 512\, 513\, 516\, 906/g;
	#$next_loci_renamed =~ s/rrs2/rrs 1400/g;
	
	
	if ($next_loci =~ 'rrs') {
		if (exists $hash_failed{$next_loci}) {
			if ($hash_failed{$next_loci} =~ 'or512|or513|or516|or906') {
				print "rrs 512, 513, 516, 906: Failed ";
				$loci_pass{'rrs1'} = "FAIL";
				@split_final_annotation_line = split ('separator',$hash_failed{'rrs'});
				shift @split_final_annotation_line;
				@split_final_annotation_line = uniq(@split_final_annotation_line);
				@split_final_annotation_line = sort { $a <=> $b } @split_final_annotation_line;
				$resis{'Streptomycin'} = "Streptomycin\tFAIL\tND";
				print join ";",@split_final_annotation_line;
				print "\n";
			} else {
				print "rrs 512, 513, 516, 906: PASS\n";
				$loci_pass{'rrs1'} = "PASS";
			}
			if ($hash_failed{$next_loci} =~ 'or1400') {
				print "rrs 1400: Failed ";
				$loci_pass{'rrs2'} = "FAIL";
				@split_final_annotation_line = split ('separator',$hash_failed{'rrs'});
				
				shift @split_final_annotation_line;
				@split_final_annotation_line = uniq(@split_final_annotation_line);
				@split_final_annotation_line = sort { $a <=> $b } @split_final_annotation_line;
				$resis{'Kanamycin/Amikacin'} = "Kanamycin/Amikacin\tFAIL\tND";
				print join ";",@split_final_annotation_line;
				print "\n";
			} else {
				print "rrs 1400: PASS\n";
				$loci_pass{'rrs2'} = "PASS";
			}
		
		} else {		
		print "rrs 512, 513, 516, 906: PASS\n";
		print "rrs 1400: PASS\n";
			$loci_pass{'rrs1'} = "PASS";
			$loci_pass{'rrs2'} = "PASS";
			$loci_pass{'rrs'} = "PASS";
		}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	} else {
	
		if (exists $hash_failed{$next_loci}) {
			if ($next_loci =~ 'katG') {  ## special case loop for katG to only call FAIL resistance if a high confidence site 
				if ($hash_failed{$next_loci} =~ 'or279|or315|or525|or394') {
					$loci_pass{$next_loci} = "FAIL";
					print "$next_loci_renamed: Failed "; 
					@split_final_annotation_line = split ('separator',$hash_failed{$next_loci});
					shift @split_final_annotation_line;
					@split_final_annotation_line = uniq(@split_final_annotation_line);
					@split_final_annotation_line = sort { $a <=> $b } @split_final_annotation_line;
					$resis{'Isoniazid'} = "Isoniazid\tFAIL\tND";
				} else {
					$loci_pass{$next_loci} = "PASS";
					print "$next_loci_renamed: Failed "; 
					@split_final_annotation_line = split ('separator',$hash_failed{$next_loci});
					shift @split_final_annotation_line;
					@split_final_annotation_line = uniq(@split_final_annotation_line);
					@split_final_annotation_line = sort { $a <=> $b } @split_final_annotation_line;
				}
			} else {
				$loci_pass{$next_loci} = "FAIL";
				print "$next_loci_renamed: Failed "; 
				@split_final_annotation_line = split ('separator',$hash_failed{$next_loci});
				shift @split_final_annotation_line;
				@split_final_annotation_line = uniq(@split_final_annotation_line);
				@split_final_annotation_line = sort { $a <=> $b } @split_final_annotation_line;
	
			}
			if ($next_loci =~ 'promoter|rrs|rss') {
				print "position ";
			} else {
				print "codon ";
			}
				print join ";",@split_final_annotation_line;
				print "\n";
		} else {
			print "$next_loci_renamed: PASS\n";
			$loci_pass{$next_loci} = "PASS";
		}




	}


}









print "\n-------------------------------------------------------------------------\n";

if ($hash_loci{'rrs'} =~ '512|513|516|906') {
	$hash_loci{'rrs1'} = $hash_loci{'rrs'};
}
if ($hash_loci{'rrs'} =~ '1400') {
	$hash_loci{'rrs2'} = $hash_loci{'rrs'};
	
}

foreach my $key (keys %neutral_pos) {
#	print "$key\t$neutral_pos{$key}\n";
}
#exit;

print "CLIMS SECTION:\n\n";

print "Sample ID\t$prefix\nAnalysis Date\t$date\n\n";
print "Test\tProperty\tNotation\tResult\n";

my @loci_list=('rpoB','katG','oxyR-ahpC promoter region','inhA','mabA-inhA promoter region','mabA','pncA','pncA promoter region','embB','embC-embA promoter region','rpsL','rrs1','rrs2','eis promoter region','gyrA','gyrB','ethA');


foreach my $next_loci(@loci_list) {


	my $next_loci_renamed = $next_loci;
	$next_loci_renamed =~ s/rrs1/rrs 512\, 513\, 516\, 906/g;
	$next_loci_renamed =~ s/rrs2/rrs 1400/g;
	if (exists $hash_loci{$next_loci}) {
		
		@split_final_annotation_line = split ('separator',$hash_loci{$next_loci});			
		shift @split_final_annotation_line;
		@split_final_annotation_line = uniq(@split_final_annotation_line);
		print "WGS_HI_CONF_MUTATION\t$next_loci_renamed\tPASS\t";
		print join ";",@split_final_annotation_line;
		print "\n";
	} else {
		if ($loci_pass{$next_loci} =~ 'FAIL') {
			print "WGS_HI_CONF_MUTATION\t$next_loci_renamed\t$loci_pass{$next_loci}\tFAILED QC. Could not determine presence/absence of mutation in this region\n";
		} else {

			print "WGS_HI_CONF_MUTATION\t$next_loci_renamed\t$loci_pass{$next_loci}\tNo high-confidence mutation\n"; 
		}
	}
}

my $number_of_frameshifts = @frameshift;
if  ($number_of_frameshifts > 0) {
	foreach my $in_frameshift_list (@frameshift) {	
		if ($in_frameshift_list =~ 'pncA|rpoB|katG|ethA') { ## only flag as resistant if frameshift in these loci
			my @split_list = split ('\t',$in_frameshift_list);
			if (!exists $neutral_pos{$split_line[1].$split_line[2].$split_line[3]}) {  ## disregard Neutral mutation
				$resis{$split_list[0]} = "$split_list[0]\t\t\tPASS\tResistant";
			}
		}
	}
}

print "\n";


## same antimicrobial list as before except order is changed
@antimicrobial = ('Rifampin','Isoniazid','Pyrazinamide','Ethambutol','Streptomycin','Kanamycin/Amikacin','Kanamycin','Fluoroquinolones','Ethionamid');
foreach $screen(@antimicrobial) {
	if (($resis{$screen} =~ 'Susceptible') && (exists $unknown_per_antibiotics{$screen}) && ($screen =~ 'Rifampin|Isoniazid|Pyrazinamide|Ethambutol')) {
	#if (($resis{$screen} =~ 'Susceptible') && (exists $unknown_per_antibiotics{$screen})) {
		$resis{$screen} =~ s/Susceptible/Unknown/g;
	}
	$resis{$screen} =~ s/Resistant/RESISTANT (predicted)/g;
	#$resis{$screen} =~ s/Susceptible/Susceptible (sugg/g;
	$resis{$screen} =~ s/\t+/\t/;
	$resis{$screen} =~ s/ND/Not Determined/g;
	print "WGS_TB_RESISTANCE\t$resis{$screen}\n";
#	} 
}
print "\n";

my $arrSize = @lineage;
if ($arrSize==1) {
	print "WGS_LINEAGE\tLINEAGE\t\t$lineage[0]\n";
} else {
	print "WGS_LINEAGE\tLINEAGE\t\tUndetermined\n";
}



exit;



##################################################################
### Subroutines ###

sub average {
	my @array = @_; # save the array passed to this function
	my $sum; # create a variable to hold the sum of the array's values
	foreach (@array) { $sum += $_; } # add each element of the array 
	# to the sum
	return $sum/@array; # divide sum by the number of elements in the
	# array to find the mean
}


### translate Codon to aminoacids
sub codon2aa{
	my($codon)=@_;
	$codon=uc $codon;
	if (($codon =~ 'N') || (!defined($codon))) {
		$hash_codon{$codon} = "XXX";
	} else {
		my(%hash_codon)=('TCA'=>'Ser','TCC'=>'Ser','TCG'=>'Ser','TCT'=>'Ser','TTC'=>'Phe','TTT'=>'Phe','TTA'=>'Leu','TTG'=>'Leu','TAC'=>'Tyr','TAT'=>'Tyr','TAA'=>'Stop','TAG'=>'Stop','TGC'=>'Cys','TGT'=>'Cys','TGA'=>'Stop','TGG'=>'Trp','CTA'=>'Leu','CTC'=>'Leu','CTG'=>'Leu','CTT'=>'Leu','CCA'=>'Pro','CCC'=>'Pro','CCG'=>'Pro','CCT'=>'Pro','CAC'=>'His','CAT'=>'His','CAA'=>'Gln','CAG'=>'Gln','CGA'=>'Arg','CGC'=>'Arg','CGG'=>'Arg','CGT'=>'Arg','ATA'=>'Ile','ATC'=>'Ile','ATT'=>'Ile','ATG'=>'Met','ACA'=>'Thr','ACC'=>'Thr','ACG'=>'Thr','ACT'=>'Thr','AAC'=>'Asn','AAT'=>'Asn','AAA'=>'Lys','AAG'=>'Lys','AGC'=>'Ser','AGT'=>'Ser','AGA'=>'Arg','AGG'=>'Arg','GTA'=>'Val','GTC'=>'Val','GTG'=>'Val','GTT'=>'Val','GCA'=>'Ala','GCC'=>'Ala','GCG'=>'Ala','GCT'=>'Ala','GAC'=>'Asp','GAT'=>'Asp','GAA'=>'Glu','GAG'=>'Glu','GGA'=>'Gly','GGC'=>'Gly','GGG'=>'Gly','GGT'=>'Gly');
		return $hash_codon{$codon};
	}
}

