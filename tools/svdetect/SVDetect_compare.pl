#!/usr/bin/perl -w

=pod

=head1 NAME

SVDetect Compare for Galaxy

Version: 0.8b for Galaxy

=head1 SYNOPSIS

SVDetect_compare.pl links2compare -conf <configuration_file> [-help] [-man]

=cut

# -------------------------------------------------------------------

use strict;
use warnings;

use Pod::Usage;
use Getopt::Long;

use Config::General;
use Tie::IxHash;

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#PARSE THE COMMAND LINE
my %OPT;
GetOptions(\%OPT,
	   'conf=s',
	   'out1=s', #GALAXY
	   'out2=s', #GALAXY
	   'out3=s', #GALAXY
	   'out4=s', #GALAXY
	   'out5=s', #GALAXY
	   'out6=s', #GALAXY
	   'out7=s', #GALAXY
	   'out8=s', #GALAXY
	   'out9=s', #GALAXY
	   'l=s', #GALAXY
	   'N=s', #GALAXY
	   'help',
           'man'
	  );

pod2usage() if $OPT{help};
pod2usage(-verbose=>2) if $OPT{man};
pod2usage(-message=> "$!", -exitval => 2) if (!defined $OPT{conf});


pod2usage() if(@ARGV<1);

tie (my %func, 'Tie::IxHash',links2compare=>\&links2compare);

foreach my $command (@ARGV){
    pod2usage(-message=> "Unknown command \"$command\"", -exitval => 2) if (!defined($func{$command}));
}
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#READ THE CONFIGURATION FILE
my $conf=Config::General->new(    -ConfigFile        => $OPT{conf},
                                  -Tie => "Tie::IxHash",
                                  -AllowMultiOptions => 1,
				  -LowerCaseNames    => 1,
				  -AutoTrue => 1);
my %CONF= $conf->getall;
validateconfiguration(\%CONF);							#validation of the configuration parameters


my $SAMTOOLS_BIN_DIR="/bioinfo/local/samtools"; #GALAXY
my $BEDTOOLS_BIN_DIR="/bioinfo/local/BEDTools/bin"; #GALAXY

my $pt_log_file=$OPT{l}; #GALAXY
my $log_file=$CONF{general}{output_dir}.$OPT{N}.".svdetect_compare.log"; #GALAXY
open LOG,">$log_file" or die "$0: can't open ".$log_file.":$!\n";#GALAXY

my @pt_sv_file=($OPT{out1},$OPT{out2},$OPT{out3}) if($OPT{out1}); #GALAXY common,sample,reference
my @pt_circos_file=($OPT{out4},$OPT{out5},$OPT{out6}) if($OPT{out4}); #GALAXY common,sample,reference
my @pt_bed_file=($OPT{out7},$OPT{out8},$OPT{out9}) if($OPT{out7}); #GALAXY common,sample,reference

$CONF{compare}{sample_link_file}=readlink($CONF{compare}{sample_link_file});#GALAXY
$CONF{compare}{sample_link_file}=~s/.sv.txt//; #GALAXY

$CONF{compare}{reference_link_file}=readlink($CONF{compare}{reference_link_file});#GALAXY
$CONF{compare}{reference_link_file}=~s/.sv.txt//; #GALAXY

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#COMMAND EXECUTION
foreach my $command (@ARGV){
    &{$func{$command}}();
}
print LOG "-- end\n";

close LOG;#GALAXY
system "rm $pt_log_file ; ln -s $log_file $pt_log_file"; #GALAXY

exit(0);
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#FUNCTIONS

# -----------------------------------------------------------------------------#
#MAIN FUNCTION number 5:Comparison between samples, common or specific links
sub links2compare{
    
    my @compare_files;
    
    compareSamples($CONF{general}{output_dir},
		   $CONF{compare}{list_samples},
		   $CONF{compare}{sample_link_file},
		   $CONF{compare}{reference_link_file},
		   $CONF{compare}{min_overlap},
		   $CONF{compare}{same_sv_type},
		   \@compare_files);

    my $pt_ind=0;
 
    for my $input_file (@compare_files){
	
	$input_file=$CONF{general}{output_dir}.$input_file;
	
	my $output_file=$input_file;
	$output_file=~s/unique$/compared/;
	
	sortLinks($input_file, $output_file,1);
	
	if($CONF{compare}{circos_output}){
	    links2segdup($CONF{circos}{organism_id},
			 $CONF{circos}{colorcode},
			 $output_file,
			 $output_file.".segdup.txt");
	    system "rm $pt_circos_file[$pt_ind]; ln -s $output_file.segdup.txt $pt_circos_file[$pt_ind]" if (defined $pt_circos_file[$pt_ind]); #GALAXY
	}
	
	if($CONF{compare}{bed_output}){
	links2bedfile($CONF{compare}{read_lengths},
		      $CONF{bed}{colorcode},
		      $output_file,
		      $output_file.".bed");
	system "rm $pt_bed_file[$pt_ind]; ln -s $output_file.bed $pt_bed_file[$pt_ind]" if (defined $pt_bed_file[$pt_ind]); #GALAXY
	}
	
	if($CONF{compare}{sv_output}){
	    
	    links2SVfile ($output_file, $output_file.".sv.txt");
	    system "rm $pt_sv_file[$pt_ind]; ln -s $output_file.sv.txt $pt_sv_file[$pt_ind]" if (defined $pt_sv_file[$pt_ind]); #GALAXY
	}
	$pt_ind++;
	
    }
    unlink(@compare_files);

}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub compareSamples{
    
    my ($dir,$list_samples,$sample_file,$reference_file,$min_overlap,$same_sv_type,$file_names)=@_;
    
    my @bedpefiles;
    my @list=split(",",$list_samples);
    my @list_files=($sample_file,$reference_file);

    print LOG "\# Comparison procedure...\n";
    print LOG "-- samples=$list_samples\n".
	 "-- minimum overlap=$min_overlap\n".
	 "-- same SV type=$same_sv_type\n";

    #conversion of links to bedPE format file
    print LOG "-- Conversion of links.filtered files to bedPE format\n";
    for my $s (0..$#list) {

	links2bedPElinksfile($list[$s],$list_files[$s],$list_files[$s].".bedpe.txt");
	push(@bedpefiles,$list_files[$s].".bedpe.txt");

    }

    #get common links between all samples compared
    print LOG "-- Getting common links between all samples with BEDTools\n";
    my $common_name=join(".",@list);
    
    my $nb=scalar @list;
    my $command="";
    my $prompt=">";
    
    while ($nb>0){
	
	for my $i (0..$#list_files){
		
	    $command.="$BEDTOOLS_BIN_DIR/pairToPair -type both -f $min_overlap -a ".$list_files[$i].".bedpe.txt";
	    my $pipe=0;
	    
	    for my $j ($i+1..$#list_files){
		
		$command.="| $BEDTOOLS_BIN_DIR/pairToPair -type both -f $min_overlap -a stdin" if($pipe);
		$command.=" -b ".$list_files[$j].".bedpe.txt";
		$pipe=1;
		
	    }

	    $command.=$prompt.$dir.$common_name.".bedpe.tmp;";
	    $prompt=">>";
	    
	    my $first=shift(@list_files); push(@list_files,$first);
	    last;
	}
	$nb--;
    }
    
    system ($command);
    
    push(@bedpefiles,$dir.$common_name.".bedpe.tmp");
    
    #Post comparison to get common links if same type only (as an option)
    open( FILE, "<".$dir.$common_name.".bedpe.tmp") or die "Can't open".$dir.$common_name.".bedpe.tmp : $!";
    open( OUT, ">".$dir.$common_name.".bedpe.unique") or die "Can't write in ".$dir.$common_name.".bedpe.unique : $!";
    
    while(<FILE>){
	my @t=split("\t",$_);
	my $s=(split("_",$t[6]))[0];
	my ($sv1,$sv2)=($t[7],$t[18]);
	splice(@t,11,$#t);
	
	if($same_sv_type){
	    print OUT join("\t",@t)."\n" if($sv1 eq $sv2);
	}else{
	    print OUT join("\t",@t)."\n";
	}
    }
    close FILE;
    close OUT;
    
    bedPElinks2linksfile($dir.$common_name.".bedpe.unique", $dir.$common_name.".unique");
    push(@bedpefiles,$dir.$common_name.".bedpe.unique");
    push(@$file_names,$common_name.".unique");
    print LOG "-- output created: ".$dir.$common_name.".compared\n";
    
    #get specific links for each sample
    print LOG "-- Getting specific links for each sample\n";
    for my $s (0..$#list) {
	system("grep -Fxv -f ".$dir.$common_name.".bedpe.unique ".$list_files[$s].".bedpe.txt >".$dir.$list[$s].".bedpe.unique");
	bedPElinks2linksfile($dir.$list[$s].".bedpe.unique",$dir.$list[$s].".unique");
	push(@bedpefiles,$dir.$list[$s].".bedpe.unique");
	push(@$file_names,$list[$s].".unique");
	print LOG "-- output created: ".$dir.$list[$s].".compared\n";
    }
    
    unlink(@bedpefiles);

}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#convert the links file to the circos format
sub links2segdup{
    
    my($id,$color_code,$links_file,$segdup_file)=@_;
    
    print LOG "\# Converting to the circos format...\n";
    
    tie (my %hcolor,'Tie::IxHash');						#color-code hash table
    foreach my $col (keys %{$color_code}){
	my ($min_links,$max_links)=split(",",$color_code->{$col});
	$hcolor{$col}=[$min_links,$max_links];
    }
    
    open LINKS, "<$links_file" or die "$0: can't open $links_file :$!\n";
    open SEGDUP, ">$segdup_file" or die "$0: can't write in the output: $segdup_file :$!\n";
    
    my $index=1;
    while(<LINKS>){
	
	my ($chr1,$start1,$end1,$chr2,$start2,$end2,$count)=(split)[0,1,2,3,4,5,6];
	
	my $color=getColor($count,\%hcolor,"circos");				#get the color-code according the number of links
	
	print SEGDUP "$index\t$id$chr1\t$start1\t$end1\tcolor=$color\n".	#circos output
		     "$index\t$id$chr2\t$start2\t$end2\tcolor=$color\n";
	$index++;
    }
    
    close LINKS;
    close SEGDUP;
    print LOG "-- output created: $segdup_file\n";
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#convert the links file to the bedPE format for BEDTools usage
sub links2bedPElinksfile{

    my ($sample,$links_file,$bedpe_file)=@_;
    
    open LINKS, "<$links_file" or die "$0: can't open $links_file :$!\n";
    open BEDPE, ">$bedpe_file" or die "$0: can't write in the output: $bedpe_file :$!\n";
    
    my $nb_links=1;
    
    while(<LINKS>){
	
	chomp;
	my @t=split("\t",$_);
	my ($chr1,$start1,$end1,$chr2,$start2,$end2)=splice(@t,0,6);
	my $type=($chr1 eq $chr2)? "INTRA":"INTER";
	$type.="_".$t[10];
	
	$start1--; $start2--;
	
	print BEDPE "$chr1\t$start1\t$end1\t$chr2\t$start2\t$end2".
	"\t$sample"."_link$nb_links\t$type\t.\t.".
	"\t".join("|",@t)."\n";
	
	$nb_links++;
    }
    
    close LINKS;
    close BEDPE;

}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub bedPElinks2linksfile{

    my ($bedpe_file,$links_file)=@_;
 
    open BEDPE, "<$bedpe_file" or die "$0: can't open: $bedpe_file :$!\n";
    open LINKS, ">$links_file" or die "$0: can't write in the output $links_file :$!\n";
    
    while(<BEDPE>){
	
	chomp;
	my $sample=(split("_",(split("\t",$_))[6]))[0];
	my @t1=(split("\t",$_))[0,1,2,3,4,5];
	my @t2=split(/\|/,(split("\t",$_))[10]);
	push(@t2,$sample);
	
	print LINKS join("\t",@t1)."\t".join("\t",@t2)."\n";
	
    }
    close BEDPE;
    close LINKS;

}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#convert the links file to the bed format
sub links2bedfile{
    
    my ($tag_length,$color_code,$links_file,$bed_file)=@_;
    
    print LOG "\# Converting to the bed format...\n";
    
    my $compare=1;
    if($links_file!~/compared$/){
	$compare=0;
	$tag_length->{none}->{1}=$tag_length->{1};
	$tag_length->{none}->{2}=$tag_length->{2};
    }
    
    #color-code hash table
    tie (my %hcolor,'Tie::IxHash');
    my %color_order;
    $color_order{"255,255,255"}=0;
    my $n=1;
    foreach my $col (keys %{$color_code}){
	my ($min_links,$max_links)=split(",",$color_code->{$col});
	$hcolor{$col}=[$min_links,$max_links];
	$color_order{$col}=$n;
	$n++;
    }
    
    my %pair;
    my %pt;
    $n=1;
    open LINKS, "<$links_file" or die "$0: can't open $links_file:$!\n";
    
    my %str=( "F"=>"+", "R"=>"-" );

    while(<LINKS>){
	
	my @t=split;
	my $sample=($compare)? pop(@t):"none";
	
	my $chr1=$t[0]; 
	my $chr2=$t[3];
	$chr1 = "chr".$chr1 unless ($chr1 =~ m/chr/i);
	$chr2 = "chr".$chr2 unless ($chr2 =~ m/chr/i);
	my $same_chr=($chr1 eq $chr2)? 1:0;
	
	my $count=$t[6];
	my $color=getColor($count,\%hcolor,"bed");
	
	my @pairs=deleteBadOrderSensePairs(split(",",$t[7]));
	my @strand1=deleteBadOrderSensePairs(split(",",$t[8]));
	my @strand2=deleteBadOrderSensePairs(split(",",$t[9]));
	my @ends_order1=deleteBadOrderSensePairs(split(",",$t[10]));
	my @ends_order2=deleteBadOrderSensePairs(split(",",$t[11]));
	my @position1=deleteBadOrderSensePairs(split(",",$t[14]));
	my @position2=deleteBadOrderSensePairs(split(",",$t[15]));
	my @start1; my @end1; getCoordswithLeftMost(\@start1,\@end1,\@position1,\@strand1,\@ends_order1,$tag_length->{$sample});
	my @start2; my @end2; getCoordswithLeftMost(\@start2,\@end2,\@position2,\@strand2,\@ends_order2,$tag_length->{$sample});

	
	for my $p (0..$#pairs){						
	    
	    if (!exists $pair{$pairs[$p]}){
		
		if($same_chr){
		    
		    $pair{$pairs[$p]}->{0}=[ $chr1, $start1[$p]-1, $end2[$p], $pairs[$p], 0, $str{$strand1[$p]},
				    $start1[$p]-1, $end2[$p], $color,
				    2, $tag_length->{$sample}->{$ends_order1[$p]}.",".$tag_length->{$sample}->{$ends_order2[$p]}, "0,".($start2[$p]-$start1[$p]) ];    
		    $pt{$n}=$pair{$pairs[$p]}->{0};
		    $n++;
		    
		}else{
		    
		    $pair{$pairs[$p]}->{1}=[ $chr1, $start1[$p]-1, $end1[$p] , $pairs[$p]."/1", 0, $str{$strand1[$p]},
						$start1[$p]-1, $end1[$p], $color,
						1, $tag_length->{$sample}->{$ends_order1[$p]}, 0];
		    $pt{$n}=$pair{$pairs[$p]}->{1};
		    $n++;
		    
		    
		    $pair{$pairs[$p]}->{2}=[ $chr2, $start2[$p]-1, $end2[$p], $pairs[$p]."/2", 0, $str{$strand2[$p]},
						$start2[$p]-1, $end2[$p], $color,
						1, $tag_length->{$sample}->{$ends_order2[$p]}, 0];
		    $pt{$n}=$pair{$pairs[$p]}->{2};
		    $n++;
		}
	    }else{
		
		if($same_chr){
		    ${$pair{$pairs[$p]}->{0}}[8]=$color if($color_order{$color}>$color_order{${$pair{$pairs[$p]}->{0}}[8]});
		}else{
		    ${$pair{$pairs[$p]}->{1}}[8]=$color if($color_order{$color}>$color_order{${$pair{$pairs[$p]}->{1}}[8]});
		    ${$pair{$pairs[$p]}->{2}}[8]=$color if($color_order{$color}>$color_order{${$pair{$pairs[$p]}->{2}}[8]});
		}
	    }
	}
    }
    close LINKS;
    
    my $nb_pairs=$n-1;
    
    open BED, ">$bed_file" or die "$0: can't write in the output: $bed_file :$!\n";
    print BED "track name=\"$bed_file\" description=\"mate pairs involved in links\" ".
	      "visibility=2 itemRgb=\"On\"\n";
    
    for my $i (1..$nb_pairs){
	print BED join("\t",@{$pt{$i}})."\n";
    }
    
    close BED;
    
    print LOG "-- output created: $bed_file\n";
    
    undef %pair;
    undef %pt;
    
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub links2SVfile{
    
    my($links_file,$sv_file)=@_;
    
    print LOG "\# Converting to the sv output table...\n";
    open LINKS, "<$links_file" or die "$0: can't open $links_file:$!\n";
    open SV, ">$sv_file" or die "$0: can't write in the output: $sv_file :$!\n";
    
    my @header=qw(chr_type SV_type BAL_type chromosome1 start1-end1 average_dist
    chromosome2 start2-end2 nb_pairs score_strand_filtering score_order_filtering score_insert_size_filtering
    final_score breakpoint1_start1-end1 breakpoint2_start2-end2);
    
    my $nb_links=0;
    
    while (<LINKS>){
	
	my @t=split;
	my @sv=();
	my $sv_type="-";
	my $strand_ratio="-";
	my $eq_ratio="-";
	my $eq_type="-";
	my $insert_ratio="-";
	my $link="-";
	my ($bk1, $bk2)=("-","-");
	my $score="-";
	
	my ($chr1,$start1,$end1)=($t[0],$t[1],$t[2]); 
	my ($chr2,$start2,$end2)=($t[3],$t[4],$t[5]);
	my $nb_pairs=$t[6];
	$chr1 = "chr".$chr1 unless ($chr1 =~ m/chr/i);
	$chr2 = "chr".$chr2 unless ($chr2 =~ m/chr/i);
	my $chr_type=($chr1 eq $chr2)? "INTRA":"INTER";
	
	#if strand filtering
	if (defined $t[16]){
	    #if inter-chr link
	    $sv_type=$t[16];
	    if(defined $t[17] && $t[17]=~/^(\d+)\/(\d+)$/){
		$strand_ratio=floor($1/$2*100)."%";
		$score=$t[18];
	    }
	    if(defined $t[18] && $t[18]=~/^(\d+)\/(\d+)$/){
	    #if intra-chr link with insert size filtering
		$strand_ratio=floor($1/$2*100)."%";
		$link=floor($t[17]);
		if($sv_type!~/^INV/){
		    $insert_ratio=floor($1/$2*100)."%" if($t[19]=~/^(\d+)\/(\d+)$/);
		    $score=$t[20];
		}else{
		    $score=$t[19];
		}
	    }
	}
	
	if(defined $t[18] && ($t[18] eq "UNBAL" || $t[18] eq "BAL")){
	    
	    #if strand and order filtering only and/or interchr link
	    $eq_type=$t[18];
	    $eq_ratio=floor($1/$2*100)."%" if($t[19]=~/^(\d+)\/(\d+)$/);
	    ($bk1, $bk2)=($t[20],$t[21]);
	    foreach my $bk ($bk1, $bk2){$bk=~s/\),\(/ /g; $bk=~s/(\(|\))//g; $bk=~s/,/-/g;}
	    $score=$t[22];
	    
	}elsif(defined $t[19] && ($t[19] eq "UNBAL" || $t[19] eq "BAL")){
	    
	    #if all three filtering
	    $link=floor($t[17]);
	    $eq_type=$t[19];
	    $eq_ratio=floor($1/$2*100)."%" if($t[20]=~/^(\d+)\/(\d+)$/);
	    
	    if(defined $t[21] && $t[21]=~/^(\d+)\/(\d+)$/){
		$insert_ratio=floor($1/$2*100)."%";
		($bk1, $bk2)=($t[22],$t[23]);
		$score=$t[24];
		
	    }else{
		($bk1, $bk2)=($t[21],$t[22]);
		$score=$t[23];
	    }
	    foreach my $bk ($bk1, $bk2){$bk=~s/\),\(/ /g; $bk=~s/(\(|\))//g; $bk=~s/,/-/g;}
	    
	}
	
	
	push(@sv, $chr_type, $sv_type,$eq_type);
	push(@sv,"$chr1\t$start1-$end1");
	push(@sv, $link);
	push(@sv,"$chr2\t$start2-$end2",
	     $nb_pairs,$strand_ratio,$eq_ratio,$insert_ratio, decimal($score,4), $bk1, $bk2);
	
	
	print SV join("\t",@sv)."\n";
    }
    
    close LINKS;
    close SV;
    
    system "sort  -k 9,9nr -k 13,13nr $sv_file > $sv_file.sorted";
    
    open SV, "<".$sv_file.".sorted" or die "$0: can't open in the output: $sv_file".".sorted :$!\n";
    my @links=<SV>;
    close SV;
    
    open SV, ">$sv_file" or die "$0: can't write in the output: $sv_file :$!\n";
    
    print SV join("\t",@header)."\n";
    print SV @links;
    close SV;
    
    unlink($sv_file.".sorted");
  
    print LOG "-- output created: $sv_file\n";
    
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub deleteBadOrderSensePairs{
    
    my (@tab)=@_;
    my @tab2;

    foreach my $v (@tab){
	
	$v=~s/[\(\)]//g;
	push(@tab2,$v) if($v!~/[\$\*\@]$/);
    }
    return @tab2;
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#gets starts and ends Coords when start=leftmost given positions, directions and orders
sub getCoordswithLeftMost {
    
    my ($starts,$ends,$positions,$strand,$end_order,$tag_length) = @_;

    for my $i (0..scalar(@{$positions})-1) {

	if($strand->[$i] eq 'F'){
	    $starts->[$i]=$positions->[$i];
	    $ends->[$i]=$positions->[$i]+$tag_length->{$end_order->[$i]}-1;
	}else{
	    $starts->[$i]=$positions->[$i]-$tag_length->{$end_order->[$i]}+1;
	    $ends->[$i]=$positions->[$i];
	}
    }    
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub floor {
    my $nb = $_[0];
    $nb=~ s/\..*//;
    return $nb;
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub decimal{
    
  my $num=shift;
  my $digs_to_cut=shift;

  $num=sprintf("%.".($digs_to_cut-1)."f", $num) if ($num=~/\d+\.(\d){$digs_to_cut,}/);

  return $num;
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#Sort links according the concerned chromosomes and their coordinates
sub sortLinks{
    
    my ($links_file,$sortedlinks_file,$unique)=@_;
    
    print LOG "# Sorting links...\n";
    
    my $pipe=($unique)? "| sort -u":"";
    system "sort -k 1,1 -k 4,4 -k 2,2n -k 5,5n -k 8,8n $links_file $pipe > $sortedlinks_file";

}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub getColor{

    my($count,$hcolor,$format)=@_;
    for my $col ( keys % { $hcolor} ) {
       return $col if($count>=$hcolor->{$col}->[0] && $count<=$hcolor->{$col}->[1]);
    }
    return "white" if($format eq "circos");
    return "255,255,255" if($format eq "bed");
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#check if the configuration file is correct
sub validateconfiguration{
    
    my %conf=%{$_[0]};
    my $list_prgs="@ARGV";
    
    my @circos_params=qw(organism_id colorcode);
    my @bed_params=qw(colorcode);
    my @compare_params=qw(list_samples list_read_lengths sample_link_file reference_link_file);
    
    unless (defined($conf{general}{output_dir})) {
	$conf{general}{output_dir} = ".";
    }
    unless (-d $conf{general}{output_dir}){
	mkdir $conf{general}{output_dir} or die;
    }
    $conf{general}{output_dir}.="/" if($conf{general}{output_dir}!~/\/$/);

    
    if($list_prgs=~/links2compare/){
	foreach my $p (@compare_params) {
	    die("Error Config : The compare parameter \"$p\" is not defined\n") if (!defined $conf{compare}{$p});
	}
	
	unless (defined($conf{compare}{same_sv_type})) {
	    $conf{compare}{same_sv_type} = 0;
	}
	
	unless (defined($conf{compare}{min_overlap})) {
	    $conf{compare}{min_overlap} = 1E-9;
	}
	
	if($conf{compare}{circos_output}){
	    foreach my $p (@circos_params) {
		next if($list_prgs=~/^ratio/ && $p eq "colorcode");
		die("Error Config : The circos parameter \"$p\" is not defined\n") if (!defined $conf{circos}{$p});
	    }
	}
	if($conf{compare}{bed_output}){
	    foreach my $p (@bed_params) {
		die("Error Config : The bed parameter \"$p\" is not defined\n") if (!defined $conf{bed}{$p});
	    }
	    die("Error Config : The compare parameter \"list_read_lengths\" is not defined\n") if (!defined $conf{compare}{list_read_lengths});

	    my @samples=split(",",$conf{compare}{list_samples});
	    my @read_lengths=split(",",$conf{compare}{list_read_lengths});
	    for my $i (0..$#samples){
		my @l=split("-",$read_lengths[$i]);
		$conf{compare}{read_lengths}{$samples[$i]}={ 1=> $l[0], 2=> $l[1]};
	    }
	}
    }
   
    
}
#------------------------------------------------------------------------------#
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
