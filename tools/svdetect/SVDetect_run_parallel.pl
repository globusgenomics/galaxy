#!/usr/bin/perl -w

=pod

=head1 NAME

SVDetect - Program designed to the detection of structural variations
from paired-end/mate-pair sequencing data, compatible with SOLiD and Illumina (>=1.3) reads

Version: 0.8b for Galaxy

=head1 SYNOPSIS

SVDetect <command> -conf <configuration_file> [-help] [-man]
    
    Command:

    	linking		detection and isolation of links
        filtering	filtering of links according different parameters
        links2circos	links conversion to circos format
	links2bed 	paired-ends of links converted to bed format (UCSC)
	links2SV	formatted output to show most significant SVs
	cnv		calculate copy-number profiles
	ratio2circos	ratio conversion to circos density format
	ratio2bedgraph	ratio conversion to bedGraph density format (UCSC)
    
=head1 DESCRIPTION

This is a command-line interface to SVDetect.


=head1 AUTHORS

Bruno Zeitouni E<lt>bruno.zeitouni@curie.frE<gt>,
Valentina Boeva E<lt>valentina.boeva@curie.frE<gt>

=cut

# -------------------------------------------------------------------

use strict;
use warnings;

use Pod::Usage;
use Getopt::Long;

use Config::General;
use Tie::IxHash;
use FileHandle;
use Parallel::ForkManager;

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
	   'l=s', #GALAXY
	   'N=s',#GALAXY
	   'help',#GALAXY
           'man'
	  );

pod2usage() if $OPT{help};
pod2usage(-verbose=>2) if $OPT{man};
pod2usage(-message=> "$!", -exitval => 2) if (!defined $OPT{conf});

pod2usage() if(@ARGV<1);

tie (my %func, 'Tie::IxHash',linking=>\&createlinks,
			     filtering=>\&filterlinks,
			     links2circos=>\&links2circos,
			     links2bed=>\&links2bed,
			     links2compare=>\&links2compare,
			     links2SV=>\&links2SV,
			     cnv=>\&cnv,
			     ratio2circos=>\&ratio2circos,
			     ratio2bedgraph=>\&ratio2bedgraph);

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

my $pt_log_file=$OPT{l}; #GALAXY
my $pt_links_file=$OPT{out1} if($OPT{out1}); #GALAXY
my $pt_flinks_file=$OPT{out2} if($OPT{out2}); #GALAXY
my $pt_sv_file=$OPT{out3} if($OPT{out3}); #GALAXY
my $pt_circos_file=$OPT{out4} if($OPT{out4}); #GALAXY
my $pt_bed_file=$OPT{out5} if($OPT{out5}); #GALAXY

$CONF{general}{mates_file}=readlink($CONF{general}{mates_file});#GALAXY
$CONF{general}{cmap_file}=readlink($CONF{general}{cmap_file});#GALAXY

my $log_file=$CONF{general}{output_dir}.$OPT{N}.".svdetect_run.log"; #GALAXY
open LOG,">$log_file" or die "$0: can't open ".$log_file.":$!\n";#GALAXY
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#COMMAND EXECUTION
foreach my $command (@ARGV){
    &{$func{$command}}();
}
print LOG "-- end\n";#GALAXY

close LOG;#GALAXY
system "rm $pt_log_file ; ln -s $log_file $pt_log_file"; #GALAXY
exit(0);
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#FUNCTIONS
#------------------------------------------------------------------------------#
#MAIN FUNCTION number 1: Detection of links from mate-pairs data
sub createlinks{
    
    my %CHR;									#main hash table 1: fragments, links
    my %CHRID;
    my @MATEFILES;
    
    my $output_prefix=$CONF{general}{mates_file}.".".$CONF{general}{sv_type};
    my @path=split(/\//,$output_prefix);
    $output_prefix=$CONF{general}{output_dir}.$path[$#path];
    my $tmp_mates_prefix=$CONF{general}{tmp_dir}."mates/".$path[$#path];
    my $tmp_links_prefix=$CONF{general}{tmp_dir}."links/".$path[$#path];
    
    shearingChromosome(\%CHR, \%CHRID,						#making the genomic fragment library with the detection parameters
		     $CONF{detection}{window_size},
		     $CONF{detection}{step_length},
		     $CONF{general}{cmap_file});
    
    if($CONF{detection}{split_mate_file}){
    
	splitMateFile(\%CHR, \%CHRID, \@MATEFILES, $tmp_mates_prefix,
		      $CONF{general}{sv_type},
		      $CONF{general}{mates_file},
		      $CONF{general}{input_format},
		      $CONF{general}{read_lengths}
		      );
    }else{
	
	@MATEFILES=qx{ls $tmp_mates_prefix*} or die "# Error: No splitted mate files already created at $CONF{general}{tmp_dir} :$!";
	chomp(@MATEFILES);
	print LOG "# Splitted mate files already created.\n";
    }
    
    
    #Parallelization of the linking per chromosome for intra +  interchrs
    my $pm = new Parallel::ForkManager($CONF{general}{num_threads});
    
    foreach my $matefile (@MATEFILES){
	
	my $pid = $pm->start and next;
	getlinks(\%CHR, \%CHRID, $matefile);
	$pm->finish;
    
    }
    $pm->wait_all_children;
    
    #Merge the chromosome links file into only one
    my @LINKFILES= qx{ls $tmp_links_prefix*links} or die "# Error: No links files created at $CONF{general}{tmp_dir} :$!";
    chomp(@LINKFILES);
    catFiles( \@LINKFILES => "$output_prefix.links" );
    
    system "rm $pt_links_file; ln -s $output_prefix.links $pt_links_file" if (defined $pt_links_file); #GALAXY
    print LOG "# Linking end procedure : output created: $output_prefix.links\n";
    #unlink(@LINKFILES);
    #unlink(@MATEFILES);
	
    undef %CHR;
    undef %CHRID;
    
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub getlinks {
    
    my ($chr,$chrID,$tmp_mates_prefix)=@_;
    
    my $tmp_links_prefix=$tmp_mates_prefix;
    $tmp_links_prefix=~s/\/mates\//\/links\//;
    
    my %PAIR;									#main hash table 2: pairs
    
    linking($chr,$chrID, \%PAIR,						#creation of all links from chromosome coordinates of pairs
	  $CONF{general}{read_lengths},
	  $CONF{detection}{window_size},
	  $CONF{detection}{step_length},
	  $tmp_mates_prefix,
	  $CONF{general}{input_format},
	  $CONF{general}{sv_type},
	  "$tmp_links_prefix.links.mapped"
	  );
    
    getUniqueLinks("$tmp_links_prefix.links.mapped",					#remove the doublons
		   "$tmp_links_prefix.links.unique");
    
    defineCoordsLinks($chr,$chrID, \%PAIR,					#definition of the precise coordinates of links
		      $CONF{general}{input_format},
		      $CONF{general}{sv_type},
		      $CONF{general}{read_lengths},
		      "$tmp_links_prefix.links.unique",
		      "$tmp_links_prefix.links.unique_defined"); 
    
    sortLinks("$tmp_links_prefix.links.unique_defined",				#sorting links from coordinates
	      "$tmp_links_prefix.links.sorted");

    removeFullyOverlappedLinks("$tmp_links_prefix.links.sorted",			#remove redundant links
			       "$tmp_links_prefix.links",1);				#file output
    
   
    undef %PAIR;
    
    unlink("$tmp_links_prefix.links.mapped",
	   "$tmp_links_prefix.links.unique",
	   "$tmp_links_prefix.links.unique_defined",
	   "$tmp_links_prefix.links.sorted");
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub splitMateFile{
    
    my ($chr,$chrID,$files_list,$output_prefix,$sv_type,$mates_file,$input_format,$tag_length)=@_;
    
    print LOG "# Splitting the mate file \"$mates_file\" for parallel processing...\n";
    
    my %filesHandle;
    
    #fichier matefile inter
    if($sv_type=~/^(all|inter)$/){
	my $newFileName="$output_prefix.interchrs";
	push(@{$files_list},$newFileName);
	my $fh = new FileHandle;
	$fh->open(">$newFileName");
	$filesHandle{inter}=$fh;
    }
    
    #fichiers matefiles intra
    if($sv_type=~/^(all|intra)$/){
	foreach my $k (1..$chr->{nb_chrs}){
	    my $newFileName=$output_prefix.".".$chr->{$k}->{name};
	    push(@{$files_list},$newFileName);
	    my $fh = new FileHandle;
	    $fh->open(">$newFileName");
	    $filesHandle{$k}=$fh;
	}
    }
    
    if ($mates_file =~ /.gz$/) {
	open(MATES, "gunzip -c $mates_file |") or die "$0: can't open ".$mates_file.":$!\n"; #gzcat
    }elsif($mates_file =~ /.bam$/){
	open(MATES, "$SAMTOOLS_BIN_DIR/samtools view $mates_file |") or die "$0: can't open ".$mates_file.":$!\n";#GALAXY
    }else{
	open MATES, "<".$mates_file or die "$0: can't open ".$mates_file.":$!\n";
    }
    
    
    while(<MATES>){
	
	my @t=split;								
	my ($chr_read1, $chr_read2, $firstbase_read1, $firstbase_read2, $end_order_read1,$end_order_read2);
	
	next if (!readMateFile(\$chr_read1, \$chr_read2, \$firstbase_read1, \$firstbase_read2, \$end_order_read1, \$end_order_read2, \@t, $input_format,$tag_length));
	
	next unless (exists $chrID->{$chr_read1} && exists $chrID->{$chr_read2});
	
	($chr_read1, $chr_read2)= ($chrID->{$chr_read1},$chrID->{$chr_read2});
	
	if( ($sv_type=~/^(all|inter)$/) && ($chr_read1 ne $chr_read2) ){
	    my $fh2print=$filesHandle{inter};
	    print $fh2print join("\t",@t)."\n";
	}
	
	if( ($sv_type=~/^(all|intra)$/) && ($chr_read1 eq $chr_read2) ){
	    my $fh2print=$filesHandle{$chr_read1};
	    print $fh2print join("\t",@t)."\n";
	    
	}
    }
    
    foreach my $name (keys %filesHandle){
	my $fh=$filesHandle{$name};
	$fh->close;
    }
    
    print LOG "# Splitted mate files of \"$mates_file\" created.\n";
}


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub filterlinks{
    
    my %CHR;
    my %CHRID;
    my @LINKFILES;
    my @FLINKFILES;
    
    my $output_prefix=$CONF{general}{mates_file}.".".$CONF{general}{sv_type};
    my @path=split(/\//,$output_prefix);
    $output_prefix=$CONF{general}{output_dir}.$path[$#path];
    my $tmp_links_prefix=$CONF{general}{tmp_dir}."links/".$path[$#path];
    
    createChrHashTables(\%CHR,\%CHRID,
			$CONF{general}{cmap_file});
    
    if($CONF{filtering}{split_link_file}){
    
	splitLinkFile(\%CHR, \%CHRID, \@LINKFILES,
		      $tmp_links_prefix,
		      $CONF{general}{sv_type},
		      "$output_prefix.links",
		      );
    }else{
	
	@LINKFILES=qx{ls $tmp_links_prefix*links} or die "# Error: No splitted link files already created\n";
	chomp(@LINKFILES);
	print LOG "# Splitted link files already created.\n";
    }
    
    #Parallelization of the filtering per chromosome for intra +  interchrs
    my $pm = new Parallel::ForkManager($CONF{general}{num_threads});
    
    foreach my $linkfile (@LINKFILES){
	
	my $pid = $pm->start and next;
	getFilteredlinks(\%CHR, \%CHRID, $linkfile);
	$pm->finish;
    
    }
    $pm->wait_all_children;
    
    #Merge the chromosome links file into only one
    @FLINKFILES= qx{ls $tmp_links_prefix*filtered} or die "# Error: No links files created\n";
    chomp(@FLINKFILES);
    catFiles( \@FLINKFILES => "$output_prefix.links.filtered" );
    
    system "rm $pt_flinks_file; ln -s $output_prefix.links.filtered $pt_flinks_file" if (defined $pt_flinks_file); #GALAXY
    print LOG"# Filtering end procedure : output created: $output_prefix.links.filtered\n";
    
    undef %CHR;
    undef %CHRID;
    
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub splitLinkFile{
    
    my ($chr,$chrID,$files_list,$input_prefix,$sv_type,$link_file)=@_;
    
    print LOG "# Splitting the link file for parallel processing...\n";
    
    my %filesHandle;
    
    #fichier matefile inter
    if($sv_type=~/^(all|inter)$/){
	my $newFileName="$input_prefix.interchrs.links";
	push(@{$files_list},$newFileName);
	my $fh = new FileHandle;
	$fh->open(">$newFileName");
	$filesHandle{inter}=$fh;
    }
    
    #fichiers matefiles intra
    if($sv_type=~/^(all|intra)$/){
	foreach my $k (1..$chr->{nb_chrs}){
	    my $newFileName=$input_prefix.".".$chr->{$k}->{name}.".links";
	    push(@{$files_list},$newFileName);
	    my $fh = new FileHandle;
	    $fh->open(">$newFileName");
	    $filesHandle{$k}=$fh;
	}
    }
    
    open LINKS, "<".$link_file or die "$0: can't open ".$link_file.":$!\n";
    while(<LINKS>){
	
	my @t=split;
	my ($chr_read1,$chr_read2)=($t[0],$t[3]);							
	
	next unless (exists $chrID->{$chr_read1} && exists $chrID->{$chr_read2});
	
	($chr_read1, $chr_read2)= ($chrID->{$chr_read1},$chrID->{$chr_read2});
	
	if( ($sv_type=~/^(all|inter)$/) && ($chr_read1 ne $chr_read2) ){
	    my $fh2print=$filesHandle{inter};
	    print $fh2print join("\t",@t)."\n";
	}
	
	if( ($sv_type=~/^(all|intra)$/) && ($chr_read1 eq $chr_read2) ){
	    my $fh2print=$filesHandle{$chr_read1};
	    print $fh2print join("\t",@t)."\n";
	    
	}
    }
    
    foreach my $name (keys %filesHandle){
	my $fh=$filesHandle{$name};
	$fh->close;
    }
    
    print LOG "# Splitted link files created.\n";
}


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#MAIN FUNCTION number 2: Filtering processing
sub getFilteredlinks {
    
    my ($chr,$chrID,$tmp_links_prefix)=@_;
    my %PAIR;

    strandFiltering($chr,$chrID,
		     $CONF{filtering}{nb_pairs_threshold},			#filtering of links
		     $CONF{filtering}{strand_filtering},
		     $CONF{filtering}{chromosomes},
		     $CONF{general}{input_format},
		     $CONF{general}{cmap_file},
		     $CONF{general}{mates_orientation},
		     $CONF{general}{read_lengths},
		    $tmp_links_prefix,
		    "$tmp_links_prefix.filtered",
		    ); 						

    if($CONF{filtering}{strand_filtering}){					#re-definition of links coordinates with strand filtering

	my @tmpfiles;
	
	rename("$tmp_links_prefix.filtered","$tmp_links_prefix.filtered_unique");
	
	getUniqueLinks("$tmp_links_prefix.filtered_unique",
		       "$tmp_links_prefix.filtered");
	
	push(@tmpfiles,"$tmp_links_prefix.filtered_unique");
	
	if($CONF{filtering}{order_filtering}){					#filtering using the order
	    
	    rename("$tmp_links_prefix.filtered","$tmp_links_prefix.filtered_ordered");
	    
	    orderFiltering($chr,$chrID,
			      $CONF{filtering}{nb_pairs_threshold},
			      $CONF{filtering}{nb_pairs_order_threshold},			   
			      $CONF{filtering}{mu_length},
			      $CONF{filtering}{sigma_length},
			      $CONF{general}{mates_orientation},	      
			      $CONF{general}{read_lengths},
			      "$tmp_links_prefix.filtered_ordered", 
			      "$tmp_links_prefix.filtered",				
			     );
	    
	    push(@tmpfiles,"$tmp_links_prefix.filtered_ordered");
	}
	
	if (($CONF{filtering}{insert_size_filtering})&&
	    ($CONF{general}{sv_type} ne 'inter')){ 				 
	    
	    rename("$tmp_links_prefix.filtered","$tmp_links_prefix.filtered_withoutIndelSize");
	    
	    addInsertionInfo($chr,$chrID,
			     $CONF{filtering}{nb_pairs_threshold},
			     $CONF{filtering}{order_filtering},
			     $CONF{filtering}{indel_sigma_threshold},
			     $CONF{filtering}{dup_sigma_threshold},
			     $CONF{filtering}{singleton_sigma_threshold}, 
			     $CONF{filtering}{mu_length},
			     $CONF{filtering}{sigma_length},
			     $CONF{general}{mates_orientation},
			     $CONF{general}{read_lengths},
			     "$tmp_links_prefix.filtered_withoutIndelSize", 
			     "$tmp_links_prefix.filtered"			
			     );
	    
	    push(@tmpfiles,"$tmp_links_prefix.filtered_withoutIndelSize");    
	}
	
	sortLinks("$tmp_links_prefix.filtered",
		  "$tmp_links_prefix.filtered_sorted");
	
	removeFullyOverlappedLinks("$tmp_links_prefix.filtered_sorted",
			           "$tmp_links_prefix.filtered_nodup",
				   );				
	
	postFiltering("$tmp_links_prefix.filtered_nodup",
		      "$tmp_links_prefix.filtered",
		       $CONF{filtering}{final_score_threshold});
	
	push(@tmpfiles,"$tmp_links_prefix.filtered_sorted","$tmp_links_prefix.filtered_nodup");
	
	unlink(@tmpfiles);
	
	
    }
    undef %PAIR;
    
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#MAIN FUNCTION number 3: Circos format conversion for links
sub links2circos{
    
    my $input_file=$CONF{general}{mates_file}.".".$CONF{general}{sv_type}.".links.filtered";
    my @path=split(/\//,$input_file);
    $input_file=$CONF{general}{output_dir}.$path[$#path];
    
    my $output_file.=$input_file.".segdup.txt";
    
    links2segdup($CONF{circos}{organism_id},
		 $CONF{circos}{colorcode},
		 $input_file,
		 $output_file);							#circos file output
    
    system "rm $pt_circos_file; ln -s $output_file $pt_circos_file" if (defined $pt_circos_file); #GALAXY
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#MAIN FUNCTION number 4: Bed format conversion for links
sub links2bed{
    
    my $input_file=$CONF{general}{mates_file}.".".$CONF{general}{sv_type}.".links.filtered";
    my @path=split(/\//,$input_file);
    $input_file=$CONF{general}{output_dir}.$path[$#path];
    
    my $output_file.=$input_file.".bed.txt";
    
    links2bedfile($CONF{general}{read_lengths},
		  $CONF{bed}{colorcode},
		  $input_file,
		  $output_file);						#bed file output
    
    system "rm $pt_bed_file; ln -s $output_file $pt_bed_file" if (defined $pt_bed_file); #GALAXY
    
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#MAIN FUNCTION number 6: Bed format conversion for links
sub links2SV{
    
    my $input_file=$CONF{general}{mates_file}.".".$CONF{general}{sv_type}.".links.filtered";
    
    my @path=split(/\//,$input_file);
    $input_file=$CONF{general}{output_dir}.$path[$#path];
    
    my $output_file.=$input_file.".sv.txt";
    
    
    links2SVfile( $input_file,
		  $output_file);
    
     system "rm $pt_sv_file; ln -s $output_file $pt_sv_file" if (defined $pt_sv_file); #GALAXY	
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#MAIN FUNCTION number 7: copy number variations, coverage ratio calculation
sub cnv{
    
    my %CHR;
    my %CHRID;
    my @MATEFILES;
    my @MATEFILES_REF;
    
    my $output_prefix=$CONF{general}{mates_file};
    my $output_prefix_ref=$CONF{detection}{mates_file_ref};
    my @path=split(/\//,$output_prefix);
    my @path_ref=split(/\//,$output_prefix_ref);
    $output_prefix=$CONF{general}{output_dir}.$path[$#path];
    $output_prefix_ref=$CONF{general}{output_dir}.$path_ref[$#path_ref];
    my $tmp_mates_prefix=$CONF{general}{tmp_dir}."mates/".$path[$#path];
    my $tmp_mates_prefix_ref=$CONF{general}{tmp_dir}."mates/".$path_ref[$#path_ref];
    my $tmp_density_prefix=$CONF{general}{tmp_dir}."density/".$path[$#path];
    
    shearingChromosome(\%CHR, \%CHRID,						#making the genomic fragment library with the detection parameters
		     $CONF{detection}{window_size},
		     $CONF{detection}{step_length},
		     $CONF{general}{cmap_file});
    
    if($CONF{detection}{split_mate_file}){
    
	splitMateFile(\%CHR, \%CHRID, \@MATEFILES, $tmp_mates_prefix,
		      "intra",
		      $CONF{general}{mates_file},
		      $CONF{general}{input_format},
		      $CONF{general}{read_lengths}
		      );
	
	splitMateFile(\%CHR, \%CHRID, \@MATEFILES_REF, $tmp_mates_prefix_ref,
		      "intra",
		      $CONF{detection}{mates_file_ref},
		      $CONF{general}{input_format},
		      $CONF{general}{read_lengths}
		      );
	
	
    }else{
	
	@MATEFILES=qx{ls $tmp_mates_prefix*} or die "# Error: No splitted sample mate files of \"$CONF{general}{mates_file}\" already created at $CONF{general}{tmp_dir} :$!";
	chomp(@MATEFILES);
	@MATEFILES_REF=qx{ls $tmp_mates_prefix_ref*} or die "# Error: No splitted reference mate files of \"$CONF{detection}{mates_file_ref}\" already created at $CONF{general}{tmp_dir} :$!";
	chomp(@MATEFILES_REF);
	print LOG "# Splitted sample and reference mate files already created.\n";
    }
    
    #Parallelization of the cnv per chromosome
    my $pm = new Parallel::ForkManager($CONF{general}{num_threads});
    
    foreach my $file (0..$#MATEFILES){
	
	my $pid = $pm->start and next;

	densityCalculation(\%CHR, \%CHRID, $file,
	  $CONF{general}{read_lengths},
	  $CONF{detection}{window_size},
	  $CONF{detection}{step_length},
	  \@MATEFILES,
	  \@MATEFILES_REF,
	  $MATEFILES[$file].".density",
	  $CONF{general}{input_format});

	$pm->finish;
	
    }
    $pm->wait_all_children;
    
    #Merge the chromosome links file into only one
    my @DENSITYFILES= qx{ls $tmp_density_prefix*density} or die "# Error: No density files created at $CONF{general}{tmp_dir} :$!";
    chomp(@DENSITYFILES);
    catFiles( \@DENSITYFILES => "$output_prefix.density" );

    print LOG "# cnv end procedure : output created: $output_prefix.density\n";
    
    
    undef %CHR;
    undef %CHRID;
    
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#MAIN FUNCTION number 8: Circos format conversion for cnv ratios
sub ratio2circos{
    
    my $input_file=$CONF{general}{mates_file}.".density";
    
    my @path=split(/\//,$input_file);
    $input_file=$CONF{general}{output_dir}.$path[$#path];
    
    my $output_file.=$input_file.".segdup.txt";
    
    ratio2segdup($CONF{circos}{organism_id},		
		 $input_file,
		 $output_file);
}
#------------------------------------------------------------------------------#
#MAIN FUNCTION number 9: BedGraph format conversion for cnv ratios
sub ratio2bedgraph{
    
    my $input_file=$CONF{general}{mates_file}.".density";
    
    my @path=split(/\//,$input_file);
    $input_file=$CONF{general}{output_dir}.$path[$#path];
    
    my $output_file.=$input_file.".bedgraph.txt";
    
    ratio2bedfile($input_file,
		  $output_file);						#bed file output
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#Creation of the fragment library
sub shearingChromosome{
    
    print LOG "# Making the fragments library...\n";
    
    my ($chr,$chrID,$window,$step,$cmap_file)=@_;				#window and step sizes parameters

    createChrHashTables($chr,$chrID,$cmap_file);				#hash tables: chromosome ID <=> chromsomes Name
    
    foreach my $k (1..$chr->{nb_chrs}){						
	
	print LOG"-- $chr->{$k}->{name}\n";
	
	my $frag=1;
	for (my $start=0; $start<$chr->{$k}->{length}; $start+=$step){
	    
	    my $end=($start<($chr->{$k}->{length})-$window)? $start+$window-1:($chr->{$k}->{length})-1;
	    $chr->{$k}->{$frag}=[$start,$end];					#creation of fragments, coordinates storage
	    
	    if($end==($chr->{$k}->{length})-1){
		$chr->{$k}->{nb_frag}=$frag;					#nb of fragments per chromosome
		last;
	    }
	    $frag++;
	}
    }
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#Creation of chromosome hash tables from the cmap file
sub createChrHashTables{
    
    my ($chr,$chrID,$cmap_file)=@_;
    $chr->{nb_chrs}=0;
    
    open CMAP, "<".$cmap_file or die "$0: can't open ".$cmap_file.":$!\n";
    while(<CMAP>){

	if(/^\s+$/){ next;}
	my ($k,$name,$length) = split;
	$chr->{$k}->{name}=$name;
	$chr->{$k}->{length}=$length;
	$chrID->{$name}=$k;
	$chr->{nb_chrs}++;
	
    }
    close CMAP;
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#Read the mate file according the input format file (solid, eland or sam)
sub readMateFile{
    
    my ($chr1,$chr2,$pos1,$pos2,$order1,$order2,$t,$file_type,$tag_length)=@_;
    my ($strand1,$strand2);
    
    if($file_type eq "solid"){
	
	($$chr1,$$chr2,$$pos1,$$pos2,$$order1,$$order2)=($$t[6],$$t[7],$$t[8]+1,$$t[9]+1,1,2); #0-based
	
    }else{
	my ($tag_length1,$tag_length2);
	($$chr1,$$chr2,$$pos1,$strand1,$$pos2,$strand2,$$order1,$$order2,$tag_length1,$tag_length2)=($$t[11],$$t[12],$$t[7],$$t[8],$$t[9],$$t[10],1,2,length($$t[1]),length($$t[2])) #1-based
	if($file_type eq "eland");
	
	if($file_type eq "sam"){
	    
	    return 0 if ($$t[0]=~/^@/);						#header sam filtered out
	   
	    ($$chr1,$$chr2,$$pos1,$$pos2)=($$t[2],$$t[6],$$t[3],$$t[7]);
	   
	    return 0 if (($$t[1]&0x0004) || ($$t[1]&0x0008));
	   
	    $$chr2=$$chr1 if($$chr2 eq "=");
	   
	    $strand1 = (($$t[1]&0x0010))? 'R':'F';
	    $strand2 = (($$t[1]&0x0020))? 'R':'F';
	    
	    $$order1=  (($$t[1]&0x0040))? '1':'2';
	    $$order2=  (($$t[1]&0x0080))? '1':'2';
	    $tag_length1 = $tag_length->{$$order1};
	    $tag_length2 = $tag_length->{$$order2};
	}
	
	$$pos1 = -($$pos1+$tag_length1) if ($strand1 eq "R");		#get sequencing starts
	$$pos2 = -($$pos2+$tag_length2) if ($strand2 eq "R");
    }
    return 1;
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#Parsing of the mates files and creation of links between 2 chromosomal fragments
sub linking{
    
    my ($chr,$chrID,$pair,$tag_length,$window_dist,$step,$mates_file,$input_format,$sv_type,$links_file)=@_;
    my %link;
    
    my $record=0;
    my $nb_links=0;
    my $warn=10000;
    
    my @sfile=split(/\./,$mates_file);
    my $fchr=$sfile[$#sfile];
    
    my $fh = new FileHandle;
   
    print LOG "# $fchr : Linking procedure...\n";
    print LOG "-- file=$mates_file\n".
	 "-- chromosome=$fchr\n".
	 "-- input format=$input_format\n".
	 "-- type=$sv_type\n".
	 "-- read1 length=$tag_length->{1}, read2 length=$tag_length->{2}\n".
	 "-- window size=$window_dist, step length=$step\n";
    
    if ($mates_file =~ /.gz$/) {
	$fh->open("gunzip -c $mates_file |") or die "$0: can't open ".$mates_file.":$!\n"; #gzcat
    }elsif($mates_file =~ /.bam$/){
	$fh->open("$SAMTOOLS_BIN_DIR/samtools view $mates_file |") or die "$0: can't open ".$mates_file.":$!\n";#GALAXY
    }else{
	$fh->open("<".$mates_file) or die "$0: can't open ".$mates_file.":$!\n";
    }
    
    
    while(<$fh>){
	
	my @t=split;								#for each mate-pair
	my $mate=$t[0];
	my ($chr_read1, $chr_read2, $firstbase_read1, $firstbase_read2, $end_order_read1,$end_order_read2);
	
	next if(exists $$pair{$mate});
	
	next if (!readMateFile(\$chr_read1, \$chr_read2, \$firstbase_read1, \$firstbase_read2, \$end_order_read1, \$end_order_read2, \@t, $input_format,$tag_length));
	
	next unless (exists $chrID->{$chr_read1} && exists $chrID->{$chr_read2});
	
	($chr_read1, $chr_read2)= ($chrID->{$chr_read1},$chrID->{$chr_read2});
	
	if($sv_type ne "all"){
	    if( ($sv_type eq "inter") && ($chr_read1 ne $chr_read2) ||
	        ($sv_type eq "intra") && ($chr_read1 eq $chr_read2) ){
	    }else{
		next;
	    }
	}
	
	$$pair{$mate}=[$chr_read1, $chr_read2, $firstbase_read1, $firstbase_read2, $end_order_read1, $end_order_read2 ]; 	#fill out the hash pair table (ready for the defineCoordsLinks function)
	
	$record++; 
	
	my ($coord_start_read1,$coord_end_read1,$coord_start_read2,$coord_end_read2);	#get the coordinates of each read
	
	recupCoords($firstbase_read1,\$coord_start_read1,\$coord_end_read1,$tag_length->{$end_order_read1},$input_format);
	recupCoords($firstbase_read2,\$coord_start_read2,\$coord_end_read2,$tag_length->{$end_order_read2},$input_format);
	
	for(my $i=1;$i<=$chr->{$chr_read1}->{'nb_frag'};$i++){			#fast genome parsing for link creation
	    
	    if (abs ($coord_start_read1-${$chr->{$chr_read1}->{$i}}[0]) <= $window_dist){
		
		if(overlap($coord_start_read1,$coord_end_read1,${$chr->{$chr_read1}->{$i}}[0],${$chr->{$chr_read1}->{$i}}[1])){
		    
		    for(my $j=1;$j<=$chr->{$chr_read2}->{'nb_frag'};$j++){
			
			if (abs ($coord_start_read2-${$chr->{$chr_read2}->{$j}}[0]) <= $window_dist) {
			    
			    if(overlap($coord_start_read2,$coord_end_read2,${$chr->{$chr_read2}->{$j}}[0],${$chr->{$chr_read2}->{$j}}[1])){
				
				makeLink(\%link,$chr_read1,$i,$chr_read2,$j,$mate,\$nb_links);	#make the link
			    }
			    
			}else{
			    
			    $j=getNextFrag($coord_start_read2,$j,${$chr->{$chr_read2}->{$j}}[0],$chr->{$chr_read2}->{nb_frag},$window_dist,$step);
			}				
		    }
		}
		
	    }else{
		
		$i=getNextFrag($coord_start_read1,$i,${$chr->{$chr_read1}->{$i}}[0],$chr->{$chr_read1}->{nb_frag},$window_dist,$step);
	    }
	}
	
	if($record>=$warn){
	    print LOG "-- $fchr : $warn mate-pairs analysed - $nb_links links done\n";
	    $warn+=10000;
	}
    }
    $fh->close;
    
    if(!$nb_links){
	print LOG "-- $fchr : No mate-pairs !\n".
	 "-- $fchr : No links have been found with the selected type of structural variations \($sv_type\)\n";
    }
    
    print LOG "-- $fchr : Total : $record mate-pairs analysed - $nb_links links done\n";
    
    print LOG "-- $fchr : writing...\n";
    
    $fh = new FileHandle;
    
    $fh->open(">".$links_file) or die "$0: can't write in the output ".$links_file." :$!\n";
    
    foreach my $chr1 ( sort { $a <=> $b} keys %link){				#Sorted links output
	
	foreach my $chr2 ( sort { $a <=> $b} keys %{$link{$chr1}}){
	    
	    foreach my $frag1 ( sort { $a <=> $b} keys %{$link{$chr1}{$chr2}}){
		
		foreach my $frag2 ( sort { $a <=> $b} keys %{$link{$chr1}{$chr2}{$frag1}}){
		
		    my @count=split(",",$link{$chr1}{$chr2}{$frag1}{$frag2});
		    print $fh "$chr->{$chr1}->{name}\t".(${$chr->{$chr1}->{$frag1}}[0]+1)."\t".(${$chr->{$chr1}->{$frag1}}[1]+1)."\t".
				"$chr->{$chr2}->{name}\t".(${$chr->{$chr2}->{$frag2}}[0]+1)."\t".(${$chr->{$chr2}->{$frag2}}[1]+1)."\t".
				scalar @count."\t".				#nb of read
				$link{$chr1}{$chr2}{$frag1}{$frag2}."\n";	#mate list
		}  
	    }
	}
    }
    
    $fh->close;
    
    undef %link;

}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#remove exact links doublons according to the mate list
sub getUniqueLinks{
    
    my ($links_file,$nrlinks_file)=@_;
    my %links;
    my %pt;
    my $nb_links;
    my $n=1;
    
    my $record=0;
    my $warn=300000;
    
    my @sfile=split(/\./,$links_file);
    my $fchr=$sfile[$#sfile-2];
    
    my $fh = new FileHandle;

    print LOG "# $fchr : Getting unique links...\n";
    $fh->open("<$links_file") or die "$0: can't open $links_file :$!\n";
   
    while(<$fh>){
	
	my @t=split;
	my $mates=$t[7];
	$record++;
	
	if(!exists $links{$mates}){						#Unique links selection
	       
	    $links{$mates}=[@t];
	    $pt{$n}=$links{$mates};
	    $n++;
	    
	
	}else{									#get the link coordinates from the mate-pairs list
	    
	    for my $i (1,2,4,5){						#get the shortest regions
		
		$links{$mates}->[$i]=($t[$i]>$links{$mates}->[$i])? $t[$i]:$links{$mates}->[$i]	#maximum start
		if($i==1 || $i==4);
		$links{$mates}->[$i]=($t[$i]<$links{$mates}->[$i])? $t[$i]:$links{$mates}->[$i]	#minimum end
		if($i==2 || $i==5);
	    }
	}
	if($record>=$warn){
	    print LOG "-- $fchr : $warn links analysed - ".($n-1)." unique links done\n";
	    $warn+=300000;
	}
    }
    $fh->close;
    
    $nb_links=$n-1;
    print LOG "-- $fchr : Total : $record links analysed - $nb_links unique links done\n";
    
    $fh = new FileHandle;
    $fh->open(">$nrlinks_file") or die "$0: can't write in the output: $nrlinks_file :$!\n";
    print LOG "-- $fchr : writing...\n";
    for my $i (1..$nb_links){
	
	print $fh join("\t",@{$pt{$i}})."\n";				#all links output
    }
    
    $fh->close;
    
    undef %links;
    undef %pt;
    
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#get the new coordinates of each link from the mate list
sub defineCoordsLinks{
    
    my ($chr,$chrID,$pair,$input_format,$sv_type,$tag_length,$links_file,$clinks_file)=@_;
    
    my @sfile=split(/\./,$links_file);
    my $fchr=$sfile[$#sfile-2];
    
    my $fh = new FileHandle;
    my $fh2 = new FileHandle;
    
    $fh->open("<$links_file") or die "$0: can't open $links_file :$!\n";
    $fh2->open(">$clinks_file") or die "$0: can't write in the output: $clinks_file :$!\n";
    
    print LOG "# $fchr : Defining precise link coordinates...\n";
    
    my $record=0;
    my $warn=100000;
    
    my %coords;
    my %strands;
    my %order;
    my %ends_order;
    
    while(<$fh>){
	

	my ($col1,$col2)=(1,2);							#for an intrachromosomal link
	my $diffchr=0;								#difference between chr1 and chr2
	my ($chr1,$chr2,$mates_list,$npairs)=(split)[0,3,7,8];
	($chr1,$chr2) = ($chrID->{$chr1},$chrID->{$chr2});
	if ($chr1 != $chr2){							#for an interchromosomal link
	    $col1=$col2=0;							#no distinction
	    $diffchr=1;
	}
	
	my @pairs=split(",",$mates_list);
	
	$coords{$col1}{$chr1}->{start}=undef;
	$coords{$col1}{$chr1}->{end}=undef;
	$coords{$col2}{$chr2}->{start}=undef;
	$coords{$col2}{$chr2}->{end}=undef;
	$strands{$col1}{$chr1}=undef;
	$strands{$col2}{$chr2}=undef;
	$ends_order{$col1}{$chr1}=undef;
	$ends_order{$col2}{$chr2}=undef;
	
	
	$order{$col1}{$chr1}->{index}->{1}=undef;
	$order{$col1}{$chr1}->{index}->{2}=undef;
	$order{$col2}{$chr2}->{index}->{1}=undef;
	$order{$col2}{$chr2}->{index}->{2}=undef;
	$order{$col1}{$chr1}->{order}=undef;
	$order{$col2}{$chr2}->{order}=undef;
	
	$record++;

	for my $p (0..$#pairs){							#for each pair
	    
	    my ($coord_start_read1,$coord_end_read1);
	    my ($coord_start_read2,$coord_end_read2);
	    my $strand_read1=recupCoords(${$$pair{$pairs[$p]}}[2],\$coord_start_read1,\$coord_end_read1,$tag_length->{${$$pair{$pairs[$p]}}[4]},$input_format);
	    my $strand_read2=recupCoords(${$$pair{$pairs[$p]}}[3],\$coord_start_read2,\$coord_end_read2,$tag_length->{${$$pair{$pairs[$p]}}[5]},$input_format);

	    if(!$diffchr){							#for a intrachromosomal link
		if($coord_start_read2<$coord_start_read1){				#get the closer start coordinate for each column
		    ($col1,$col2)=(2,1);
		}else{
		    ($col1,$col2)=(1,2);
		}
	    }
	    
	    push(@{$coords{$col1}{${$$pair{$pairs[$p]}}[0]}->{start}},$coord_start_read1);	#get coords and strands of f3 and r3 reads
	    push(@{$coords{$col1}{${$$pair{$pairs[$p]}}[0]}->{end}},$coord_end_read1);
	    push(@{$coords{$col2}{${$$pair{$pairs[$p]}}[1]}->{start}},$coord_start_read2);
	    push(@{$coords{$col2}{${$$pair{$pairs[$p]}}[1]}->{end}},$coord_end_read2);
	    push(@{$strands{$col1}{${$$pair{$pairs[$p]}}[0]}},$strand_read1);
	    push(@{$strands{$col2}{${$$pair{$pairs[$p]}}[1]}},$strand_read2);
	    push(@{$ends_order{$col1}{${$$pair{$pairs[$p]}}[0]}},${$$pair{$pairs[$p]}}[4]);
	    push(@{$ends_order{$col2}{${$$pair{$pairs[$p]}}[1]}},${$$pair{$pairs[$p]}}[5]);
	}
	
	($col1,$col2)=(1,2) if(!$diffchr);
	
	my $coord_start_chr1=min(min(@{$coords{$col1}{$chr1}->{start}}),min(@{$coords{$col1}{$chr1}->{end}}));		#get the biggest region
	my $coord_end_chr1=max(max(@{$coords{$col1}{$chr1}->{start}}),max(@{$coords{$col1}{$chr1}->{end}}));		
	my $coord_start_chr2=min(min(@{$coords{$col2}{$chr2}->{start}}),min(@{$coords{$col2}{$chr2}->{end}}));	
	my $coord_end_chr2=max(max(@{$coords{$col2}{$chr2}->{start}}),max(@{$coords{$col2}{$chr2}->{end}}));
	
	@{$order{$col1}{$chr1}->{index}->{1}}= sort {${$coords{$col1}{$chr1}->{start}}[$a] <=>  ${$coords{$col1}{$chr1}->{start}}[$b]} 0 .. $#{$coords{$col1}{$chr1}->{start}};
	@{$order{$col2}{$chr2}->{index}->{1}}= sort {${$coords{$col2}{$chr2}->{start}}[$a] <=>  ${$coords{$col2}{$chr2}->{start}}[$b]} 0 .. $#{$coords{$col2}{$chr2}->{start}};
	
	foreach my $i (@{$order{$col1}{$chr1}->{index}->{1}}){			#get the rank of the chr2 reads according to the sorted chr1 reads (start coordinate sorting)
	    foreach my $j (@{$order{$col2}{$chr2}->{index}->{1}}){
		
		if(${$order{$col1}{$chr1}->{index}->{1}}[$i] == ${$order{$col2}{$chr2}->{index}->{1}}[$j]){
		    ${$order{$col1}{$chr1}->{index}->{2}}[$i]=$i;
		    ${$order{$col2}{$chr2}->{index}->{2}}[$i]=$j;
		    last;
		}
	    }
	}
	
	foreach my $i (@{$order{$col1}{$chr1}->{index}->{2}}){			#use rank chr1 as an ID
	    foreach my $j (@{$order{$col2}{$chr2}->{index}->{2}}){
		
		if(${$order{$col1}{$chr1}->{index}->{2}}[$i] == ${$order{$col2}{$chr2}->{index}->{2}}[$j]){
		    ${$order{$col1}{$chr1}->{order}}[$i]=$i+1;
		    ${$order{$col2}{$chr2}->{order}}[$i]=$j+1;
		    last;
		}
	    }
	}
	
	@pairs=sortTablebyIndex(\@{$order{$col1}{$chr1}->{index}->{1}},\@pairs);#sorting of the pairs, strands, and start coords from the sorted chr2 reads
	@{$strands{$col1}{$chr1}}=sortTablebyIndex(\@{$order{$col1}{$chr1}->{index}->{1}},$strands{$col1}{$chr1});
	@{$strands{$col2}{$chr2}}=sortTablebyIndex(\@{$order{$col1}{$chr1}->{index}->{1}},$strands{$col2}{$chr2});
	@{$ends_order{$col1}{$chr1}}=sortTablebyIndex(\@{$order{$col1}{$chr1}->{index}->{1}},$ends_order{$col1}{$chr1});
	@{$ends_order{$col2}{$chr2}}=sortTablebyIndex(\@{$order{$col1}{$chr1}->{index}->{1}},$ends_order{$col2}{$chr2});
	@{$coords{$col1}{$chr1}->{start}}=sortTablebyIndex(\@{$order{$col1}{$chr1}->{index}->{1}},$coords{$col1}{$chr1}->{start});
	@{$coords{$col2}{$chr2}->{start}}=sortTablebyIndex(\@{$order{$col1}{$chr1}->{index}->{1}},$coords{$col2}{$chr2}->{start});
	
	
	my @link=($chr->{$chr1}->{name}, $coord_start_chr1 , $coord_end_chr1,	#all information output
		  $chr->{$chr2}->{name}, $coord_start_chr2 , $coord_end_chr2,
		  scalar @pairs,
		  join(",",@pairs),
		  join(",",@{$strands{$col1}{$chr1}}),
		  join(",",@{$strands{$col2}{$chr2}}),
		  join(",",@{$ends_order{$col1}{$chr1}}),
		  join(",",@{$ends_order{$col2}{$chr2}}),
		  join(",",@{$order{$col1}{$chr1}->{order}}),
		  join(",",@{$order{$col2}{$chr2}->{order}}),
		  join(",",@{$coords{$col1}{$chr1}->{start}}),
		  join(",",@{$coords{$col2}{$chr2}->{start}}));
	
	print $fh2 join("\t",@link)."\n";
	
	if($record>=$warn){
	    print LOG "-- $fchr : $warn links processed\n";
	    $warn+=100000;
	}
    }
    $fh->close;
    $fh2->close;
    
    print LOG "-- $fchr : Total : $record links processed\n";
    
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#Sort links according the concerned chromosomes and their coordinates
sub sortLinks{
    
    my ($links_file,$sortedlinks_file,$unique)=@_;
    
    my @sfile=split(/\./,$links_file);
    my $fchr=$sfile[$#sfile-2];
    
    
    print LOG "# $fchr : Sorting links...\n";
    
    my $pipe=($unique)? "| sort -u":"";
    system "sort -k 1,1 -k 4,4 -k 2,2n -k 5,5n -k 8,8n $links_file $pipe > $sortedlinks_file";

}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#removal of fully overlapped links
sub removeFullyOverlappedLinks{
    
    my ($links_file,$nrlinks_file,$warn_out)=@_;
    
    my %pt;
    my $n=1;
    
    my @sfile=split(/\./,$links_file);
    my $fchr=$sfile[$#sfile-2];
    
    my $fh = new FileHandle;
    
    $fh->open("<$links_file") or die "$0: can't open $links_file :$!\n";
    while(<$fh>){
	
	my @t=split("\t",$_);
	$pt{$n}=[@t];
	$n++;
    }
    $fh->close;
    
    my $nb_links=$n-1;
    my $nb=$nb_links;

    my %pt2;
    my $nb2=1;
    my $record=0;
    my $warn=10000;
    
    print LOG "# $fchr : Removing fully overlapped links...\n";
    
    LINK:

    for my $i (1..$nb){
	    
	my @link=();
	my @next_link=();
	my $ind1=$i;
	
	$record++;
	if($record>=$warn){
	    print LOG "-- $fchr : $warn unique links analysed - ".($nb2-1)." non-overlapped links done\n";
	    $warn+=10000;
	}
	
	if(exists $pt{$ind1}){
	    @link=@{$pt{$ind1}};                                                #link1
	}else{
	    next LINK;
	}
	
	my ($chr1,$start1,$end1,$chr2,$start2,$end2)=($link[0],$link[1],$link[2],$link[3],$link[4],$link[5]);	#get info of link1
	my @mates=deleteBadOrderSensePairs(split(",",$link[7]));
    
	my $ind2=$ind1+1;
	$ind2++ while (!exists $pt{$ind2}&& $ind2<=$nb);			#get the next found link
	
	if($ind2<=$nb){
	    
	    @next_link=@{$pt{$ind2}};					        #link2
	    my ($chr3,$start3,$end3,$chr4,$start4,$end4)=($next_link[0],$next_link[1],$next_link[2],$next_link[3],$next_link[4],$next_link[5]); #get info of link2
	    my @next_mates=deleteBadOrderSensePairs(split(",",$next_link[7]));
	  
	    while(($chr1 eq $chr3 && $chr2 eq $chr4) && overlap($start1,$end1,$start3,$end3)){	#loop here according to the chr1 coordinates, need an overlap between links to enter
		    
		if(!overlap($start2,$end2,$start4,$end4)){			#if no overlap with chr2 coordinates ->next link2
		    
		    $ind2++;
		    $ind2++ while (!exists $pt{$ind2}&& $ind2<=$nb);
		   
		    if($ind2>$nb){						#if no more link in the file -> save link1
			
			$pt2{$nb2}=\@link;
			$nb2++;
			next LINK;
		    }
		    
		    @next_link=@{$pt{$ind2}};
		    ($chr3,$start3,$end3,$chr4,$start4,$end4)=($next_link[0],$next_link[1],$next_link[2],$next_link[3],$next_link[4],$next_link[5]);
		    @next_mates=deleteBadOrderSensePairs(split(",",$next_link[7]));
		    next;
		}
		
		my %mates=map{$_ =>1} @mates;					#get the equal number of mates
		my @same_mates = grep( $mates{$_}, @next_mates );
		my $nb_mates= scalar @same_mates;
	      
		if($nb_mates == scalar @mates){
		    
		    delete $pt{$ind1};						#if pairs of link 1 are all included in link 2 -> delete link1
		    next LINK;							#go to link2, link2 becomes link1
		    
		}else{
			delete $pt{$ind2} if($nb_mates == scalar @next_mates);	#if pairs of link2 are all included in link 1 -> delete link2
			$ind2++;							#we continue by checking the next link2
			$ind2++ while (!exists $pt{$ind2}&& $ind2<=$nb);
			
			if($ind2>$nb){						#if no more link in the file -> save link1
			  
			    $pt2{$nb2}=\@link;
			    $nb2++;
			    next LINK;
			}
		    
			@next_link=@{$pt{$ind2}};					#get info of link2
			($chr3,$start3,$end3,$chr4,$start4,$end4)=($next_link[0],$next_link[1],$next_link[2],$next_link[3],$next_link[4],$next_link[5]);
			@next_mates=deleteBadOrderSensePairs(split(",",$next_link[7]));
			
		}
	    }
	}
	$pt2{$nb2}=\@link;							#if no (more) link with chr1 coordinates overlap -> save link1
	$nb2++;
    }																																				
   
    print LOG "-- $fchr : Total : $nb_links unique links analysed - ".($nb2-1)." non-overlapped links done\n";
    
    #OUTPUT
    
    $fh = new FileHandle;
    $fh->open(">$nrlinks_file") or die "$0: can't write in the output: $nrlinks_file :$!\n";
    print LOG "-- $fchr :  writing...\n";
    for my $i (1..$nb2-1){
	
	print $fh join("\t",@{$pt2{$i}});				#all links output
    }
    
    close $fh;
    
    print LOG "-- $fchr : output created: $nrlinks_file\n" if($warn_out);
    
    undef %pt;
    undef %pt2;
}
#------------------------------------------------------------------------------#
sub postFiltering {
    
    my ($links_file,$pflinks_file, $finalScore_thres)=@_;
    
    my @sfile=split(/\./,$links_file);
    my $fchr=$sfile[$#sfile-2];
    
    
    my ($nb,$nb2)=(0,0);
    
    print LOG "# $fchr : Post-filtering links...\n";
    print LOG "-- $fchr : final score threshold = $finalScore_thres\n";
    
    my $fh = new FileHandle;
    my $fh2 = new FileHandle;
    
    $fh->open("<$links_file") or die "$0: can't open $links_file :$!\n";
    $fh2->open(">$pflinks_file") or die "$0: can't write in the output: $pflinks_file :$!\n";
    
    
    while(<$fh>){
	
	my @t=split("\t",$_);
	my $score=$t[$#t-1];
	
	if($score >= $finalScore_thres){
	    print $fh2 join("\t", @t);
	    $nb2++;
	}
	$nb++;
    }
    $fh->close;
    $fh2->close;
    
    print LOG "-- $fchr : Total : $nb unique links analysed - $nb2 links kept\n";
    print LOG "-- $fchr : output created: $pflinks_file\n";
}



#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#Filtering of the links
sub strandFiltering{
    
    my($chr,$chrID,$pairs_threshold,$strand_filtering,$chromosomes,
       $input_format,$cmap_file,$mate_sense, $tag_length,$links_file,$flinks_file)=@_;
    
    my @sfile=split(/\./,$links_file);
    my $fchr=$sfile[$#sfile-1];
    
    
    my %chrs;
    my %chrs1;
    my %chrs2;
    my $nb_chrs;
    my $exclude;
    
    if($chromosomes){
	my @chrs=split(",",$chromosomes);
	$nb_chrs=scalar @chrs;
	$exclude=($chrs[0]=~/^\-/)? 1:0;
	for my $chrName (@chrs){
	    $chrName=~s/^(\-)//;
	    my $col=($chrName=~s/_(1|2)$//);
	    
	    if(!$col){
		$chrs{$chrID->{$chrName}}=undef 
	    }else{
		$chrs1{$chrID->{$chrName}}=undef if($1==1);
		$chrs2{$chrID->{$chrName}}=undef if($1==2);
	    }
	}
    }
 
    my $record=0;
    my $nb_links=0;
    my $warn=10000;
    
    my $sens_ratio_threshold=0.6;
    
    print LOG "\# Filtering procedure...\n";
    print LOG "\# Number of pairs and strand filtering...\n";
    print LOG "-- file=$links_file\n";
    print LOG "-- nb_pairs_threshold=$pairs_threshold, strand_filtering=".(($strand_filtering)? "yes":"no").
	 ", chromosomes=".(($chromosomes)? "$chromosomes":"all")."\n";
    
    
    
    my $fh = new FileHandle;
    my $fh2 = new FileHandle;
    
    $fh->open("<$links_file") or die "$0: can't open $links_file :$!\n";
    $fh2->open(">$flinks_file") or die "$0: can't write in the output: $flinks_file :$!\n";

    while(<$fh>){
	
	my @t=split;	#for each link
	my $is_good=1;
	$record++;
	
	
	if($chromosomes){
	    
	    my ($chr1,$chr2)=($chrID->{$t[0]},$chrID->{$t[3]});
	    
	    if(!$exclude){
		$is_good=(exists $chrs{$chr1} && exists $chrs{$chr2})? 1:0;
		$is_good=(exists $chrs1{$chr1} && exists $chrs2{$chr2})? 1:0 if(!$is_good);
		$is_good=($nb_chrs==1 && (exists $chrs1{$chr1} || exists $chrs2{$chr2}))? 1:0 if(!$is_good);
	    }else{
		$is_good=(exists $chrs{$chr1} || exists $chrs{$chr2})? 0:1;
		$is_good=(exists $chrs1{$chr1} || exists $chrs2{$chr2})? 0:1 if($is_good);
	    }
	}
	
	$is_good = ($is_good && $t[6] >= $pairs_threshold)? 1 :0;		#filtering according the number of pairs
	if($is_good && $strand_filtering){					#if filtering according the strand sense
	    
	    my @mates=split(/,/,$t[7]);						#get the concordant pairs in the strand sense
	    my @strands1=split(/,/,$t[8]);
	    my @strands2=split(/,/,$t[9]);
	    
	    my %mate_class=( 'FF' => 0, 'RR' => 0, 'FR' => 0,  'RF' => 0);
	    
	    my %mate_reverse=( 'FF' => 'RR', 'RR' => 'FF',			#group1: FF,RR
		       'FR' => 'RF', 'RF' => 'FR');				#group2: FR,RF
	    
	    my %mate_class2=( $mate_sense=>"NORMAL_SENSE", inverseSense($mate_sense)=>"NORMAL_SENSE",
		    substr($mate_sense,0,1).inverseSense(substr($mate_sense,1,1))=>"REVERSE_SENSE",
		    inverseSense(substr($mate_sense,0,1)).substr($mate_sense,1,1)=>"REVERSE_SENSE");
	    
	    if($t[6] == 1){
		
		push(@t,$mate_class2{$strands1[0].$strands2[0]},"1/1",1,1);
		
	    }else{
		
		tie (my %class,'Tie::IxHash');
		my $split;
		
		foreach my $i (0..$#mates){
		    $mate_class{$strands1[$i].$strands2[$i]}++;		#get the over-represented group
		}
		
		my $nb_same_sens_class=$mate_class{FF}+$mate_class{RR};
		my $nb_diff_sens_class=$mate_class{FR}+$mate_class{RF};
		my $sens_ratio=max($nb_same_sens_class,$nb_diff_sens_class)/($nb_same_sens_class+$nb_diff_sens_class);
		
		if($sens_ratio < $sens_ratio_threshold){
		    %class=(1=>'FF', 2=>'FR');
		    $split=1;
		}else{
		    $class{1}=($nb_same_sens_class >  $nb_diff_sens_class)? 'FF':'FR';	#if yes get the concerned class
		    $split=0;
		}
		
		$is_good=getConsistentSenseLinks(\@t,\@mates,\@strands1,\@strands2,$tag_length,$mate_sense,\%mate_reverse,\%mate_class2,\%class,$split,$pairs_threshold);
	    }
	}
	
	if($is_good){								#PRINT
	    
	    my $nb=scalar @t;
	    if($nb > 20){
		my @t2=splice(@t,0,20);
		print $fh2 join("\t",@t2)."\n";
		$nb_links++;
	    }
	    $nb_links++;
	    print $fh2 join("\t",@t)."\n";
	}
	
	if($record>=$warn){
	   print LOG "-- $fchr : $warn links analysed - $nb_links links kept\n";
	    $warn+=10000;
	}
    }
    $fh->close;
    $fh2->close;
    
    print LOG "-- $fchr : No links have been found with the selected filtering parameters\n" if(!$nb_links);
    
    print LOG "-- $fchr : Total : $record links analysed - $nb_links links kept\n";
    

}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub getConsistentSenseLinks{
    
    my ($t,$mates,$strands1,$strands2,$tag_length,$mate_sense, $mate_reverse,$mate_class2, $class, $split,$thres)=@_;
    
    my $npairs=scalar @$mates;

    my @ends_order1 = split (/,/,$$t[10]);
    my @ends_order2 = split (/,/,$$t[11]);
    my @order1 = split (/,/,$$t[12]);
    my @order2 = split (/,/,$$t[13]);
    my @positions1 = split (/,/,$$t[14]);
    my @positions2 = split (/,/,$$t[15]);
    
    my @newlink;

    foreach my $ind (keys %{$class} ){
	
	tie (my %flink,'Tie::IxHash');
	my @orders2remove=();
	
	foreach my $i (0..$#{$mates}){					#get the pairs belonging the over-represented group
	
	    if((($$strands1[$i].$$strands2[$i]) eq $$class{$ind}) || (($$strands1[$i].$$strands2[$i]) eq $$mate_reverse{$$class{$ind}})){
		push(@{$flink{mates}},$$mates[$i]);
		push(@{$flink{strands1}},$$strands1[$i]);
		push(@{$flink{strands2}},$$strands2[$i]);
		push(@{$flink{ends_order1}},$ends_order1[$i]);
		push(@{$flink{ends_order2}},$ends_order2[$i]);
		push(@{$flink{positions1}},$positions1[$i]);
		push(@{$flink{positions2}},$positions2[$i]);
			    
	    }else{
		push(@orders2remove,$order1[$i]);
	    }
	}
	
	@{$flink{order1}}=();
	@{$flink{order2}}=();
	if(scalar @orders2remove > 0){
	    getNewOrders(\@order1,\@order2,\@orders2remove,$flink{order1},$flink{order2})
	}else{
	    @{$flink{order1}}=@order1;
	    @{$flink{order2}}=@order2;
	}
	
	my @ends1; getEnds(\@ends1,$flink{positions1},$flink{strands1},$flink{ends_order1},$tag_length);
	my @ends2; getEnds(\@ends2,$flink{positions2},$flink{strands2},$flink{ends_order2},$tag_length);
	
	my $fnpairs=scalar @{$flink{mates}};
	my $strand_filtering_ratio=$fnpairs."/".$npairs;
	my $real_ratio=$fnpairs/$npairs;
	
	if($fnpairs>=$thres){							#filtering according the number of pairs
	    
	    push(@newlink,
		$$t[0],
		min(min(@{$flink{positions1}}),min(@ends1)),
		max(max(@{$flink{positions1}}),max(@ends1)),
		$$t[3],
		min(min(@{$flink{positions2}}),min(@ends2)),
		max(max(@{$flink{positions2}}),max(@ends2)),
		$fnpairs,
		join(",",@{$flink{mates}}),
		join(",",@{$flink{strands1}}),
		join(",",@{$flink{strands2}}),
		join(",",@{$flink{ends_order1}}),
		join(",",@{$flink{ends_order2}}),
		join(",",@{$flink{order1}}),
		join(",",@{$flink{order2}}),
		join(",",@{$flink{positions1}}),
		join(",",@{$flink{positions2}}),
		$$mate_class2{${$flink{strands1}}[0].${$flink{strands2}}[0]},
		$strand_filtering_ratio,
		$real_ratio,
		$npairs
	    );
	}
    }
    
    if (grep {defined($_)} @newlink) {
	@$t=@newlink;
	return 1
    }
    return 0;
    
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub getNewOrders{
    
    my($tab1,$tab2,$list,$newtab1,$newtab2)=@_;
    my $j=1;
    my $k=1;
    for my $i (0..$#{$tab2}){
	my $c=0;
	for my $j (0..$#{$list}){
	    $c++ if(${$list}[$j] < ${$tab2}[$i]);
	    if(${$list}[$j] == ${$tab2}[$i]){
		$c=-1; last;
	    }
	}
	if($c!=-1){
	    push(@{$newtab2}, ${$tab2}[$i]-$c);
	    push(@{$newtab1}, $k);
	    $k++;
	}
    }
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#Filtering of the links using their order
sub orderFiltering {

    my ($chr,$chrID,$nb_pairs_threshold,$nb_pairs_order_threshold,$mu,$sigma,$mate_sense,$tag_length,$links_file,$flinks_file)=@_;
    
    my @sfile=split(/\./,$links_file);
    my $fchr=$sfile[$#sfile-2];
    
    
    my $diff_sense_ends=(($mate_sense eq "FR") || ($mate_sense eq "RF"))? 1:0; 
    
    my $record=0;
    my $warn=10000;
    my $nb_links=0;
    
    my $quant05 = 1.644854;
    my $quant001 = 3.090232;
    my $alphaDist = $quant05 * 2 * $sigma;
    my $maxFragmentLength = &floor($quant001 * $sigma + $mu);
    
    print LOG "\# Filtering by order...\n";
    print LOG "-- mu length=$mu, sigma length=$sigma, nb pairs order threshold=$nb_pairs_order_threshold\n";
    print LOG "-- distance between comparable pairs was set to $alphaDist\n";
    print LOG "-- maximal fragment length was set to $maxFragmentLength\n";

    
    my $fh = new FileHandle;
    my $fh2 = new FileHandle;
    
    $fh->open("<$links_file") or die "$0: can't open $links_file :$!\n";
    $fh2->open(">$flinks_file") or die "$0: can't write in the output: $flinks_file :$!\n";
    
    while(<$fh>){
	
	$record++; 
	my @t = split; 
	my ($chr1,$chr2,$mates_list)=@t[0,3,7];
	my @pairs=split(",",$mates_list);
	($chr1,$chr2) = ($chrID->{$chr1},$chrID->{$chr2});
	my ($coord_start_chr1,$coord_end_chr1,$coord_start_chr2,$coord_end_chr2) = @t[1,2,4,5];
	my $numberOfPairs = $t[6];
	my @strand1 = split (/,/,$t[8]);
	my @strand2 = split (/,/,$t[9]);
	my @ends_order1 = split (/,/,$t[10]);
	my @ends_order2 = split (/,/,$t[11]);
	my @order1 = split (/,/,$t[12]);
	my @order2 = split (/,/,$t[13]);
	my @positions1 = split (/,/,$t[14]); 
	my @positions2 = split (/,/,$t[15]);
	my @ends1; getEnds(\@ends1,\@positions1,\@strand1,\@ends_order1,$tag_length);
	my @ends2; getEnds(\@ends2,\@positions2,\@strand2,\@ends_order2,$tag_length);
	my $clusterCoordinates_chr1;
	my $clusterCoordinates_chr2;
	my $reads_left = 0;
	
	my $ifRenv = $t[16];
	my $strand_ratio_filtering=$t[17]; 
	
	#kind of strand filtering. For example, will keep only FFF-RRR from a link FFRF-RRRF if <F-R> orientation is correct
	my ($singleBreakpoint, %badInFRSense) = findBadInFRSenseSOLiDSolexa(\@strand1,\@strand2,\@ends_order1,\@ends_order2,\@order1,\@order2,$mate_sense);
										#find pairs type F-RRRR or FFFF-R in the case if <R-F> orientation is correct
										#These pairs are annotated as BED pairs forever! They won't be recycled!
	my $table;
	for my $i (0..$numberOfPairs-1) {					#fill the table with non adequate pairs: pairID	numberOfNonAdPairs nonAdPairIDs
	    my $nonAdeq = 0;
	    for my $j (0..$i-1) {
		if (exists($table->{$j}->{$i})) {		    
		    $nonAdeq++;
		    $table->{$i}->{$j} = 1;		    
		}
	    }
	    for my $j ($i+1..$numberOfPairs-1) {
		if ($positions1[$j]-$positions1[$i]>$alphaDist) {		    
		    if (&reversed ($i,$j,$ifRenv,\@positions2)) {
			$nonAdeq++;
			$table->{$i}->{$j} = 1;
		    }
		}
	    }
	    $table->{$i}->{nonAdeq} = $nonAdeq;   	    
	}	
							
	for my $bad (keys %badInFRSense) { #remove pairs type F-RRRR or FFFF-R in the case of <R-F> orientation
	    &remove($bad,$table);	    
	}	
	
	my @falseReads;	
							#RRRR-F -> RRRR or R-FFFF -> FFFF	
	@falseReads = findBadInRFSenseSOLiDSolexa(\@strand1,\@ends_order1,$mate_sense, keys %{$table});
										#these pairs will be recycled later as $secondTable
	for my $bad (@falseReads) {
	    &remove($bad,$table);
	}			
		
	my $bad = &check($table); 
	while ($bad ne "OK") {							#clear the table to reject non adequate pairs in the sense of ORDER
	   # push (@falseReads, $bad);  remove completely!!!
	    &remove($bad,$table);  
	    $bad = &check($table);	    
	}
	
	$reads_left = scalar keys %{$table}; 
	my $coord_start_chr1_cluster1 = min(min(@positions1[sort {$a<=>$b} keys %{$table}]),min(@ends1[sort {$a<=>$b} keys %{$table}]));
	my $coord_end_chr1_cluster1 = max(max(@positions1[sort {$a<=>$b} keys %{$table}]),max(@ends1[sort {$a<=>$b} keys %{$table}]));   
	my $coord_start_chr2_cluster1 = min(min(@positions2[sort {$a<=>$b} keys %{$table}]),min(@ends2[sort {$a<=>$b} keys %{$table}]));
	my $coord_end_chr2_cluster1 = max(max(@positions2[sort {$a<=>$b} keys %{$table}]),max(@ends2[sort {$a<=>$b} keys %{$table}]));  	
	
	$clusterCoordinates_chr1 = '('.$coord_start_chr1_cluster1.','.$coord_end_chr1_cluster1.')';
	$clusterCoordinates_chr2 = '('.$coord_start_chr2_cluster1.','.$coord_end_chr2_cluster1.')';
	
	my $ifBalanced = 'UNBAL';
	my $secondTable;
	my $clusterCoordinates;
	my ($break_pont_chr1,$break_pont_chr2);
	
	my $signatureType="";
	
	my $maxCoord1 =$chr->{$chr1}->{length};
	my $maxCoord2 =$chr->{$chr2}->{length};
	
	if (scalar @falseReads) {
	    @falseReads = sort @falseReads; 
	    #now delete FRFR choosing the majority
	    my @newfalseReads;	   							#find and remove pairs type RRRR-F or R-FFFF	
	    @newfalseReads = findBadInRFSenseSOLiDSolexa(\@strand1,\@ends_order1,$mate_sense,@falseReads);	     #these @newfalseReads  won't be recycled
	    my %hashTmp;
	    for my $count1 (0..scalar(@falseReads)-1) {
		my $i = $falseReads[$count1];
		$hashTmp{$i} = 1;
		for my $bad (@newfalseReads) {
		    if ($bad == $i) {
			delete $hashTmp{$i};
			next;
		    }
		}
	    }
	    @falseReads = sort keys %hashTmp;  #what is left
	    for my $count1 (0..scalar(@falseReads)-1) {				#fill the table for reads which were previously rejected  
		my $nonAdeq = 0;
		my $i = $falseReads[$count1];

		for my $count2 (0..$count1-1) {
		    my $j = $falseReads[$count2];
		    if (exists($secondTable->{$j}->{$i})) {		    
			$nonAdeq++;
			$secondTable->{$i}->{$j} = 1;		    
		    }
		}
		for my $count2 ($count1+1..scalar(@falseReads)-1) {
		    my $j = $falseReads[$count2];
		    if ($positions1[$j]-$positions1[$i]>$alphaDist) {
			if (&reversed ($i,$j,$ifRenv,\@positions2)) {
			    $nonAdeq++;
			    $secondTable->{$i}->{$j} = 1;
			}
		    }
		}
		$secondTable->{$i}->{nonAdeq} = $nonAdeq;
	    }
	    
	    my @falseReads2;
	    my $bad = &check($secondTable);
	    while ($bad ne "OK") {						#clear the table to reject non adequate pairs
		push (@falseReads2, $bad);
		&remove($bad,$secondTable);
		$bad = &check($secondTable);
	    }
	    if (scalar keys %{$secondTable} >= $nb_pairs_order_threshold) {
		my $coord_start_chr1_cluster2 = min(min(@positions1[sort {$a<=>$b} keys %{$secondTable}]),min(@ends1[sort {$a<=>$b} keys %{$secondTable}]));
		my $coord_end_chr1_cluster2 = max(max(@positions1[sort {$a<=>$b} keys %{$secondTable}]),max(@ends1[sort {$a<=>$b} keys %{$secondTable}]));   
		my $coord_start_chr2_cluster2 = min(min(@positions2[sort {$a<=>$b} keys %{$secondTable}]),min(@ends2[sort {$a<=>$b} keys %{$secondTable}]));
		my $coord_end_chr2_cluster2 = max(max(@positions2[sort {$a<=>$b} keys %{$secondTable}]),max(@ends2[sort {$a<=>$b} keys %{$secondTable}]));   		 
		
		$ifBalanced = 'BAL';
		
		if ($ifBalanced eq 'BAL') {
		    
		    if (scalar keys %{$table} < $nb_pairs_order_threshold) {
			$ifBalanced = 'UNBAL';					#kill cluster 1!
			($table,$secondTable)=($secondTable,$table);		#this means that one needs to exchange cluster1 with cluster2
			$reads_left = scalar keys %{$table};
			$coord_start_chr1_cluster1 = $coord_start_chr1_cluster2;
			$coord_end_chr1_cluster1 = $coord_end_chr1_cluster2;   
			$coord_start_chr2_cluster1 = $coord_start_chr2_cluster2;
			$coord_end_chr2_cluster1 = $coord_end_chr2_cluster2;
			$clusterCoordinates_chr1 = '('.$coord_start_chr1_cluster1.','.$coord_end_chr1_cluster1.')';
			$clusterCoordinates_chr2 = '('.$coord_start_chr2_cluster1.','.$coord_end_chr2_cluster1.')';
			
		    } else {
			
			$reads_left += scalar keys %{$secondTable};
			next if ($reads_left < $nb_pairs_threshold);
			
			if ($coord_end_chr1_cluster2 < $coord_start_chr1_cluster1) {
			    ($table,$secondTable)=($secondTable,$table);	 #this means that one needs to exchange cluster1 with cluster2
			    
			    ($coord_start_chr1_cluster1,$coord_start_chr1_cluster2) = ($coord_start_chr1_cluster2,$coord_start_chr1_cluster1);
			    ($coord_end_chr1_cluster1,$coord_end_chr1_cluster2)=($coord_end_chr1_cluster2,$coord_end_chr1_cluster1);   
			    ($coord_start_chr2_cluster1,$coord_start_chr2_cluster2)=($coord_start_chr2_cluster2,$coord_start_chr2_cluster1);
			    ($coord_end_chr2_cluster1 , $coord_end_chr2_cluster2)=($coord_end_chr2_cluster2 , $coord_end_chr2_cluster1);
			    
			    $clusterCoordinates_chr1 = '('.$coord_start_chr1_cluster1.','.$coord_end_chr1_cluster1.'),'.$clusterCoordinates_chr1;
			    $clusterCoordinates_chr2 = '('.$coord_start_chr2_cluster1.','.$coord_end_chr2_cluster1.'),'.$clusterCoordinates_chr2;			
			}else {
			    $clusterCoordinates_chr1 .= ',('.$coord_start_chr1_cluster2.','.$coord_end_chr1_cluster2.')';
			    $clusterCoordinates_chr2 .= ',('.$coord_start_chr2_cluster2.','.$coord_end_chr2_cluster2.')';			
			}
			$coord_start_chr1 = min($coord_start_chr1_cluster1,$coord_start_chr1_cluster2);
			$coord_end_chr1 = max($coord_end_chr1_cluster1,$coord_end_chr1_cluster2);
			$coord_start_chr2 =  min($coord_start_chr2_cluster1,$coord_start_chr2_cluster2);
			$coord_end_chr2 = max($coord_end_chr2_cluster1,$coord_end_chr2_cluster2);
			#to calculate breakpoints one need to take into account read orientation in claster..
			my $leftLetterOk = substr($mate_sense, 0, 1);  #R
			my $rightLetterOk = substr($mate_sense, 1, 1); #F
			
			
			my @index1 = keys %{$table};
			my @index2 = keys %{$secondTable};
			
			my (@generalStrand1,@generalStrand2) = 0;
			
			if ($leftLetterOk eq $rightLetterOk) { #SOLID mate-pairs
			    $leftLetterOk = 'R';
			    $rightLetterOk = 'F';
			    @generalStrand1 = translateSolidToRF(\@strand1,\@ends_order1);	
			    @generalStrand2 = translateSolidToRF(\@strand2,\@ends_order2);
			} else {
			    @generalStrand1 = @strand1;
			    @generalStrand2 = @strand2; # TODO check if it is correct
			}
			if ($generalStrand1[$index1[0]] eq $leftLetterOk && $generalStrand1[$index2[0]] eq $rightLetterOk) { #(R,F)
			    $break_pont_chr1 = '('.$coord_end_chr1_cluster1.','.$coord_start_chr1_cluster2.')';
			    
			    if ($generalStrand2[$index1[0]] eq $rightLetterOk && $generalStrand2[$index2[0]] eq $leftLetterOk) {
				if ($coord_end_chr2_cluster1 >= $coord_end_chr2_cluster2) {
				    $break_pont_chr2 = '('.$coord_end_chr2_cluster2.','.$coord_start_chr2_cluster1.')';
				    $signatureType = "TRANSLOC";
				} else {
				    $break_pont_chr2 = '('.max(($coord_end_chr2_cluster1-$maxFragmentLength),1).','.$coord_start_chr2_cluster1.')';	
				    $break_pont_chr2 .= ',('.$coord_end_chr2_cluster2.','.min(($coord_start_chr2_cluster2+$maxFragmentLength),$maxCoord2).')';
				    $signatureType = "INS_FRAGMT";
				}		    
				
			    } elsif ($generalStrand2[$index1[0]] eq $leftLetterOk && $generalStrand2[$index2[0]] eq $rightLetterOk) {
				if ($coord_end_chr2_cluster1 >= $coord_end_chr2_cluster2) {
				    $break_pont_chr2 = '('.max(($coord_end_chr2_cluster2-$maxFragmentLength),1).','.$coord_start_chr2_cluster2.')';	
				    $break_pont_chr2 .= ',('.$coord_end_chr2_cluster1.','.min(($coord_start_chr2_cluster1+$maxFragmentLength),$maxCoord2).')';
				    $signatureType = "INV_INS_FRAGMT";
				} else {
				    $break_pont_chr2 = '('.$coord_end_chr2_cluster1.','.$coord_start_chr2_cluster2.')';
				    $signatureType = "INV_TRANSLOC";
				}	
			    } else {
				#should not occur
				print STDERR "\nError in orderFiltering\n\n";				    
			    }
			}
			
			elsif ($generalStrand1[$index1[0]] eq $rightLetterOk && $generalStrand1[$index2[0]] eq $leftLetterOk) { #(F,R)
			    $break_pont_chr1 = '('.max(($coord_end_chr1_cluster1-$maxFragmentLength),1).','.$coord_start_chr1_cluster1.')';	
			    $break_pont_chr1 .= ',('.$coord_end_chr1_cluster2.','.min(($coord_start_chr1_cluster2+$maxFragmentLength),$maxCoord1).')';				
			    if ($generalStrand2[$index1[0]] eq $rightLetterOk && $generalStrand2[$index2[0]] eq $leftLetterOk) { 
				if ($coord_end_chr2_cluster1 >= $coord_end_chr2_cluster2) {
				    $break_pont_chr2 = '('.$coord_end_chr2_cluster2.','.$coord_start_chr2_cluster1.')';
				    $signatureType = "INV_INS_FRAGMT";
				} else {
				    $break_pont_chr2 = '('.max(($coord_end_chr2_cluster1-$maxFragmentLength),1).','.$coord_start_chr2_cluster1.')';	
				    $break_pont_chr2 .= ',('.$coord_end_chr2_cluster2.','.min(($coord_start_chr2_cluster2+$maxFragmentLength),$maxCoord2).')';
				    $signatureType = "INV_COAMPLICON";
				}		    
				
			    } elsif ($generalStrand2[$index1[0]] eq $leftLetterOk && $generalStrand2[$index2[0]] eq $rightLetterOk) {
				if ($coord_end_chr2_cluster1 >= $coord_end_chr2_cluster2) {
				    $break_pont_chr2 = '('.max(($coord_end_chr2_cluster2-$maxFragmentLength),1).','.$coord_start_chr2_cluster2.')';	
				    $break_pont_chr2 .= ',('.$coord_end_chr2_cluster1.','.min(($coord_start_chr2_cluster1+$maxFragmentLength),$maxCoord2).')';
				    $signatureType = "COAMPLICON";
				} else {
				    $break_pont_chr2 = '('.$coord_end_chr2_cluster1.','.$coord_start_chr2_cluster2.')';
				    $signatureType = "INS_FRAGMT";
				}	
			    } else {
				#should not occur
				$signatureType = "UNDEFINED";
			    }
			}			    
			else { # (F,F) or (R,R)  something strange. We will discard the smallest cluster
			    $ifBalanced = 'UNBAL';				
			    if (scalar keys %{$secondTable} > scalar keys %{$table}) {
				($table,$secondTable)=($secondTable,$table);		#this means that one needs to exchange cluster1 with cluster2
				
				$coord_start_chr1_cluster1 = $coord_start_chr1_cluster2;
				$coord_end_chr1_cluster1 = $coord_end_chr1_cluster2;   
				$coord_start_chr2_cluster1 = $coord_start_chr2_cluster2;
				$coord_end_chr2_cluster1 = $coord_end_chr2_cluster2;
				$clusterCoordinates_chr1 = '('.$coord_start_chr1_cluster1.','.$coord_end_chr1_cluster1.')';
				$clusterCoordinates_chr2 = '('.$coord_start_chr2_cluster1.','.$coord_end_chr2_cluster1.')';			    
			    }
			    $reads_left = scalar keys %{$table};
			}
			if ($ifBalanced eq 'BAL') {
			    $ifRenv = $signatureType;
			}	
		    }
		}
	    }   
	}	
	if ($ifBalanced ne 'BAL') {
	    #define possible break point
	    $coord_start_chr1 = $coord_start_chr1_cluster1;
	    $coord_end_chr1 = $coord_end_chr1_cluster1;
	    $coord_start_chr2 = $coord_start_chr2_cluster1;
	    $coord_end_chr2 = $coord_end_chr2_cluster1;
	    
	    my $region_length_chr1 = $coord_end_chr1-$coord_start_chr1;
	    my $region_length_chr2 = $coord_end_chr2-$coord_start_chr2;

	    my $leftLetterOk = substr($mate_sense, 0, 1);  #R
	    my $rightLetterOk = substr($mate_sense, 1, 1); #F
	    
	    my @index = keys %{$table};			    
	    unless ($diff_sense_ends) {
		my $firstEndOrder1 = $ends_order1[$index[0]];
		my $firstEndOrder2 = $ends_order2[$index[0]];
	        $break_pont_chr1 = (($strand1[$index[0]] eq 'R' && $firstEndOrder1 == 2) || ($strand1[$index[0]] eq 'F' && $firstEndOrder1 == 1))?'('.$coord_end_chr1.','.min(($coord_start_chr1+$maxFragmentLength),$maxCoord1).')':'('.max(($coord_end_chr1-$maxFragmentLength),1).','.$coord_start_chr1.')';
	        $break_pont_chr2 = (($strand2[$index[0]] eq 'R' && $firstEndOrder2 == 2) || ($strand2[$index[0]] eq 'F' && $firstEndOrder2 == 1))?'('.$coord_end_chr2.','.min(($coord_start_chr2+$maxFragmentLength),$maxCoord2).')':'('.max(($coord_end_chr2-$maxFragmentLength),1).','.$coord_start_chr2.')';	
	    } else {
	        $break_pont_chr1 = ($strand1[$index[0]] eq $leftLetterOk )?'('.$coord_end_chr1.','.min(($coord_start_chr1+$maxFragmentLength),$maxCoord1).')':'('.max(($coord_end_chr1-$maxFragmentLength),1).','.$coord_start_chr1.')';
	        $break_pont_chr2 = ($strand2[$index[0]] eq $leftLetterOk )?'('.$coord_end_chr2.','.min(($coord_start_chr2+$maxFragmentLength),$maxCoord2).')':'('.max(($coord_end_chr2-$maxFragmentLength),1).','.$coord_start_chr2.')';	
	    }
	    
	    if ($chr1 ne $chr2){
		$ifRenv="INV_TRANSLOC" if($ifRenv eq "REVERSE_SENSE");
		$ifRenv="TRANSLOC" if($ifRenv eq "NORMAL_SENSE");
	    }  
	}
	
	if (($ifBalanced eq 'BAL')&&( (scalar keys %{$table}) + (scalar keys %{$secondTable}) < $nb_pairs_threshold)) {
	    next; 								#discard the link
	}	
 	if (($ifBalanced eq 'UNBAL')&&(scalar keys %{$table} < $nb_pairs_threshold)) {
	    next; 								#discard the link
	}
	my $ratioTxt = "$reads_left/".(scalar @pairs);
	my ($n1,$nTot) = split ("/",$strand_ratio_filtering);
	my $ratioReal = $reads_left/$nTot;    
	
	if ($coord_start_chr1<=0) {
	    $coord_start_chr1=1;
	}
	if ($coord_start_chr2<=0) {
	    $coord_start_chr2=1;
	}
	#create output	
	my @link=($chr->{$chr1}->{name}, $coord_start_chr1 , $coord_end_chr1,	#all information output
		  $chr->{$chr2}->{name}, $coord_start_chr2 , $coord_end_chr2,
		  $reads_left,
		  &redraw(1,$table,$secondTable,\%badInFRSense,$ifBalanced,\@pairs),
		  &redraw(1,$table,$secondTable,\%badInFRSense,$ifBalanced,\@strand1),
		  &redraw(1,$table,$secondTable,\%badInFRSense,$ifBalanced,\@strand2),
		  &redraw(1,$table,$secondTable,\%badInFRSense,$ifBalanced,\@ends_order1),
		  &redraw(1,$table,$secondTable,\%badInFRSense,$ifBalanced,\@ends_order2),
		  &redraw(2,$table,$secondTable,\%badInFRSense,$ifBalanced,\@order1),
		  &redraw(2,$table,$secondTable,\%badInFRSense,$ifBalanced,\@order2),
		  &redraw(1,$table,$secondTable,\%badInFRSense,$ifBalanced,\@positions1),
		  &redraw(1,$table,$secondTable,\%badInFRSense,$ifBalanced,\@positions2),
		  $ifRenv,
		  $strand_ratio_filtering,
		  $ifBalanced, $ratioTxt, $break_pont_chr1, $break_pont_chr2,
		  $ratioReal, $nTot);
	
	$nb_links++;
	print $fh2 join("\t",@link)."\n";
	
	if($record>=$warn){
	    print LOG "-- $fchr : $warn links analysed - $nb_links links kept\n";
	    $warn+=10000;
	}
	
    }
    $fh->close;
    $fh2->close;
    
    print LOG "-- $fchr : Total : $record links analysed - $nb_links links kept\n";

}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#gets information about ends positions given start, direction and order
sub getEnds {
    my ($ends,$starts,$strand,$end_order,$tag_length) = @_;
    for my $i (0..scalar(@{$starts})-1) {
	$ends->[$i] = getEnd($starts->[$i],$strand->[$i],$end_order->[$i],$tag_length);
    }    
}
sub getEnd {
    my ($start,$strand, $end_order,$tag_length) = @_;       
    return ($strand eq 'F')? $start+$tag_length->{$end_order}-1:$start-$tag_length->{$end_order}+1;
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
sub addInsertionInfo { 								#add field with INS,DEL,NA and distance between clusters and performs filtering 
    
    my ($chr,$chrID,$nb_pairs_threshold,$order_filtering,$indel_sigma_threshold,$dup_sigma_threshold,$singleton_sigma_threshold,$mu,$sigma,$mate_sense,$tag_length,$links_file,$flinks_file)=@_;

    my @sfile=split(/\./,$links_file);
    my $fchr=$sfile[$#sfile-2];
    

    my $diff_sense_ends=(($mate_sense eq "FR") || ($mate_sense eq "RF"))? 1:0; 
    
    my $record=0;
    my $nb_links=0;
    my $warn=10000;
    
    print LOG "\# Filtering out normal pairs using insert size...\n";
    print LOG "-- mu length=$mu, sigma length=$sigma, indel sigma threshold=$indel_sigma_threshold, dup sigma threshold=$dup_sigma_threshold\n";
    print LOG "-- using ".($mu-$indel_sigma_threshold*$sigma)."-".
		      ($mu+$indel_sigma_threshold*$sigma)." as normal range of insert size for indels\n";
    print LOG "-- using ".($mu-$dup_sigma_threshold*$sigma)."-".
		      ($mu+$dup_sigma_threshold*$sigma)." as normal range of insert size for duplications\n";
    print LOG "-- using ".($mu-$singleton_sigma_threshold*$sigma)." as the upper limit of insert size for singletons\n" if($mate_sense eq "RF");

    my $fh = new FileHandle;
    my $fh2 = new FileHandle;
    
    $fh->open("<$links_file") or die "$0: can't open $links_file :$!\n";
    $fh2->open(">$flinks_file") or die "$0: can't write in the output: $flinks_file :$!\n";
    
    while(<$fh>){
	
	$record++; 
	my @t = split;
	my ($chr1,$chr2,$mates_list)=@t[0,3,7];
	
	if($chrID->{$chr1} ne $chrID->{$chr2}) {				#if inter-chromosomal link here (because sv_type=all), 
	    $nb_links++;
	    
	    $t[16]="INV_TRANSLOC" if($t[16] eq "REVERSE_SENSE");
	    $t[16]="TRANSLOC" if($t[16] eq "NORMAL_SENSE");
	    
	    $t[16].= "\t";
	    $t[19].= "\t";
	    
	    print $fh2 join("\t",@t)."\n";
	    
	    if($record>=$warn){
		print LOG "-- $fchr : $warn links processed - $nb_links links kept\n";
		$warn+=10000;
	    }								
	    next;								
	}									
	
	my $ifRenv = $t[16];
	my $ifBalanced = "UNBAL";
	$ifBalanced = $t[18] if ($order_filtering);
		
	my $numberOfPairs = $t[6];
	my @positions1 = deleteBadOrderSensePairs(split (/,/,$t[14]));
	my @positions2 = deleteBadOrderSensePairs(split (/,/,$t[15]));
	
	if ($ifBalanced eq "BAL") {
	    
	    if ($ifRenv eq "INV_TRANSLOC") {
		$ifRenv = "INV_FRAGMT"; #for intrachromosomal inverted translocation is the same as inverted fragment
	    }
	    if ($ifRenv eq "NORMAL_SENSE") {
		$ifRenv = "TRANSLOC"; 
	    }
	    if ($ifRenv eq "REVERSE_SENSE") {
		$ifRenv = "INV_FRAGMT"; #for intrachromosomal inverted translocation is the same as inverted fragment
	    }
	    $t[19].= "\t";
	    
	    my $meanDistance = 0;
	    
	    for my $i (0..$numberOfPairs-1) {
		$meanDistance += $positions2[$i]-$positions1[$i];
	    }
	    $meanDistance /= $numberOfPairs;      
	    
	    $t[16] = $ifRenv."\t".$meanDistance;
	    #dont touch the annotation. It should be already OK.
	    
	} else {
	    #only for unbalanced
	    
	    my $ifoverlap=overlap($t[1],$t[2],$t[4],$t[5]);
	    
	    my $ends_sense_class = (deleteBadOrderSensePairs(split (/,/,$t[8])))[0].
				   (deleteBadOrderSensePairs(split (/,/,$t[9])))[0];
	    my $ends_order_class = (deleteBadOrderSensePairs(split (/,/,$t[10])))[0].
				   (deleteBadOrderSensePairs(split (/,/,$t[11])))[0];
	    
	    my $indel_type = $ifRenv;
	    
	    my $meanDistance = "N/A";
	    
	    ($meanDistance, $indel_type) = checkIndel ($numberOfPairs,		#identify insertion type for rearrangments without inversion, calculates distance between cluster
						       \@positions1,		#assign N/A to $indel_type if unknown
						       \@positions2,
						       $ifRenv,
						       $ifoverlap,
						       $indel_sigma_threshold,
						       $dup_sigma_threshold,
						       $singleton_sigma_threshold,
						       $mu,
						       $sigma,
						       $ifBalanced,
						       $ends_sense_class,
						       $ends_order_class,
						       $mate_sense,
						       $diff_sense_ends,
						       );
	    
	    #filtering of pairs with distance inconsistant with the SV
	    if ($ifRenv ne "REVERSE_SENSE") {
		my $maxCoord1 =$chr->{$chrID->{$chr1}}->{length};	
		my $maxCoord2 =$chr->{$chrID->{$chr2}}->{length};
		$meanDistance = recalc_t_usingInsertSizeInfo(\@t,$mu,$sigma,$meanDistance,$tag_length,$diff_sense_ends,$mate_sense,
							     $maxCoord1,$maxCoord2,$ends_sense_class,$ends_order_class,$nb_pairs_threshold,$order_filtering);
		next if ($t[6] < $nb_pairs_threshold);
	    }else{
		 $t[19].= "\t";
	    }
	    $t[16] = $indel_type."\t".$meanDistance;
	}
	
	$nb_links++;
		
	print $fh2 join("\t",@t)."\n";
	if($record>=$warn){
	     print LOG "-- $fchr : $warn links processed - $nb_links links kept\n";
	    $warn+=10000;
	}
    }
    $fh->close;
    $fh2->close; 
    
    print LOG "-- $fchr : Total : $record links analysed - $nb_links links kept\n";
    
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub checkIndel  {
    
    my ($numberOfPairs,	$positions1, $positions2, $ifRenv, $ifoverlap, $indel_sigma_threshold, $dup_sigma_threshold, $singleton_sigma_threshold,
	$mu, $sigma, $ifBalanced,$ends_sense_class,$ends_order_class,$mate_sense,$diff_sense_ends) = @_;
    
    my $meanDistance = 0;
    
    for my $i (0..$numberOfPairs-1) {
	$meanDistance += $positions2->[$i]-$positions1->[$i];
    }
    $meanDistance /= $numberOfPairs;
    
    return ($meanDistance,"INV_DUPLI") if (($ifRenv eq "REVERSE_SENSE") && ($meanDistance<$mu+$dup_sigma_threshold*$sigma) );
    
    return ($meanDistance,"INVERSION") if ($ifRenv eq "REVERSE_SENSE");
    
    if($diff_sense_ends){
	return ($meanDistance, "LARGE_DUPLI") if ($ends_sense_class ne $mate_sense) && ($meanDistance>$mu+$dup_sigma_threshold*$sigma) ;
	return ($meanDistance, "SINGLETON") if (($meanDistance<$mu-$singleton_sigma_threshold*$sigma) && $mate_sense eq "RF" && ($ends_sense_class eq inverseSense($mate_sense)));
    }else{
	return ($meanDistance, "LARGE_DUPLI") if (($ends_sense_class eq $mate_sense) && ($ends_order_class eq "12") || ($ends_sense_class eq inverseSense($mate_sense)) && ($ends_order_class eq "21")) &&
	($meanDistance>$mu+$dup_sigma_threshold*$sigma) ; 
    }
 
    return ($meanDistance, "SMALL_DUPLI") if (($meanDistance<$mu-$dup_sigma_threshold*$sigma) && $ifoverlap);
    
    return ($meanDistance, "DUPLICATION") if ($diff_sense_ends && ($ends_sense_class ne $mate_sense) && ($meanDistance<$mu-$dup_sigma_threshold*$sigma) ) ;
    
    return ($meanDistance, "INSERTION") if ($meanDistance<$mu -$indel_sigma_threshold*$sigma);
    return ($meanDistance, "DELETION") if ($meanDistance>$mu+$indel_sigma_threshold*$sigma);

    return ($meanDistance, "UNDEFINED");
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#sub reacalulate @t so that get rid of unconsistent pairs (unconsistent insert size )
sub recalc_t_usingInsertSizeInfo {
    my($t,$mu,$sigma,$meanDistance,$tag_length,$diff_sense_ends,$mate_sense,$maxCoord1,$maxCoord2,$ends_sense_class,$ends_order_class,$nb_pairs_threshold,$order_filtering) = @_;
 	    
    my @badPairs;
	    
    my @positions1 = getAllEntries($t->[14]);
    my @positions2 = getAllEntries($t->[15]);	    
	    
    if ($meanDistance < $mu) {
	for my $i (0..scalar(@positions1)-1) {
	    if (substr($positions2[$i],-1,1) ne '$' && substr($positions2[$i],-1,1) ne '*' && $positions2[$i]-$positions1[$i]>=$mu) {
	    	push(@badPairs,$i);		    
	    }
    	}
    } else {
	for my $i (0..scalar(@positions1)-1) {
	    if (substr($positions2[$i],-1,1) ne '$' && substr($positions2[$i],-1,1) ne '*' && $positions2[$i]-$positions1[$i]<=$mu) {
		push(@badPairs,$i);
	    }
	}	    
    }

    if (scalar (@badPairs)>0) {
	#print join("\t",@badPairs).": ".join("\t",@t)."\n";	        
	#remove these inconsistant links
	$t->[6] -= scalar(@badPairs); #numberOfPairs
	return if ($t->[6] < $nb_pairs_threshold); 	

	$t->[7] = mark_values(\@badPairs, $t->[7]);
	$t->[8] = mark_values(\@badPairs, $t->[8]);
	$t->[9] = mark_values(\@badPairs, $t->[9]);
	$t->[10] = mark_values(\@badPairs, $t->[10]);
	$t->[11] = mark_values(\@badPairs, $t->[11]);
	
	$t->[12] = mark_indexes(\@badPairs, $t->[12]);
	$t->[13] = mark_indexes(\@badPairs, $t->[13]);	    
		    
	$t->[14] = mark_values(\@badPairs, $t->[14]);
	$t->[15] = mark_values(\@badPairs, $t->[15]);
	$t->[19] = recalculate_ratio($t->[6],$t->[19]) if ($order_filtering); #add the second ratio
	$t->[17] = recalculate_ratio($t->[6],$t->[17]) unless ($order_filtering);
	($t->[1],$t->[2]) = recalculate_boundaries($t->[14],$t->[8],$t->[10],$tag_length);
	($t->[4],$t->[5]) = recalculate_boundaries($t->[15],$t->[9],$t->[11],$tag_length);	 
		    
	#recalc breakpoints:
	my $quant001 = 3.090232;
	my $maxFragmentLength = &floor($quant001 * $sigma + $mu);
	$t->[20] = recalc_breakpoints($mate_sense,$maxCoord1,$t->[14],substr($ends_sense_class,0,1),substr($ends_order_class,0,1),$t->[1],$t->[2],$maxFragmentLength,$diff_sense_ends ) if ($order_filtering);
	$t->[21] = recalc_breakpoints($mate_sense,$maxCoord2,$t->[15],substr($ends_sense_class,1,1),substr($ends_order_class,1,1),$t->[4],$t->[5],$maxFragmentLength,$diff_sense_ends ) if ($order_filtering);
	#recalc total ratio
	$t->[22] = $t->[6] / $t->[23] if ($order_filtering);
	$t->[18] = $t->[6] / $t->[19] unless ($order_filtering);
	
    	@positions1 = deleteBadOrderSensePairs(split (/,/,$t->[14]));
    	@positions2 = deleteBadOrderSensePairs(split (/,/,$t->[15]));		
	
	$meanDistance = 0;
	
	for my $i (0..scalar(@positions1)-1) {
	    $meanDistance += $positions2[$i]-$positions1[$i];
	}
	$meanDistance /= scalar(@positions1);
 	
    } else {
	$t->[17] = recalculate_ratio((split(/\//,$t->[17]))[0],$t->[17]) unless ($order_filtering);
	$t->[19] = recalculate_ratio((split(/\//,$t->[19]))[0],$t->[19]) if ($order_filtering);
    
    } #nothing has been filtered
    return $meanDistance;
}

sub recalculate_ratio {
    my ($left, $ratio) = @_;
    my @elements = split (/\//,$ratio);
    $elements[1]= $elements[0];
    $elements[0]=$left;
    return $ratio."\t".join("/",@elements);    
}

sub recalc_breakpoints {
    my ($mate_sense,$maxCoord,$startString,$strand,$firstEndOrder,$coord_start_chr,$coord_end_chr,$maxFragmentLength,$diff_sense_ends ) = @_;
    my $break_pont_chr;
    
    my $leftLetterOk = substr($mate_sense, 0, 1);  #R
    my $rightLetterOk = substr($mate_sense, 1, 1); #F	    

	
    my @positions = deleteBadOrderSensePairs(split (/,/,$startString));

    unless ($diff_sense_ends) {
	$break_pont_chr = (($strand eq 'R' && $firstEndOrder == 2) || ($strand eq 'F' && $firstEndOrder == 1))?'('.$coord_end_chr.','.min(($coord_start_chr+$maxFragmentLength),$maxCoord).')':'('.max(($coord_end_chr-$maxFragmentLength),1).','.$coord_start_chr.')';
    } else {
	$break_pont_chr = ($strand eq $leftLetterOk)?'('.$coord_end_chr.','.min(($coord_start_chr+$maxFragmentLength),$maxCoord).')':'('.max(($coord_end_chr-$maxFragmentLength),1).','.$coord_start_chr.')';
    }    
    return $break_pont_chr;
}
sub recalculate_boundaries {
    my ($startString,$senseString,$endsOrderString,$tag_length) = @_;
    my @positions = deleteBadOrderSensePairs(split (/,/,$startString));
    my @strands = deleteBadOrderSensePairs(split (/,/,$senseString));
    my @ends_orders = deleteBadOrderSensePairs(split (/,/,$endsOrderString));
    my @ends; getEnds(\@ends,\@positions,\@strands,\@ends_orders,$tag_length);   
    my $coord_start_cluster = min(min(@positions),min(@ends));
    my $coord_end_cluster = max(max(@positions),max(@ends));
    return ($coord_start_cluster,$coord_end_cluster);
}

sub remove_indexes {
    my ($bads, $string) = @_;
    my @elements = deleteBadOrderSensePairs(split (/,/,$string));
    for my $i (reverse sort %{$bads}) {
	delete $elements[$i];
    }
    return "(".join(",",@elements).")";
}
##add @ to to elements 
sub mark_values {
    my ($bads, $string) = @_;
    my @elements = getAllEntries($string);
    for my $i (@{$bads}) {
	$elements[$i] .= "@";
    }    
    return "(".join(",",@elements).")";
}
##add @ to to indexes
sub mark_indexes {
    my ($bads, $string) = @_;
    my @elements = getAllEntries($string);
    for my $i ((0..scalar(@elements)-1)) {
	for my $j (@{$bads}) {	    
	    $elements[$i] .= "@" if ($elements[$i] eq ($j+1));
	}	   
    }

    return "(".join(",",@elements).")";    
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub redraw {
    
    my ($type,$table,$secondTable,$badInFRSense,$ifBalanced,$arr) = @_;	

    my $out;
    my @first_arr;
    if ($ifBalanced eq 'BAL') {
	my @second_arr;
	my $lastPushed = 1;
	if ($type == 1) {
	    for my $i (0 .. scalar(@{$arr})-1) {
		if (exists ($table->{$i})) {
	    	    push(@first_arr,$arr->[$i]);
		    $lastPushed = 1;
		}elsif (exists ($secondTable->{$i})) {
		    push(@second_arr,$arr->[$i]);
		    $lastPushed = 2;
		} elsif ($lastPushed == 1) {
		    if (exists ($badInFRSense->{$i})) {
		        push(@first_arr,$arr->[$i]."\$");
		    }else {
		        push(@first_arr,$arr->[$i]."*");
		    }
		} elsif ($lastPushed == 2) {
		    if (exists ($badInFRSense->{$i})) {
		        push(@second_arr,$arr->[$i]."\$");
		    }else {
			push(@second_arr,$arr->[$i]."*");
		    }
		} else {print "Error!";exit;}
	    }
	} else {
	     for my $i (@{$arr}) {
		if (exists ($table->{$i-1})) {
	    	    push(@first_arr,$i);
		    $lastPushed = 1;
		}elsif (exists ($secondTable->{$i-1})) {
		    push(@second_arr,$i);
		    $lastPushed = 2;
		} elsif ($lastPushed == 1) {
		    if (exists ($badInFRSense->{$i-1})) {
		        push(@first_arr,$i."\$");
		    }else {
			push(@first_arr,$i."*");
		    }
		} elsif ($lastPushed == 2) {
		    if (exists ($badInFRSense->{$i-1})) {
		        push(@second_arr,$i."\$");
		    }else {
			push(@second_arr,$i."*");
		    }
		} else {print "Error!";exit;}
	    }
	}		
	$out = '('.join(",",@first_arr).'),('.join(",",@second_arr).')';
    }
    else {
	if ($type == 1) {
	    for my $i (0 .. scalar(@{$arr})-1) { 
		if (exists ($table->{$i})) {
	    	    push(@first_arr,$arr->[$i]);
		} else {
		    if (exists ($badInFRSense->{$i})) {
			push(@first_arr,$arr->[$i]."\$");
		    }else {
		        push(@first_arr,$arr->[$i]."*");
		    }
		}
	    }
	} else {
	    for my $i (@{$arr}) {
		if (exists ($table->{$i-1})) {
	    	    push(@first_arr,$i);
		} else {
		    if (exists ($badInFRSense->{$i-1})) {
			push(@first_arr,$i."\$");
		    }else {
			push(@first_arr,$i."*");
		    }
		}
	    }	    
	}	
	$out = '('.join(",",@first_arr).')';
    }    
    return $out;
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub check {
    
    my $table = $_[0];
    my $bad = 'OK';
    my $max = 0;  
    for my $i (sort {$a<=>$b} keys %{$table}) { 				
	unless ($table->{$i}->{nonAdeq} == 0) {
	    if ($max<$table->{$i}->{nonAdeq}) {
		$max=$table->{$i}->{nonAdeq};
		$bad = $i;
	    }
	}	    
    } 
    return $bad;
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub reversed {
    
    my ($i,$j,$ifRenv,$positions) = @_;
    if (($ifRenv eq 'REVERSE_SENSE' && $positions->[$i]<$positions->[$j]) || ($ifRenv ne 'REVERSE_SENSE' && $positions->[$i]>$positions->[$j])){
	return 1;
    }
    return 0;    
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub remove {
    
    my ($bad,$table) = @_;
    for my $i (sort {$a<=>$b} keys %{$table}) {
	if ($bad == $i) {
	    delete($table->{$i});;
	} else {
	    if (exists($table->{$i}->{$bad})) {
		delete($table->{$i}->{$bad});
		$table->{$i}->{nonAdeq}--;
	    }
	}
    }
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub findBadInRFSenseSOLiDSolexa { #choose maximum: FFFFs or RRRRs
    
    my ($strand,$ends_order,$mate_sense,@keysLeft) = @_;
    
    my $leftLetterOk = substr($mate_sense, 0, 1);  #R
    my $rightLetterOk = substr($mate_sense, 1, 1); #F
    
    my (@standardArray);
    if ($leftLetterOk eq $rightLetterOk) { #SOLID mate-pairs
	$leftLetterOk = 'R';
	$rightLetterOk = 'F';
	@standardArray = translateSolidToRF($strand,$ends_order);
    } else {
	@standardArray = @{$strand};	
    }
   
    my $ifR = 0;
    my @Rs;
    
    for my $i (@keysLeft) {	
	if ($standardArray[$i] eq $leftLetterOk) {
	    $ifR++;
	    push(@Rs,$i);
	}
    }
 

    my $ifF = 0;
    my @Fs;
    
    for my $i (@keysLeft) {	
	if ($standardArray[$i] eq $rightLetterOk) {
	    $ifF++;
	    push(@Fs,$i);
	}
    }
    
    if($ifR>=$ifF) {
	return @Fs;
    }
    return @Rs;     
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub findBadInFRSenseSOLiDSolexa { #should work both for SOLiD and Solexa 
    
    my ($strand1,$strand2,$ends_order1,$ends_order2,$order1,$order2) = ($_[0],$_[1],$_[2],$_[3],$_[4],$_[5]);
    my $mate_sense =  $_[6];
      
    my $leftLetterOk = substr($mate_sense, 0, 1);  #R
    my $rightLetterOk = substr($mate_sense, 1, 1); #F
    
    my (@standardArray1,@standardArray2);
   
    if ($leftLetterOk eq $rightLetterOk) { #SOLID mate-pairs
	$leftLetterOk = 'R';
	$rightLetterOk = 'F';
	@standardArray1 = translateSolidToRF($strand1,$ends_order1);
	my @arr = getOrderedStrands($strand2,$order2); 
	my @ends2 = getOrderedStrands($ends_order2,$order2); 
	@standardArray2 = translateSolidToRF(\@arr,\@ends2);  

    } else {
	@standardArray1 = @{$strand1};
	@standardArray2 = getOrderedStrands($strand2,$order2);
    }
   
    #we will try 4 possibilities, 2 for each end of the link: RFRR-FFF->RFFFF , RFRR-FFF->RRRFFF
   
    #for the first end:   
   
    my @array = @standardArray1; 
    my %badInFRSense1;
    for my $i (1..scalar (@array)-1){ # FRFRFFFF -> FFFFFF and RRFRFRFFFF -> RRFFFFFF
	if ($array[$i-1] eq $rightLetterOk && $array[$i] eq $leftLetterOk) {
	    $badInFRSense1{$i}=1;
	    $array[$i] = $rightLetterOk;
	}
    }    
    my $numberRRRFFF_or_FFF_1 = scalar(@array)-scalar(keys %badInFRSense1);
    @array = @standardArray1;    
    my %badInFRSense0;
    for my $i (reverse(1..scalar (@array)-1)){ # FRFRFFFFRR -> FFFFFFRR
	if ($array[$i-1] eq $rightLetterOk && $array[$i] eq $leftLetterOk) {
	    $badInFRSense0{$i-1}=1;
	    $array[$i-1] = $leftLetterOk;

	}
    }
    my $numberRRF1 = scalar(@array)-scalar(keys %badInFRSense0);
    
    #for the second end:   
    @array = @standardArray2;
 
    my %badInFRSense3;
    for my $i (1..scalar(@array)-1){
	if ($array[$i-1] eq $rightLetterOk && $array[$i] eq $leftLetterOk) {
	    $badInFRSense3{$order2->[$i]}=1;
	    $array[$i] = $rightLetterOk;
	}
    }
    my $numberRRRFFF_or_FFF_2 = scalar(@array)-scalar(keys %badInFRSense3);

    @array = @standardArray2;
    my %badInFRSense5;
    for my $i (reverse(1..scalar (@array)-1)){ # FRFRFFFF -> FFFFFF
	if ($array[$i-1] eq $rightLetterOk && $array[$i] eq $leftLetterOk) {
	    $badInFRSense5{$i-1}=1;
	    $array[$i-1] = $leftLetterOk;
	}
    }
    my $numberRRF2 = scalar(@array)-scalar(keys %badInFRSense5);  
    
    if ($numberRRF1>=$numberRRRFFF_or_FFF_1 && $numberRRF1 >= $numberRRRFFF_or_FFF_2 && $numberRRF1 >=$numberRRF2) {
	return (1,%badInFRSense0);
    }   
    
    if ($numberRRRFFF_or_FFF_1 >=$numberRRF1 && $numberRRRFFF_or_FFF_1 >= $numberRRRFFF_or_FFF_2 && $numberRRRFFF_or_FFF_1 >= $numberRRF2) {
	return (1,%badInFRSense1);
    }
    
    if ($numberRRRFFF_or_FFF_2 >= $numberRRF1 && $numberRRRFFF_or_FFF_2 >= $numberRRRFFF_or_FFF_1 && $numberRRRFFF_or_FFF_2 >=$numberRRF2) {
	return (2,%badInFRSense3);
    }
    
    if ($numberRRF2 >=  $numberRRF1 && $numberRRF2 >= $numberRRRFFF_or_FFF_1 && $numberRRF2 >= $numberRRRFFF_or_FFF_2 ) {
	return (2,%badInFRSense5);
    }  
    
    #should not get here:
    print STDERR "Error in findBadInFRSenseSOLiDSolexa()!\n";
    return (1,%badInFRSense1); 
}

sub getOrderedStrands {
    my ($strand,$order) = ($_[0],$_[1]);
    my @arr;
    for my $i (0..scalar(@{$strand})-1) {
	push(@arr,$strand->[$order->[$i]-1]); 
    }    
    return @arr; 
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub checkClusters {
    
    my ($ifRenv,$coord_start_chr1_cluster1,$coord_start_chr1_cluster2,$coord_start_chr2_cluster1,$coord_start_chr2_cluster2) = @_;
    if ($ifRenv eq 'REVERSE_SENSE') {
	if ($coord_start_chr1_cluster1 <= $coord_start_chr1_cluster2) {
	    return ($coord_start_chr2_cluster1 <= $coord_start_chr2_cluster2)?1:0; 
	}
	return ($coord_start_chr2_cluster1 >= $coord_start_chr2_cluster2)?1:0; 
    }
    #if NORM
    if ($coord_start_chr1_cluster1 <= $coord_start_chr1_cluster2) {
        return ($coord_start_chr2_cluster1 >= $coord_start_chr2_cluster2)?1:0; 
    }
    return ($coord_start_chr2_cluster1 <= $coord_start_chr2_cluster2)?1:0; 
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub translateSolidToRF {
    my ($strandArr,$ends_orderArr)=@_;
    my @array;
    for my $i (0..scalar(@{$strandArr})-1) {					
	if ($ends_orderArr->[$i]==1 && $strandArr->[$i] eq 'F') {
	    push(@array,'F');
	}
	if ($ends_orderArr->[$i]==2 && $strandArr->[$i] eq 'F') {
	    push(@array,'R');
	}
	if ($ends_orderArr->[$i]==1 && $strandArr->[$i] eq 'R') {
	    push(@array,'R');   
	}
	if ($ends_orderArr->[$i]==2 && $strandArr->[$i] eq 'R') {
	    push(@array,'F');    
	}
    }
    return @array;
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#convert the links file to the circos format
sub links2segdup{
    
    my($id,$color_code,$links_file,$segdup_file)=@_;
    
    print LOG "# Converting to the circos format...\n";
    
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
    
    print LOG "# Converting to the bed format...\n";
    
    my $compare=1;
    if($links_file!~/compared$/){
	$compare=0;
	$tag_length->{none}->{1}=$tag_length->{1};
	$tag_length->{none}->{2}=$tag_length->{2};
    }
    
    #color-code hash table
    tie (my %hcolor,'Tie::IxHash');
    my %color_order;
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
sub getAllEntries{    
    my (@tab)=split (/,/,$_[0]);
    my @tab2;

    foreach my $v (@tab){
	
	$v=~s/[\(\)]//g;
	push(@tab2,$v);
    }
    return @tab2;
}#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub getAllEntriesWOspecialChar{    
    my (@tab)=split (/,/,$_[0]);
    my @tab2;

    foreach my $v (@tab){
	
	$v=~s/[\(\)\$\*\@]//g;
	push(@tab2,$v);
    }
    return @tab2;
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub links2SVfile{
    
    my($links_file,$sv_file)=@_;
    
    print LOG "# Converting to the sv output table...\n";
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
sub densityCalculation{
    
    my ($chr,$chrID,$file,$tag_length,$window_dist,$step,$mates_file,$mates_file_ref,$density_file,$input_format)=@_;
    
    my @sfile=split(/\./,$$mates_file[$file]);
    my $fchr=$sfile[$#sfile];
    
    my $fh = new FileHandle;
    
    my %density;
    my %density_ref;
    my @ratio;
    my ($cov,$cov_ref);
  
    #FREQUENCY CALCULATION PROCEDURE
    print LOG "# $fchr : Frequency calculation procedure...\n";
    &FreqCalculation(\%density,$chr,$chrID,$tag_length,$window_dist,$step,$$mates_file[$file],$input_format);
    &FreqCalculation(\%density_ref,$chr,$chrID,$tag_length,$window_dist,$step,$$mates_file_ref[$file],$input_format);
   
    #MAKING RATIO AND OUTPUT
    print LOG "\# Ratio calculation procedure...\n";
    $density_file=~s/\/mates\//\/density\//;
    $fh->open(">".$density_file) or die "$0: can't write in the output ".$density_file." :$!\n";
    
    foreach my $k (1..$chr->{nb_chrs}){
	foreach my $frag (1..$chr->{$k}->{nb_frag}){
		
	    @ratio= ($chr->{$k}->{name},
		    (${$chr->{$k}->{$frag}}[0]+1),
		    (${$chr->{$k}->{$frag}}[1]+1));
		
	    $cov=(exists $density{$k}{$frag}->{count})? $density{$k}{$frag}->{count}:0;
	    $cov_ref=(exists $density_ref{$k}{$frag}->{count})? $density_ref{$k}{$frag}->{count}:0;
		
	    push(@ratio,$cov,$cov_ref);
	    push(@ratio,log($cov/$cov_ref)) if($cov && $cov_ref);
	    push(@ratio,-log($cov_ref+1)) if(!$cov && $cov_ref);
	    push(@ratio,log($cov+1)) if($cov && !$cov_ref);
	    next if(!$cov && !$cov_ref);
		
	    print $fh join("\t",@ratio)."\n";
	}
    }
    
    $fh->close;
    print LOG "-- output created: $density_file\n";
    
    undef %density;
    undef %density_ref;
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub FreqCalculation{
    
    my ($density,$chr,$chrID,$tag_length,$window_dist,$step,$mates_file,$input_format) = @_;
    
    my @sfile=split(/\./,$mates_file);
    my $fchr=$sfile[$#sfile];
    my $fh = new FileHandle;
    
    my $nb_windows=0;
    my $warn=100000;
    my $record=0;
    my %pair;
 
    my ($sumX,$sumX2) = (0,0);
 
    print LOG "\# Frequency calculation for $mates_file...\n";
 
     if ($mates_file =~ /.gz$/) {
	$fh->open("gunzip -c $mates_file |") or die "$0: can't open ".$mates_file.":$!\n";
    }elsif($mates_file =~ /.bam$/){
	o$fh->open("$SAMTOOLS_BIN_DIR/samtools view $mates_file |") or die "$0: can't open ".$mates_file.":$!\n";#GALAXY
    }else{
	$fh->open("<".$mates_file) or die "$0: can't open ".$mates_file.":$!\n";
    }
    
    while(<$fh>){
	
	my @t=split;
	my $mate=$t[0];
	    
	my ($chr_read1, $chr_read2, $firstbase_read1, $firstbase_read2, $end_order_read1, $end_order_read2,);
	
	next if(exists $pair{$mate});
	
	next if (!readMateFile(\$chr_read1, \$chr_read2, \$firstbase_read1, \$firstbase_read2,\$end_order_read1, \$end_order_read2, \@t, $input_format,$tag_length));
	
	next unless (exists $chrID->{$chr_read1} || exists $chrID->{$chr_read2}); 
	($chr_read1, $chr_read2)= ($chrID->{$chr_read1},$chrID->{$chr_read2});
	
	$pair{$mate}=undef;
	$record++;
	
	my ($coord_start_read1,$coord_end_read1, $coord_start_read2,$coord_end_read2);
	
	recupCoords($firstbase_read1,\$coord_start_read1,\$coord_end_read1,$tag_length->{$end_order_read1},$input_format);
	recupCoords($firstbase_read2,\$coord_start_read2,\$coord_end_read2,$tag_length->{$end_order_read2},$input_format);
	
	my $length = abs($coord_start_read1-$coord_start_read2);
	$sumX += $length;							#add to sum and sum^2 for mean and variance calculation
	$sumX2 += $length*$length;	
	
	for(my $i=1;$i<=$chr->{$chr_read1}->{'nb_frag'};$i++){
	    
	    if (abs ($coord_start_read1-${$chr->{$chr_read1}->{$i}}[0]) <= $window_dist){
		
		&addToDensity($density,$chr_read1,$i,\$nb_windows)
		if(overlap($coord_start_read1,$coord_end_read2,${$chr->{$chr_read1}->{$i}}[0],${$chr->{$chr_read1}->{$i}}[1]));
		
	    }else{
		
		$i=getNextFrag($coord_start_read1,$i,${$chr->{$chr_read1}->{$i}}[0],$chr->{$chr_read1}->{nb_frag},$window_dist,$step);
	    }
	}
	
	if($record>=$warn){
	    print LOG "-- $warn mate-pairs analysed - $nb_windows points created\n";
	    $warn+=100000;
	}
    }
    $fh->close;

    print LOG "-- $fchr : Total : $record mate-pairs analysed - $nb_windows points created\n";

    if($record>0){
	
	my $mu = $sumX/$record;
	my $sigma = sqrt($sumX2/$record - $mu*$mu);
	print LOG "-- $fchr : mu length = $mu, sigma length = $sigma\n";
    }

}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub ratio2segdup{
    
    my($id,$density_file,$segdup_file)=@_;
    
    print LOG "# Converting to circos format...\n";
        
    open RATIO, "<$density_file" or die "$0: can't open $density_file :$!\n";
    open SEGDUP, ">$segdup_file" or die "$0: can't write in the output: $segdup_file :$!\n";
  
    while(<RATIO>){
	chomp;
	my ($chr1,$start1,$end1,$ratio)=(split /\t/)[0,1,2,5];	
	print SEGDUP "$id$chr1\t$start1\t$end1\t$ratio\n";
    }
    
    close RATIO;
    close SEGDUP;
    print LOG "-- output created: $segdup_file\n";
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub ratio2bedfile{
    
    my($density_file,$bed_file)=@_;
    
    print LOG "# Converting to bedGraph format...\n";
        
    open RATIO, "<$density_file" or die "$0: can't open $density_file :$!\n";
    open BED, ">$bed_file" or die "$0: can't write in the output: $bed_file :$!\n";
    print BED "track type=bedGraph name=\"$bed_file\" description=\"log ratios for cnv detection\" ".
	      "visibility=2 color=255,0,0 alwaysZero=\"On\"\n";
    
    while(<RATIO>){
	chomp;
	my ($chr1,$start1,$end1,$ratio)=(split /\t/)[0,1,2,5];
	$chr1 = "chr".$chr1 unless ($chr1 =~ m/chr/);
	print BED "$chr1\t".($start1-1)."\t$end1\t$ratio\n";
    }
    
    close RATIO;
    close BED;
    print LOG "-- output created: $bed_file\n";
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub inverseSense{
    
    my $mate_sense=$_[0];
    my %reverse=( 'F'  => 'R' , 'R'  => 'F' ,
		  'FF' => 'RR', 'RR' => 'FF',			
		  'FR' => 'RF', 'RF' => 'FR');
    return $reverse{$mate_sense};
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub getNextFrag{
    
    my ($read_start,$frag_num,$frag_start,$frag_last,$window_dist,$step)=@_;
    
    my $how_far = $read_start-$frag_start;
    my $nb_windows_toskip;
    
    if($how_far>0){
	
	$nb_windows_toskip=($how_far/$step)-($window_dist/$step);
	$nb_windows_toskip=~ s/\..*//;
	$nb_windows_toskip=0 if($nb_windows_toskip<0);
	return ($frag_num + $nb_windows_toskip);
    }
    return $frag_last;
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
sub recupCoords{
    
    my($c_hit,$cs_hit,$ce_hit,$tag_length,$input_format)=@_;
    my $strand = 'F';
    
    if ($c_hit=~s/^\-//) {
	$strand='R';
	    $$cs_hit=$c_hit;
	    $$ce_hit=$c_hit-($tag_length-1); 
    }else{
	$$cs_hit=$c_hit;
	$$ce_hit=$c_hit+($tag_length-1); 
    }
    return $strand;
    
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub overlap {
    my($cs_hit,$ce_hit,$cs_region,$ce_region)=@_;
    if( (($cs_hit < $cs_region) && ($ce_hit < $cs_region )) || (($cs_hit > $ce_region) && ($ce_hit > $ce_region )) ) {
        return 0;
    }
    return 1;
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub makeLink {
    
    my ($link,$chr1,$frag1,$chr2,$frag2,$mt,$nb)=@_;

    if($chr1>$chr2){
        ($chr1,$chr2)= ($chr2,$chr1);
        ($frag1,$frag2)= ($frag2,$frag1);
    }
    
     if($chr1 == $chr2){
        if($frag1>$frag2){
	    ($frag1,$frag2)= ($frag2,$frag1);
	}
    }
  
    if(!exists $link->{$chr1}->{$chr2}->{$frag1}->{$frag2}){
        $link->{$chr1}->{$chr2}->{$frag1}->{$frag2}=$mt;
	$$nb++;
    }elsif($link->{$chr1}->{$chr2}->{$frag1}->{$frag2}!~/(^|,)$mt(,|$)/){
        $link->{$chr1}->{$chr2}->{$frag1}->{$frag2}.=",$mt";
    }
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#fonction of adding the read to the density profile
sub addToDensity {
    
    my ($density,$chr1,$frag1,$nb)=@_;
    
    if(!exists $density->{$chr1}->{$frag1}->{count}){
            $density->{$chr1}->{$frag1}->{count}=1;
	    $$nb++;
    }else{
            $density->{$chr1}->{$frag1}->{count}++;
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
sub decimal{
    
  my $num=shift;
  my $digs_to_cut=shift;

  $num=sprintf("%.".($digs_to_cut-1)."f", $num) if ($num=~/\d+\.(\d){$digs_to_cut,}/);

  return $num;
}

#------------------------------------------------------------------------------#
sub max {
    
    my($max) = shift(@_);
    foreach my $temp (@_) {
        $max = $temp if $temp > $max;
    }
    return($max);
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub min {
    
    my($min) = shift(@_);
    foreach my $temp (@_) {
        $min = $temp if $temp < $min;
    }
    return($min);
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub sortTablebyIndex{
    my ($tab1,$tab2)=@_;
    my @tab3;
    
    foreach my $i (@$tab1){
	$tab3[$i]=$$tab2[$$tab1[$i]];
    }
    return @tab3;
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub round {
  my $number = shift || 0;
  my $dec = 10 ** (shift || 0);
  return int( $dec * $number + .5 * ($number <=> 0)) / $dec;
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub getUniqueTable{
    
    my (@tab)=@_;
    my (%saw,@out)=();
    undef %saw;
    return sort(grep(!$saw{$_}++, @tab)); 
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
sub catFiles {
    
    unlink("$_[1]") if(exists $_[1]);
    system qq( cat "$_" >> "$_[1]" ) for @{$_[0]};
}
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#check if the configuration file is correct
sub validateconfiguration{
    
    my %conf=%{$_[0]};
    my $list_prgs="@ARGV";
    
    my @general_params=qw(input_format mates_orientation read1_length read2_length mates_file cmap_file);
    my @detection_params=qw(split_mate_file window_size step_length split_mate_file);
    my @filtering_params=qw(split_link_file nb_pairs_threshold strand_filtering split_link_file);
    my @circos_params=qw(organism_id colorcode);
    my @bed_params=qw(colorcode);
    my @compare_params=qw(list_samples file_suffix);
    
    foreach my $dir ($conf{general}{output_dir},$conf{general}{tmp_dir}){
	
	unless (defined($dir)) {
	    $dir = ".";
	}
	unless (-d $dir){
	    mkdir $dir or die;
	}
	$dir.="/" if($dir!~/\/$/);
    }
    
    unless (defined($conf{general}{num_threads})) {
	    $conf{general}{num_threads} = 1;
	}
    $conf{general}{num_threads}=24 if($conf{general}{num_threads}>24);
    
    if($list_prgs!~/links2compare/){
    
	foreach my $p (@general_params){
	    die("Error Config : The parameter \"$p\" is not defined\n") if (!defined $conf{general}{$p});
	}
	
	$conf{general}{input_format}="sam" if($conf{general}{input_format} eq "bam");
	
	unless (defined($conf{general}{sv_type})) {
	    $conf{general}{sv_type} = "all";
	}
	
	$conf{general}{read_lengths}={ 1=> $conf{general}{read1_length}, 2=> $conf{general}{read2_length}};
    }
    
    if($list_prgs=~/(linking|cnv)/){
	
	foreach my $p (@detection_params){
	    die("Error Config : The parameter \"$p\" is not defined\n") if (!defined $conf{detection}{$p});
	}
	
	die("Error Config : The parameter \"mates_file_ref\" is not defined\n") if($list_prgs=~/cnv/ && !defined $conf{detection}{mates_file_ref});
	
	if($conf{detection}{step_length}>$conf{detection}{window_size}){
	    die("Error Config : Parameter \"step_length\" should not exceed \"window size\"\n");
	}
	
	unless (-d $conf{general}{tmp_dir}."/mates"){
	    mkdir $conf{general}{tmp_dir}."/mates" or die;
	}
	
	if($list_prgs=~/linking/){
	    unless (-d $conf{general}{tmp_dir}."/links"){
		mkdir $conf{general}{tmp_dir}."/links" or die;
	    }
	}
	if($list_prgs=~/cnv/){
	    unless (-d $conf{general}{tmp_dir}."/density"){
		mkdir $conf{general}{tmp_dir}."/density" or die;
	    }
	}
	
    }
    
    if($list_prgs=~/filtering/){
    
	foreach my $p (@filtering_params) {
	    die("Error Config : The filtering parameter \"$p\" is not defined\n") if (!defined $conf{filtering}{$p});
	    
	}
	
	if(defined($conf{filtering}{chromosomes})) {
	    my @chrs=split(",",$conf{filtering}{chromosomes});
	    my $exclude=($chrs[0]=~/^\-/)? 1:0;
	    for my $chrName (@chrs){
	      
		die("Error Config : The filtering parameter \"chromosomes\" is not valid\n")
		if(($chrName!~/^\-/ && $exclude) || ($chrName=~/^\-/ && !$exclude));
		
	    }
	}
	
	if (( $conf{filtering}{order_filtering} )&& !$conf{filtering}{strand_filtering}) {
	    die("Error Config : The parameter strand_filtering is set to \"0\" while order_filtering is selected".
		"\nChange strand_filtering to \"1\" if you want to use the order filtering\n");
	}
	if (( !defined($conf{filtering}{mu_length}) || !defined($conf{filtering}{sigma_length}) )&& $conf{filtering}{order_filtering}) {
	    die("Error Config : You should set parameters \"mu_length\" and \"sigma_length\" to use order filtering\n");
	}
	if (( $conf{filtering}{insert_size_filtering} )&& !$conf{filtering}{strand_filtering}) {
	    die("Error Config : The parameter strand_filtering is set to \"0\" while insert_size_filtering is selected".
		"\nChange strand_filtering to \"1\" if you want to use the insert size filtering\n");
	}
	if (( !defined($conf{filtering}{mu_length}) || !defined($conf{filtering}{sigma_length}) )&& $conf{filtering}{insert_size_filtering}) {
	    die("Error Config : You should set parameters \"mu_length\" and \"sigma_length\" to use discriminate insertions from deletions\n");
	}
	
	if (!defined($conf{filtering}{indel_sigma_threshold})) {
	    $conf{filtering}{indel_sigma_threshold} = 2;
	}
	if (!defined($conf{filtering}{dup_sigma_threshold})) {
	    $conf{filtering}{dup_sigma_threshold} = 2;
	}
	if (!defined($conf{filtering}{singleton_sigma_threshold})) {
	    $conf{filtering}{singleton_sigma_threshold} = 4;
	}
	
	if (!defined($conf{filtering}{nb_pairs_order_threshold})) {
	    $conf{filtering}{nb_pairs_order_threshold} = 1;
	}
	
	if (!defined($conf{filtering}{final_score_threshold})) {
	    $conf{filtering}{final_score_threshold} = 0.8;
	}
	
	if ($conf{filtering}{nb_pairs_order_threshold}>$conf{filtering}{nb_pairs_threshold}) {
	    die("Error Config : Parameter \"nb_pairs_order_threshold\" should not exceed \"nb_pairs_threshold\"\n");
	}
	
    }
    
    if($list_prgs=~/2circos$/){
	foreach my $p (@circos_params) {
	    next if($list_prgs=~/^ratio/ && $p eq "colorcode");
	    die("Error Config : The circos parameter \"$p\" is not defined\n") if (!defined $conf{circos}{$p});
	}
    }
    
    if($list_prgs=~/2bed$/){
	foreach my $p (@bed_params) {
	    die("Error Config : The bed parameter \"$p\" is not defined\n") if (!defined $conf{bed}{$p});
	}
    }
    
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
