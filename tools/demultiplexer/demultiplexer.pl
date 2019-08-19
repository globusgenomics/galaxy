#!/usr/bin/env perl


use strict;
use warnings;

use Carp;
use File::Basename;
use Getopt::Long;
use IO::File;
use IO::Pipe;
use IO::Zlib;
use Pod::Usage;
use List::MoreUtils qw(uniq);
use List::Util qw(max);
use Cwd;
use File::Spec;


my $start_time = (time);

my ($adapter_seq, $barcode_file, $no_skip_adp, $no_compress, $compression_level, 
    $tag_end, $tag_read, $mismatches, $file_format, $extra_tag_bases, $help, $allow_Ns,
    $output_id, $output_dir, $galaxy_dummy_output_file);

my @sequence_files = ( );

$tag_end         = 5;
$tag_read        = 1;
$extra_tag_bases = 0;
$mismatches      = 0;
$allow_Ns        = 0;

$adapter_seq = 'AGATCGGAAGAGCG';

my $trim_bases = 0;

my $trim_read;

my $rc = GetOptions(
                    "a|adapter-seq=s"              => \$adapter_seq,
                    "t|tag-file=s"                 => \$barcode_file,
                    "f|file-format=s"              => \$file_format,
                    "nd|no-skip-adapter-dimers"    => \$no_skip_adp,
                    "nc|no-compress-output"        => \$no_compress,
                    "c|compression-level=i"        => \$compression_level,
                    "e|tag-end|barcode-end=s"      => \$tag_end,
                    "r|tag-read|barcode-read=s"    => \$tag_read,
                    "m|mismatches=i"               => \$mismatches,
                    "i|sequence-file=s"            => \@sequence_files, 
                    "x|extra-tag-bases=i"          => \$extra_tag_bases,
                    "tb|trim-bases=i"              => \$trim_bases,
                    "tr|trim-read=s"               => \$trim_read,
                    "N|allow-Ns"                   => \$allow_Ns,
                    "oid|output-id=i"              => \$output_id,
                    "o|output-dir=s"               => \$output_dir,
                    "g|galaxy-dummy-output-file=s" => \$galaxy_dummy_output_file,
                    "h|help"                       => \$help,
                   );

pod2usage(1) if $help;
pod2usage(1) unless ($rc);
pod2usage(2) unless (defined($barcode_file));
pod2usage(2) unless (defined($file_format));
pod2usage(2) unless (($tag_end  == 3) || ($tag_end  == 5));
pod2usage(2) unless (($tag_read eq '1') || 
                     ($tag_read eq '2') || 
                     ($tag_read eq 'both'));
pod2usage(2) unless (@sequence_files >= 1);
pod2usage(2) unless (@sequence_files <= 2);
pod2usage(2) unless (defined($galaxy_dummy_output_file));

if (defined($trim_read)) {
    pod2usage(2) unless (($trim_read eq '1') ||
                         ($trim_read eq '2') ||
                         ($trim_read eq 'both'));
}

if ($trim_bases > 0) {
    unless (defined($trim_read)) {
        die 'need --trim-read with --trim-bases';
    }
}

if (defined($compression_level)) {

    if ($no_compress) {
        warn "ignoring --compression-level with --no-compress";
    }
    else {
    
        unless (
                ($compression_level >= 1) 
                 &&
                ($compression_level <= 9)
               ) {
            pod2usage(2);
        }
 
    }

}

my $seq_read_ref;

if ($file_format eq 'scarf') {
    $seq_read_ref = \&read_seq_scarf; 
}
elsif ($file_format eq 'fastq') {
    $seq_read_ref = \&read_seq_fastq;
}
else {
    die "Can't handle file format '$file_format'";
}

my $tag_trim_ref;

if ($tag_end == 5) {
    $tag_trim_ref = \&trim_tag_five_prime;
}
else {
    $tag_trim_ref = \&trim_tag_three_prime;
}

unless (defined($output_id)) {
    $output_id = 0;
}

unless (defined($output_dir)) {
    $output_dir = getcwd;
}

my $touch_output = `touch $galaxy_dummy_output_file`;

my %Bcodes; # hash to store barcodes

my @barcode_lengths_uniq = &get_barcodes; # sub routine to read barcodes
my $max_barcode_length = max(@barcode_lengths_uniq);

my %BC_counter; # hash to count number of sequences having each barcode;

my %OUT_file_i;   # hash to store named got
my %OUT_file_j;   # hash to store named got

my %OUT_handle_i; # 
my %OUT_handle_j; # 


my ($sequence_file_i, $sequence_file_j) = @sequence_files;

# get name of sequence file

my ($output_i) = fileparse($sequence_file_i);

my $output_j;

my $input_name = 'reads';


if (defined($sequence_file_j)) { 
    ($output_j) = fileparse($sequence_file_j);
    $input_name = 'read pairs'; 
}

# set up sequence counter and output files for each barcode
(my $output_i_cleaned = $output_i) =~ s/_/-/g;
my $output_j_cleaned;
if (defined($sequence_file_j)) {
    ($output_j_cleaned = $output_j) =~ s/_/-/g;
}

foreach my $tag (keys %Bcodes){

    $BC_counter{$tag} = 0; 
    my $fn_bn_i = join('_','primary',$output_id,$output_i_cleaned . '-' . $tag,'visible','fastq');
    my $fn_i = File::Spec->catfile($output_dir,$fn_bn_i);
    my $fh_i = open_file_write($fn_i);

    $OUT_file_i{$tag}   = $fn_i;
    $OUT_handle_i{$tag} = $fh_i;
    
    if (defined($sequence_file_j)) {
    
        my $fn_bn_j = join('_','primary',$output_id,$output_j_cleaned . '-' . $tag,'visible','fastq');
        my $fn_j = File::Spec->catfile($output_dir,$fn_bn_j);
 
        my $fh_j = open_file_write($fn_j);
 
        $OUT_file_j{$tag}   = $fn_j;
        $OUT_handle_j{$tag} = $fh_j; 

    }

}

# set up output/files for sequences without a barcode

my $output_unknown_bn_i = join('_','primary',$output_id,$output_i_cleaned . '-' . 'unknown' ,'visible','fastq');
my $output_unknown_i = File::Spec->catfile($output_dir,$output_unknown_bn_i);

my $out_u_i = open_file_write($output_unknown_i);

my $output_unknown_j;
my $out_u_j;

if (defined($output_j)) {

    my $output_unknown_bn_j = join('_','primary',$output_id,$output_j_cleaned . '-' . 'unknown' ,'visible','fastq');
    my $output_unknown_j = File::Spec->catfile($output_dir,$output_unknown_bn_j);

    $out_u_j = open_file_write($output_unknown_j); 
}

my $total_counter  = 0;
my $unknown        = 0;
my $adapter_dimers = 0;
my $fuzzy_matches  = 0;

my $in_i;
my $in_j;

$in_i = open_file_read($sequence_file_i);

if (defined($sequence_file_j)) {
    $in_j = open_file_read($sequence_file_j);
}

local $SIG{PIPE} = sub { die "caught SIGPIPE in the face" };
PAIR:
while(1) {

    my $skip_adapter_dimer = 0;

    my $seq_ref_i = $seq_read_ref->($in_i);
	my $seq_ref_j;

    unless (defined($seq_ref_i)) { last PAIR; }

	# Grab the sequence
    my ($first_seq_id_i, $seq_i_master, $qual_i_master) = @{$seq_ref_i}; # header, seq, quality
	my ($first_seq_id_j, $seq_j_master, $qual_j_master);
	
	if (defined($sequence_file_j)) {
		my $seq_ref_j = $seq_read_ref->($in_j);
		($first_seq_id_j, $seq_j_master, $qual_j_master) = @{$seq_ref_j};
	}

    $total_counter++;

	# Check if adapter sequence
    unless ($no_skip_adp) {
        if ($seq_i_master =~ /^\S{0,4}$adapter_seq/) {
            $skip_adapter_dimer = 1;
        }
    }

    my $temp_barcode_i;
	my $temp_barcode_j;
	my $fuzzy_match_found = 0; # Reset for each pair of reads
	my $match_found_i = 0;
	my $match_found_j = 0;
	my $match_bc_i;
	my $match_bc_j;
	my $match_seq_i;
	my $match_seq_j;
	my $match_qual_i;
	my $match_qual_j;
	my $seq_i = $seq_i_master;
	my $seq_j = $seq_j_master;
	my $qual_i = $qual_i_master;
	my $qual_j = $qual_j_master;
	
	foreach my $bc_length_uniq (@barcode_lengths_uniq) {
	
		# Analyze barcode if present on read 1
		if (($tag_read eq '1') || ($tag_read eq 'both')) {

			($temp_barcode_i, $seq_i, $qual_i) = 
				@{$tag_trim_ref->($bc_length_uniq, $extra_tag_bases, $seq_i_master, $qual_i_master)}; # Must keep _master values the same since trimming in a for loop
				
			# If no mismatches AND count Ns as mismatches, and valid barcode, increase match counter and store barcode, seq, and qual
			if ($mismatches == 0 && $allow_Ns == 0) {
			
				if (exists $Bcodes{$temp_barcode_i}) {
						
						$match_found_i++;
						$match_bc_i = $temp_barcode_i;
						$match_seq_i = $seq_i;
						$match_qual_i = $qual_i;				
				}						
			
			# If mismatches
			} else {
			
				if (exists $Bcodes{$temp_barcode_i}) {
						
						$match_found_i++;
						$match_bc_i = $temp_barcode_i;
						$match_seq_i = $seq_i;
						$match_qual_i = $qual_i;				
				
				}	else { # If not a valid barcode
				
					# print "TRYING FUZZY...\n";
					
					my $match = fuzzy_match($temp_barcode_i, \%Bcodes, $mismatches, $allow_Ns);
					
					# If only 1 fuzzy match (fuzzy_match returns undef if > 1 matches), increase match counter and store barcode, seq, and qual
					
					# print "THIS WAS THE MATCH \"$match\"\n";
					if (defined($match)) { 
						
						$match_bc_i = $match; # This is the fuzzy matched valid barcode
						$match_seq_i = $seq_i;
						$match_qual_i = $qual_i;
						$match_found_i++;
						
						# Update fuzzy match counter if applicable
						if ($fuzzy_match_found == 0) { # Increase the counter only once for a pair of reads (might have multiple barcodes of possibly different lengths match)
							$fuzzy_matches++;
							$fuzzy_match_found++;
						}
					}   
				} 
			}
		}


	 
		if (defined($sequence_file_j)) {

			# Check if adapter sequence
			unless (defined($no_skip_adp)) {
				if ($seq_j_master =~ /^\S{0,4}$adapter_seq/) {
					$skip_adapter_dimer = 1;
				}
			}
			
			# Analyze barcode if present on read 2
			if (($tag_read eq '2') || ($tag_read eq 'both')) {
				
				($temp_barcode_j, $seq_j, $qual_j) = 
					@{$tag_trim_ref->($bc_length_uniq, $extra_tag_bases, $seq_j_master, $qual_j_master)};
				
				# If no mismatches, and valid barcode, increase match counter and store barcode, seq, and qual
				if ($mismatches == 0 && $allow_Ns == 0) {
				    if (exists $Bcodes{$temp_barcode_j}) {
					
						$match_found_j++;
						$match_bc_j = $temp_barcode_j;
						$match_seq_j = $seq_j;
						$match_qual_j = $qual_j;
					}
				
				# If mismatches
				} else {
				
					if (exists $Bcodes{$temp_barcode_j}) {
					
						$match_found_j++;
						$match_bc_j = $temp_barcode_j;
						$match_seq_j = $seq_j;
						$match_qual_j = $qual_j;
						
					} else { # If not a valid barcode
					
						my $match = fuzzy_match($temp_barcode_j, \%Bcodes, $mismatches, $allow_Ns);
						
						# If only 1 fuzzy match (fuzzy_match returns undef if > 1 matches), increase match counter and store barcode, seq, and qual
						if (defined($match)) {
						
							$match_bc_j = $match;
							$match_seq_j = $seq_j;
							$match_qual_j = $qual_j;
							$match_found_j++;
							
							# Update fuzzy match counter if applicable
							if ($fuzzy_match_found == 0) {
								$fuzzy_matches++; # Am only updating this 1X now...okay?
								$fuzzy_match_found++;
							}
						}
					}
				}
			}
		}
	}

    unless ($no_skip_adp) {
   
        ## move on to the next read(s) only after both input files have been read from 
        if ($skip_adapter_dimer) {
           $adapter_dimers++;
            next PAIR;
        }

    }
	
	# If no matches found, update seq and qual for unknown
	if (( defined $match_seq_i) && (defined $match_seq_j)) {
	
		# Both seq i and j should be set
		
	} elsif (defined $match_seq_i) { # Need to set j (if applicable)
		
		# Tag found on R1.  Update R2 fields (if R2 exists).
		if (defined $sequence_file_j) {
			
			$match_seq_j = $seq_j;
			$match_qual_j = $qual_j;
		}
		
		# Tag found on R2.  Update R1 fields.
	} elsif (defined $match_seq_j) { # Need to set i
	
		$match_seq_i = $seq_i;
		$match_qual_i = $qual_i;
	
	} else { # Tag was not found: need to set R1 fields and R2 fields (if R2 exists)

		$match_seq_i = $seq_i;
		$match_qual_i = $qual_i;
	
		if (defined $sequence_file_j) {
		
			$match_seq_j = $seq_j;
			$match_qual_j = $qual_j;
		}
	}
	
	
	# Trim off (more) 3' bases if $trim_bases (e.g., trim off the last 6 bases due to poor read quality)
	if ($trim_bases) {
	
		if (($trim_read eq '1') || ($trim_read eq 'both')) {
			$match_seq_i  = substr($match_seq_i, 0, 0 - $trim_bases);
			$match_qual_i = substr($match_qual_i, 0, 0 - $trim_bases);
		}
		
		if (defined $sequence_file_j) {
		
			if (($trim_read eq '2') || ($trim_read eq 'both')) {
				$match_seq_j  = substr($match_seq_j, 0, 0 - $trim_bases);
				$match_qual_j = substr($match_qual_j, 0, 0 - $trim_bases);
			}
		}
	}
	
			
    ## sanity check to make sure we don't drift out of sync (and detect if the input files aren't ordered properly)
    if (defined($sequence_file_j)) {
        if (substr($first_seq_id_i, 0, length($first_seq_id_i) - 1) ne substr($first_seq_id_j, 0, length($first_seq_id_j) - 1)) {
            die "OUT OF SYNC AFTER $total_counter READ PAIRS - first read is $first_seq_id_i while second read is $first_seq_id_j";
        }
    }

    unless (length($match_seq_i) == length($match_qual_i)) {
        die join(
                 "\t",
                 'read with seq length not equal to qual length',
                 $first_seq_id_i,
                 $match_seq_i,
                 $match_qual_i,
                );
    }
        
    if (defined($sequence_file_j)) {
        unless (($trim_bases > 0) && ($trim_read ne 'both')) { # Can only compare R1 and R2 lengths if the data has not been trimmed with the trim bases option (e.g., if bad quality bases).
            if (length($match_seq_i) != length($match_seq_j)) {
                die printf("read pair with different seq lengths %d (%s) and %d (%s)\n",length($match_seq_i),$first_seq_id_i,length($match_seq_j),$first_seq_id_j);
            }
            if (length($match_qual_i) != length($match_qual_j)) {
                die join(
                         "\t",
                         'read pair with different qual lengths after tag trimming',
                         $first_seq_id_i,
                         $first_seq_id_j,
                        );
            }
        }
        unless (length($match_seq_j) == length($match_qual_j)) {
            die join(
                     "\t",
                     'read with seq length not equal to qual length',
                     $first_seq_id_j,
                     $match_seq_j,
                     $match_qual_j,
                    );
        }
    }
	
	# Check that only 1 valid match was made
	if ($tag_read == 1) {
			
		unless ($match_found_i == 1) {
			
			# Add to unknown
			$unknown++;

			write_seq_fastq($out_u_i, $first_seq_id_i, $match_seq_i, $match_qual_i);

			if (defined($sequence_file_j)) {
			
				write_seq_fastq($out_u_j, $first_seq_id_j, $match_seq_j, $match_qual_j);
		   }
		   
		   next PAIR;
		   
		}
		
	} elsif (defined $sequence_file_j) {
	
		if ($tag_read == 2) {
	
			unless ($match_found_j == 1) {
				
				# Add to unknown
				$unknown++;

				write_seq_fastq($out_u_i, $first_seq_id_i, $match_seq_i, $match_qual_i);			
				write_seq_fastq($out_u_j, $first_seq_id_j, $match_seq_j, $match_qual_j);
				
				next PAIR;

			}
		
		} elsif ($tag_read eq 'both') {
	
			unless (($match_found_i == 1) || ($match_found_j == 1)) {
			
				# Add to unknown
				$unknown++;

				write_seq_fastq($out_u_i, $first_seq_id_i, $match_seq_i, $match_qual_i);			
				write_seq_fastq($out_u_j, $first_seq_id_j, $match_seq_j, $match_qual_j);
				
				next PAIR;
			}

		} else {
		
			croak "Unrecognized tag read \"$tag_read\"!  Valid values are 1, 2, both.\n\n";
		}
	} else {
	
		croak "Should not be possible: tag read is not 1 and no R2 file defined!\n\n";
	}

	# If the barocdes on R1 and R2 don't match, add to unknown
    if (defined($sequence_file_j) && $tag_read eq 'both') {
    
		if (
            (exists($Bcodes{$match_bc_i})) &&
            (exists($Bcodes{$match_bc_j})) &&
            ($match_bc_i ne $match_bc_j) 
           ) {

            warn join("\t", 
                      'found a pair with different valid barcodes',
                      $first_seq_id_i,
                      $first_seq_id_j,
                      $match_bc_i,
                      $match_bc_j);
            
            $unknown++;

            write_seq_fastq($out_u_i, $first_seq_id_i, $match_seq_i, $match_qual_i);
            write_seq_fastq($out_u_j, $first_seq_id_j, $match_seq_j, $match_qual_j);
            
            next PAIR;

        }

    } 
	
	if ((($tag_read == 1) || ($tag_read eq 'both')) &&
        (exists $Bcodes{$match_bc_i})) {

        $BC_counter{$match_bc_i}++;

        write_seq_fastq($OUT_handle_i{$match_bc_i}, $first_seq_id_i, $match_seq_i, $match_qual_i);

        if (defined($sequence_file_j)) {        
            write_seq_fastq($OUT_handle_j{$match_bc_i}, $first_seq_id_j, $match_seq_j, $match_qual_j);        
        }
        
        next PAIR;
	}

    if (defined($sequence_file_j)) {

        if (
            (($tag_read == 2) || ($tag_read eq 'both')) &&
            (exists $Bcodes{$match_bc_j})
           ) {

            $BC_counter{$match_bc_j}++;

            write_seq_fastq($OUT_handle_i{$match_bc_j}, $first_seq_id_i, $match_seq_i, $match_qual_i); 
            write_seq_fastq($OUT_handle_j{$match_bc_j}, $first_seq_id_j, $match_seq_j, $match_qual_j); 
            
            next PAIR;

        }
	}
   
      
   $unknown++;

   write_seq_fastq($out_u_i, $first_seq_id_i, $match_seq_i, $match_qual_i);

   if (defined($sequence_file_j)) {
	   write_seq_fastq($out_u_j, $first_seq_id_j, $match_seq_j, $match_qual_j);
   }
   

}

print "There were $total_counter Total $input_name\n";



foreach my $x (keys %Bcodes){   

    $OUT_handle_i{$x}->close();

    if (exists($OUT_handle_j{$x})) {
        $OUT_handle_j{$x}->close();
    }

    print "There were $BC_counter{$x} Barcode1 ($x) $input_name\n";

}

print "There were $unknown unknown $input_name\n";

unless ($no_skip_adp) {
    print "There were $adapter_dimers adapter dimers\n";
    print "The 3' adapter sequence used was $adapter_seq\n";
}

print "There were $fuzzy_matches fuzzy sequence tag matches\n";




my $stop_time = (time);

my $elapsed_time = $stop_time - $start_time;

print "\nThis program took $elapsed_time seconds to finish\n\n";



sub get_barcodes {


    my $fh = IO::File->new($barcode_file, 'r') || 
        die "\n\n Could not open barcode file: $barcode_file\n\n";

    my $barcode_length = 100;
    my $barcode_count  = 0;
	my @barcode_lengths;

    while (my $line = <$fh>){

        chomp $line;
       
        my @cols = split /\t/, $line;

        my $barcode = $cols[0];

        if ($barcode =~ /([ACTGactg]+)/){

            $barcode          = uc $barcode;
            $barcode_length   = length $barcode;
            $Bcodes{$barcode} = $barcode_length;
			push(@barcode_lengths,$barcode_length);

            $barcode_count++;

        }
        elsif ($barcode =~ /\S+/){
            print "\n\nERROR: The barcode file \"$barcode_file\" contains non canonical nucleotides (\"$barcode\") at line $.\n\n";
            exit;
        }
		else {
		
			die "\n\nERROR: The barcode file \"$barcode_file\" does not contain a valid barcode (\"$barcode\") at line $.\n\n";
		}

    }
	
	my @barcode_lengths_uniq = uniq(@barcode_lengths);

    my $warn = "no";

    if (%Bcodes){ # Print summary of barcodes

        my $version = 'V1.31';
        
        print "======================", "\n";
        print "demultiplexer.pl $version", "\n";
        print "======================", "\n";
        
        print "\n There are $barcode_count barcodes in file: \"$barcode_file\"\n";

		print "------------------------------------------\n";

		my $length_first = $barcode_lengths[0];
		
        foreach my $code (keys %Bcodes){

            print "$code\n";
			
			# Check if barcodes have the same length.  If not, warn user.  (Either there is a typo in $barcode_file or the data has different length indexes, e.g., GTAC indexes and Illumina indexes used together.)

			if ($Bcodes{$code} != $length_first) { 
			
                $warn     = "yes";

            }

        }

		print "==============================\n";

    }
    else { # No barcodes were found in the file

        print "\nThere was a problem reading the barcode file: $barcode_file\n";
        print "Please input one barcode per line\n\n";
        exit;

    }

    if ($warn eq "yes"){ # barcoces are not the same length

        print "\n\nWarning: Indexes are not the same length.  Their lengths were " . join (", ",sort (@barcode_lengths)) . "\n\n";

    }

    return @barcode_lengths_uniq;

}

sub read_seq_fastq {

    my ($fh) = @_; 


    my @data = ( );

    $#data = 2;

    for (0..2) {

        my $line = <$fh>;
        
        unless (defined($line)) { return undef; }
        
        chomp $line;
        
        ## Foul hack to skip second seq_id
        if ($_ == 2) { 
            $line = <$fh>;
            chomp $line; 
        }
        elsif ($_ == 0) {
            unless (substr($line, 0, 1) eq '@') {
                confess("this does not look like a valid fastq line: " .  $line);
            }       
            $line = substr($line, 1); 
        }
  
        $data[$_] = $line;

    }

    return \@data;
    
}

sub write_seq_fastq {

    my ($fh, $seq_id, $seq, $qual) = @_;


    my $argstring = join ' ', @_;

    unless (defined($fh))     { confess 'undefined $fh passed in';     }
    unless (defined($seq_id)) { confess 'undefined $seq_id passed in'; }
    unless (defined($seq))    { confess 'undefined $seq passed in';    }
    unless (defined($qual))   { confess 'undefined $qual passed in';   }
    
    print $fh '@', $seq_id, "\n";
    print $fh $seq, "\n";
    print $fh '+', "\n";
    print $fh $qual, "\n";

}

sub read_seq_scarf {

    my ($fh) = @_;


    unless (defined($fh)) {
        confess 'undefined $fh passed in'; 
    }
    
    my $line = <$fh>;

    unless (defined($line)) { return undef; }

    chomp $line;

    my @cols = split(':', $line);

    unless (@cols == 7) {
        confess("invalid scarf line (wrong column count) at line $.: $line");
    }

    foreach my $col (@cols) {
        unless (defined($col)) {
            confess "invalid scarf line (undefined field) at line $.: $line";
        }
    }

    unless (length($cols[5]) == length($cols[6])) {
        confess "invalid scarf line (different length seq and qual at line $.: $line";
    }

    my @data = (join(':', @cols[0..4]), @cols[5,6]);
    
    return \@data;

}

sub trim_tag_five_prime {

    my ($tag_length, $extra_tag_bases, $seq, $qual) = @_;

    
    unless (defined($tag_length))      { confess 'undefined $tag_length passed in';      }
    unless (defined($extra_tag_bases)) { confess 'undefined $extra_tag_bases passed in'; }
    unless (defined($seq))             { confess 'undefined $seq passed in';             }
    unless (defined($qual))            { confess 'undefined $qual passed in';            }

    my $tag;

    $tag  = substr($seq,  0, $tag_length);

    $seq  = substr($seq,  $max_barcode_length + $extra_tag_bases); 
    $qual = substr($qual, $max_barcode_length + $extra_tag_bases); 

    return [ $tag, $seq, $qual ];

}

sub trim_tag_three_prime {

    my ($tag_length, $extra_tag_bases, $seq, $qual) = @_;

    unless (defined($tag_length))      { confess 'undefined $tag_length passed in';      }
    unless (defined($extra_tag_bases)) { confess 'undefined $extra_tag_bases passed in'; }
    unless (defined($seq))             { confess 'undefined $seq passed in';             }
    unless (defined($qual))            { confess 'undefined $qual passed in';            }

    my $tag;
	
	# Assume that $max_barcode_length + $extra_tag_bases is the indexed read
	# RRRRIIIIII (R = read base, I = index base, A = adapter base)
	# RRRRIIIIIA
	# Therefore, must give special treatment for barcodes < max barcode length
	
    $tag  = substr($seq,  0-($max_barcode_length + $extra_tag_bases), $tag_length);

    $seq  = substr($seq,  0, 0-($max_barcode_length + $extra_tag_bases));
    $qual = substr($qual, 0, 0-($max_barcode_length + $extra_tag_bases));

    return [ $tag, $seq, $qual ];

}

sub open_file_read {

    my ($fn) = @_;


    unless (defined($fn)) { confess 'undefined $fn passed in'; }

    my $fh;

    if ($fn =~ /\.gz$/) {

        $fh = IO::Zlib->new();

        $fh->open($fn, 'rb');

        unless (defined($fh)) {
            confess "Couldn't open '$fn' for reading";
        }

    }

    else {

        $fh = IO::File->new();

        $fh->open($fn, 'r');

        unless (defined($fh)) {
            confess "Couldn't open '$fn' for reading";
        }

    }

    return $fh;

}

sub open_file_write {

    my ($fn) = @_;


    unless (defined($fn)) { confess 'undefined $fn passed in'; }

    my $fh;
 
    if ($fn =~ /\.gz$/) {

        $fh = IO::Pipe->new();

        my $cmd = 'gzip -c';

        if (defined($compression_level)) {
            $cmd .= " -$compression_level";
        }

        $fh->writer("$cmd > $fn");

        unless (defined($fh)) {
            confess "Couldn't open '$fn' for writing";
        }
 
    }

    else {

        $fh = IO::File->new();

        $fh->open($fn, 'w');

        unless (defined($fh)) {
            confess "Couldn't open '$fn' for writing";
        }

    }

    return $fh;

}

#Fuzzy or allow Ns to be matches
sub fuzzy_match {

    my ($sequence_tag, $tag_ref, $allowed_mismatches, $allow_Ns) = (@_);


    unless (defined($sequence_tag))       { confess 'undefined $sequence_tag passed in';         }
    unless (defined($tag_ref))            { confess 'undefined $tag_ref passed in';              }
    unless (defined($allowed_mismatches)) { confess 'undefined $allowed_mismatches passed in';   }
    unless (defined($allow_Ns)) { confess 'undefined $allow_Ns passed in';   }

    my @matches = ( );

    foreach my $tag (keys %{$tag_ref}) {
	
        my $mismatches  = 0;

        my @tag          = split '', $tag;
        my @sequence_tag = split '', $sequence_tag;
		
		 if ($#tag == $#sequence_tag) { # Only want to compare a 6-base sequence_tag with a 6-base barcode


			foreach my $pos (0..$#tag) {
			    my $base = $sequence_tag[$pos];
			    unless (($base eq $tag[$pos]) || ($allow_Ns && (uc($base) eq 'N' || ($base eq '.')))) { 
				$mismatches++; 
			    } 
			}

			if ($mismatches <= $allowed_mismatches) { 
			    push @matches, $tag; 
			}
		 }
	
    }

    if (@matches == 1) {
        return $matches[0];
    }
    else {
        return undef;
    }  

}

__END__

=head1 NAME

demultiplexer.pl - Demultiplex barcoded/indexed scarf/fastq sequence files to multiple fastq files.  Script can handle fuzzy matching as well as tags of different lengths.

=head1 SYNOPSIS

demultiplexer.pl --tag-file barcodes.txt --file-format scarf --tag-end 5 --tag-read both --sequence-file s_1_sequence.txt

demultiplexer.pl -t barcodes.txt -f scarf --e 5 -r both -i s_1_sequence.txt

demultiplexer.pl --tag-file indexes.txt --file-format fastq --tag-end 3 --tag-read 1 --sequence-file s_1_1_sequence.txt --sequence-file s_1_2_sequence.txt

demultiplexer.pl -t indexes.txt -f fastq -e 3 -r 1 -i s_1_1_sequence.txt -i s_1_2_sequence.txt

=head1 OPTIONS

=over 8

=item B<-h/--help>

Prints a brief help message and exits.

=item B<-f/--file-format>

Specifies the format of the input sequence file(s).  Current supported formats are scarf and fastq.

=item B<-t/--tag-file>

Specifies the file containing the sequence tags (barcodes/indexes) to use to bin the sequences from the input file(s).  Assumed to be tab delimited (to support chip_seq), with the tag in the 1st (or only) column.  Columns other than
the 1st/only are presently ignored.

The number of mismatches allowed is controlled by the -mismatches option.  

=item B<-nd/--no-skip-adapter-dimers>

Disables skipping of sequences that appear to be adapter dimers (default behaviour is to drop them on the floor).

=item B<-nc/--no-compress-output>

Disables gzip compression of the output files.

=item B<-c/--compression-level>

Compression level to pass to gzip (legal values are 1-9 for minimal to maximal
respectively).  See the gzip documentation for details.  If this is not specified, no argument is passed to gzip, resulting in the gzip default compression level.

=item B<-e/--barcode-end/--tag-end>

Specifies which end of the input sequences/reads the sequence tag (barcode/index) is at.  Legal values are 5 (5') or 3 (3').

=item B<-r/--tag-read>

Specifies which read has the tag (barcode/index).  Legal values are 1, 2 or both.  Only meaningful in paired end mode (two -i/-sequence-file args).

=item B<-m/--mismatches>

Perform fuzzy lookup on the sequence tags.  If the barcode/index for a read isn't an exact match for a sequence in the file specified by -tag-file, perform a fuzzy
match allowing up to the specified number of mismatches.  Note that the match must be unique.  That is, after allowing for N mismatches, there must be only one possible candidate barcode/index.

=item B<-N/--allow-Ns>

By default, Ns are treated as mismatches.  By specifying the -N flag, Ns in the index/barcode will not be counted towards the mismatch total.

=item B<-i/--sequence-file>

Specifies the sequence file(s) to use as input.  Use this option once for data from fragment lanes, twice for paired end lanes.

=item B<-x/--extra-tag-bases>

Specifies the number of extra bases following the sequence tag (barcode/index) that are not part of the tag, but shouldn't end up in the output sequences, either.  For example, the T overhang on 5' barcodes.  

=item B<-a/--adapter-seq>

Specifies the 3' adapter sequence used when filtering adapter-dimers.  Default is 'AGATCGGAAGAGCG'.

=item B<-tb/--trim-bases>

Trim N bases from the 3' end of the read(s) specified by -trim-read.  This is useful when you need to trim off cycles due to bad data quality.  Note: this trimming is performed *after* barcode trimming (if applicable).

=item B<-tr/--trim-read>

Legal values are 1, 2 and both.  Trims the number of bases specified by -trim bases from read 1, read 2, or both read 1 and read 2.

=item B<-g/galaxy-dummy-output-file>

Galaxy dummy output file.  Used for multiple output file workaround (see http://wiki.g2.bx.psu.edu/Admin/Tools/Multiple%20Output%20Files).  A blank file with this filename is created in the current working directory.  Required.

=item B<-oid/--output-id>

Output ID of galaxy dummy output file.  Determined by Galaxy.  The ID is used to name the output files so Galaxy knows where to find them.  Default is 0.

=item B<-o/--output-dir>

Output directory.  If run in galaxy, set this to $__new_file_path__.  Default is the current working directory.

=back

=head1 DESCRIPTION

B<demultiplex.pl> will read the input file(s) and bin sequences/reads by the sequence tag (barcode/index) present.  One output file will be created for each tag, plus one for sequence with an unknown/unrecognized match.  

=cut

