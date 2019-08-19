#!/usr/bin/perl

use strict;
use warnings;

use lib '/home/comp/gtac/gtac/lib/perl5';

use Carp;
use Email::Valid;
use File::Basename;
use File::Temp;
use Getopt::Long;
use HTML::TreeBuilder;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Pod::Usage;
use WWW::Mechanize;
use Cwd 'abs_path';


my ($dbsnp_130, $email_address, $help, 
    $input_file, $output_file, $protein_sequence,
    $gene_locations);

GetOptions(
           "db|dbsnp-130"       => \$dbsnp_130,
           "e|email-address=s", => \$email_address,
           "i|input-file=s"     => \$input_file,
           "o|output-file=s"    => \$output_file,
           "p|protein-sequence" => \$protein_sequence,
           "g|gene-locations=s" => \$gene_locations,
           "h|help" => \$help,
          ) or pod2usage(2);

pod2usage(1) if $help;
pod2usage(2) unless (defined($email_address));         
pod2usage(2) unless (defined($input_file));         
pod2usage(2) unless (defined($output_file));         


unless (-e $input_file) {
    warn "input file '$input_file' does not exist";
    exit(1);
}

unless (-f $input_file) {
    #warn "input file '$input_file' is not a file";
    print "input file '$input_file' is not a file";
}

if (-e $output_file) {
    #warn "output file '$output_file' exists - will overwite";
    print "output file '$output_file' exists - will overwite";
}

# Seattle seq does not like it if you submit an empty file.
unless (-s $input_file) {
    carp "Warning: input file '$input_file' is empty.  This might be because there are no SNPs in the targeted region or there was an error.  Input file will not be submitted to seattle seq, instead, the output file, '$output_file' file will be touched.\n\n";
    
    if (system("touch $output_file") == 0) {
        exit(0);
    } else {
        carp "Error: failed to touch $output_file\n\n";
    }
}

unless (Email::Valid->address($email_address)) {
    warn "invalid email address '$email_address'";
    exit(1);
}

$output_file = abs_path($output_file);
my ($input_basename, $input_path) = fileparse($input_file, qr/\.[^.]*/);
my ($input_basename_w_suffix) = fileparse($input_file);
chdir $input_path or croak "ERROR: could not cd to '$input_path': $!\n\n";


my $start_url = 'http://gvs.gs.washington.edu/SeattleSeqAnnotation'; 

unless (defined($dbsnp_130)) {
    $start_url = 'http://snp.gs.washington.edu/SeattleSeqAnnotation131/';
}

if (defined($gene_locations)) {

    unless (
            ($gene_locations eq 'NCBI')
            ||
            ($gene_locations eq 'CCDS')
            ||
            ($gene_locations eq 'both')
           ) {
        warn "invalid gene-location, should be NCBI, CCDS or both";
        exit(1);
    }

    if (defined($dbsnp_130)) {
        if ($gene_locations eq 'CCDS') {
            $gene_locations = 'CCDS2008'
        }
        elsif ($gene_locations eq 'both') {
            $gene_locations = 'both2008';
        }
    }
    else {
        if ($gene_locations eq 'CCDS') {
            $gene_locations = 'CCDS2010';
        }
        elsif ($gene_locations eq 'both') {
            $gene_locations = 'both2010';
        }
    }

}

unless (defined($gene_locations)) {
    if (defined($dbsnp_130)) {
        $gene_locations = '';
    }
    else {
        $gene_locations = 'both2010';
    }
}

my $busy_count = 0;
my $max_busy   = 10;

my $agent = WWW::Mechanize->new('autocheck' => 1);
$agent->agent_alias('Mac Safari');

 SUBMIT:
    while(1){
	
	$agent->get($start_url);

	$agent->form_name('GenotypeSummary');

	$agent->tick('columns', 'clinicalAssociation');
	$agent->field('EMail', $email_address);
	$agent->field('allele1Column', '4');
	$agent->field('gFetch', 'Display Genotypes');
	$agent->tick('columns', 'geneList');
	$agent->tick('columns', 'chimpAllele');
	$agent->field('allele2Column', '4');
	$agent->tick('columns', 'dbSNPValidation');
	$agent->field('fileFormat', 'custom');
	$agent->tick('columns', 'allelesDBSNP');
	$agent->field('chromosomeColumn', '1');
	$agent->field('GenotypeFile', $input_basename_w_suffix); # File must be in cwd
	$agent->field('referenceColumn', '3');
	$agent->tick('columns', 'scorePhastCons');
	$agent->tick('columns', 'consScoreGERP');
	$agent->field('locationColumn', '2');
	$agent->tick('columns', 'repeats');
	$agent->tick('columns', 'hasGenotypes');
	$agent->tick('columns', 'CNV');
	$agent->tick('columns', 'HapMapFreq');
	$agent->tick('columns', 'polyPhen');
	$agent->field('outputFileFormat', 'original');
	$agent->field('HapMapFreqType', 'HapMapFreqMinor');
	$agent->field('gFetch', 'Display Genotypes');
	
	if (defined($dbsnp_130)) {
	    $agent->tick('columns', 'allelesMaq');
	}
	else {
	    $agent->tick('columns', 'sampleAlleles');
	}
	
	if (defined($protein_sequence)) {
	    $agent->tick('columns', 'proteinSequence');
	}
	else {
	    $agent->untick('columns', 'proteinSequence');
	}

	$agent->field('geneData', $gene_locations);

	print "submitting...";
	
	$agent->click();
	
	if (server_busy($agent->content())) {
	    $busy_count++;

	    if ($busy_count <= $max_busy) {
		print "server is busy...sleeping", "\n";
		sleep(600);
		next SUBMIT;
	    }
	    else {
		print "server is busy...giving up";
		exit(1);
		last SUBMIT;
	    }
	    
	}
	
	print  "done!", "\n";

	last SUBMIT;

}

$agent->follow_link(text_regex => qr/monitor job progress/i);

my $done = 0;

while (!$done) {

    print "sleeping...";
    sleep(60);
    
    $agent->reload();

    my @links = $agent->links();
 
    @links = grep { defined($_->text()) } @links;

    my ($result_link) = grep { $_->text() =~ qr/show table/ } @links;

    if (defined($result_link)) {
 
        $agent->follow_link('text' => $result_link->text());

        if (!defined($agent->content()) || ($agent->content eq '')) {
            die "got blank page after following show table link";
        }

        my $download_url = download_url($input_basename, $agent->content(), $start_url);

        unless (defined($download_url)) {
            die "failed to parse download_url";
        }

        print "fetching result..."; 

        my $temp_fh = File::Temp->new();
        $temp_fh->close();
        my $temp_fn = $temp_fh->filename();

        $agent->get($download_url, ':content_file' => $temp_fn);

        gunzip $temp_fn => $output_file or die "gunzip failed: $GunzipError\n";        

        print "done!", "\n";

        $done = 1; 

    }
    else {

        my $progress_text = check_progress($agent->content());

        if (defined($progress_text)) {
            print $progress_text, "\n";
        }

    } 

}

sub check_progress {

    my ($content) = @_;


    my $root = HTML::TreeBuilder->new();

    $root->parse($content);
    $root->eof();

    my ($progress_tag) = $root->look_down(
                                          '_tag' => 'big',
                                         );
   
    my $progress_text;

    if (defined($progress_tag)) {

        my ($content_ref) = $progress_tag->content();

        if (defined($content_ref)) {    
            $progress_text = $content_ref->[0];
        }

    }

    $root->delete();

    return $progress_text;

}

sub download_url {

    my ($input_basename, $content, $start_url) = @_;

    unless (defined($input_basename)) {
        croak 'got an undef basename';
    }
     
    my $download_url;

    my $root = HTML::TreeBuilder->new();

    $root->parse($content);
    $root->eof();

    my @tags = $root->look_down('_tag' => 'big');

    CONTENT_TAG:
    foreach my $tag (@tags) {

        my $content_ref = $tag->content();

        if (defined($content_ref)) {

            if(@{$content_ref} == 1) {

                if ($content_ref->[0] =~ /$input_basename/) { 

                    my @tmp = split(/\//, $content_ref->[0]);
                    my $ss_basename = pop @tmp;
                    $ss_basename .= '.gz';

                    $download_url = $start_url; 
  
                    $download_url .= '/BatchFileDownloadServlet?';
                    $download_url .= 'file=';
                    $download_url .= $ss_basename;
                    $download_url .= '&download=true';              
 
                    last CONTENT_TAG;
 
                }

            }
 
        }

    }

    $root->delete();
 
    return $download_url;

}

sub server_busy {

    my ($content) = @_;


    if ($content =~ qr/server load is high/i) {
        if ($content =~ qr/try again/i) {
            return 1;
        }    
    }
 
    return 0;

}


__END__

=head1 NAME

submit_to_SS.pl - Submit snpInExon_SSformat.txt file to SeattleSeq Annotation 

=head1 SYNOPSIS

submit_to_SS.pl --input-file snpInExon_SSformat.txt  --output-file snpInExon_annotated.txt

submit_to_SS.pl --email-address you@wustl.edu --gene-location NCBI  --input-file snpInExon_SSformat.txt  --output-file snpInExon_annotated.txt

submit_to_SS.pl -i snpInExon_SSformat.txt -o snpInExon_annotated.txt

submit_to_SS.pl -e you@wustl.edu -g CCDS -i snpInExon_SSformat.txt -o snpInExon_annotated.txt

submit_to_SS.pl --dbsnp-130 --protein-sequence -i snpInExon_SSformat.txt -o snpInExon_annotated.txt

Options:
  -help          brief help message

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-i>

Alias for --input-file

=item B<--input-file>

Specifies the input file to send to the SeattleSeq Annotation server.

=item B<-db>

Alias for --dbsnp-130

=item B<--dbsnp-130>

Specifies dbsnp v130 (ncbi36/hg18) instead of the default v131 (ncbi37/hg19).

=item B<-e>

Alias for email-address 

=item B<--email-address>

Specifies the email address to send to the SeattleSeq Annotation server

=item B<-o>

Alias for --output-file

=item B<--output-file>

Specifies the name of the file that the SeattleSeq annotation will be stored in.

=item B<-p>

Alias for --protein-sequence

=item B<--protein-sequence>

Specifies the option to include protein sequence in the annotation result (default is no protein sequence).

=item B<-gene-location>

Specifies whether to use NCBI or CCDS gene locations (or both).  Default is both.  For dbsnp-130, CCDS 2008 is used (sorry, no way to select CCDS 2007).  For dbsnp-131, CCDS2010 is used (the only option at the moment).

=item B<-g>

Alias for --gene-location

=head1 DESCRIPTION

B<submit_to_SS.pl> will POST the input file to the SeattleSeq Annotation Server.

=cut
