#!/usr/bin/perl
 
use strict;
use warnings;
use Carp;
use Getopt::Long;

##############################################################################
# Check correct number of arguments given
##############################################################################
my $usage ="$0 - Converts annotated snpInExon and annotated indelsInExon tab-separated text files to XML Spreadsheet 2003 format.  The ExonicSNPs worksheet contains ALL exonic SNPs, the NovelSNPs worksheet contains only novel SNPs, and the Indels worksheet contains all exonic indels.  Novel snps are printed to a text file.

Usage: $0 --annotated-snp-in-exon-file <annotated snpInExon file> --annotated-indels-in-exon-file <annotated indelsInExon file> --novel-snps <novel snps file> --output-file <output file>
";

my $help;
my $snpInExon; # Worksheet 1
my $NovelSNPs; # Worksheet 2
my $Indels_filename; # Worksheet 3
my $output; # Output is in XML Spreadsheet 2003 format.  3/2/11 NR: For user friendliness, extension is renamed to .xls

GetOptions(
          "h|help|?"                          => \$help,
          "s|annotated-snp-in-exon-file=s"    => \$snpInExon,
          "i|annotated-indels-in-exon-file=s" => \$Indels_filename,
          "n|novel-snps=s"                    => \$NovelSNPs,
          "o|output-file=s"                   => \$output,
);

if ($help) {

	die $usage;
}

unless ((defined $snpInExon)       &&
        (defined $Indels_filename) &&
        (defined $output)) {

	croak "ERROR: incorrect number of arguments.\n\n$usage";
}


##############################################################################
# Set up
##############################################################################

# Check if files exist and keep track of which ones dont
my @files = ($snpInExon, $Indels_filename);
my @goahead = qw(1 1);
for (my $i = 0; $i< scalar(@files); $i++) {
	if (! -e $files[$i]) {
		print "\nWarning: $files[$i] does not exist and won't be converted to XML\n";
		$goahead[$i] = 0;
	}
}
		
##############################################################################
# Convert snpInExon_Annotated to Exonic and Novel SNP textfiles
##############################################################################

# If Database fields = "none," SNP is novel

# Open files
open (EXON, $snpInExon) or die $!;
open (NEW_EXON, ">",  $NovelSNPs) or die $!;

if (-s $snpInExon == 150) { # 150 = filesize snpInExon_Annotated.txt w/ just header and newline
	carp "Warning: annotated snp in exon file '$snpInExon' does not have any SNPs.  Therefore, there won't be any novel SNPs\n";
	# Copy header to $NovelSNPs.  This will allow a novel SNPs sheet to be made in the xls file
	unless (system("cp $snpInExon $NovelSNPs") == 0) {
		croak "Failed to 'cp $snpInExon $NovelSNPs'\n\n";
	}
} else {
	# Sort
	while (<EXON>) {
		my $text = $_;
		chomp $text;
		my @datafields = split(/\t/,$text);
		my $database = $datafields[15];
		
		if (($database eq "none") || ($. == 1)) { # Novel SNP or header
			print NEW_EXON join("\t",@datafields),"\n";
		} 
	}
}

# Close files
close (EXON) or die $!;
close (NEW_EXON) or die $!;		

		
##############################################################################
# Convert to XML
##############################################################################

# Pseudocode
# Add XML beginning lines
# For all worksheets
# 		Parse TAB seperated file
# 		Delete first line
# 		Add XML markup to eachline by printing to new file
# Add XML last lines


# Open files
open (OUTPUT, ">", $output) or die $!;


# Print beginning XML lines to OUTPUT
print OUTPUT 
'<?xml version="1.0"?>
<?mso-application progid="Excel.Sheet"?>
<Workbook xmlns="urn:schemas-microsoft-com:office:spreadsheet"
 xmlns:o="urn:schemas-microsoft-com:office:office"
 xmlns:x="urn:schemas-microsoft-com:office:excel"
 xmlns:ss="urn:schemas-microsoft-com:office:spreadsheet"
 xmlns:html="http://www.w3.org/TR/REC-html40">
 <DocumentProperties xmlns="urn:schemas-microsoft-com:office:office">
  <Author>GTAC</Author>
  <LastAuthor>GTAC</LastAuthor>
  <Created>2010-11-09T15:42:40Z</Created>
  <Version>12.00</Version>
 </DocumentProperties>
 <ExcelWorkbook xmlns="urn:schemas-microsoft-com:office:excel">
  <WindowHeight>14085</WindowHeight>
  <WindowWidth>28515</WindowWidth>
  <WindowTopX>240</WindowTopX>
  <WindowTopY>375</WindowTopY>
  <ProtectStructure>False</ProtectStructure>
  <ProtectWindows>False</ProtectWindows>
 </ExcelWorkbook>
 <Styles>
  <Style ss:ID="Default" ss:Name="Normal">
   <Alignment ss:Vertical="Bottom"/>
   <Borders/>
   <Font ss:FontName="Calibri" x:Family="Swiss" ss:Size="11" ss:Color="#000000"/>
   <Interior/>
   <NumberFormat/>
   <Protection/>
  </Style>
  <Style ss:ID="s62" ss:Name="Normal 2">
   <Alignment ss:Vertical="Bottom"/>
   <Borders/>
   <Font ss:FontName="Verdana" x:Family="Swiss"/>
   <Interior/>
   <NumberFormat/>
   <Protection/>
  </Style>
  <Style ss:ID="s64" ss:Parent="s62">
   <Alignment ss:Horizontal="Center" ss:Vertical="Bottom" ss:Rotate="75"
    ss:WrapText="1"/>
   <Borders>
    <Border ss:Position="Bottom" ss:LineStyle="Continuous" ss:Weight="1"
     ss:Color="#000080"/>
    <Border ss:Position="Left" ss:LineStyle="Continuous" ss:Weight="2"
     ss:Color="#000080"/>
    <Border ss:Position="Right" ss:LineStyle="Continuous" ss:Weight="2"
     ss:Color="#000080"/>
    <Border ss:Position="Top" ss:LineStyle="Continuous" ss:Weight="3"
     ss:Color="#000080"/>
   </Borders>
   <Font ss:FontName="Verdana" x:Family="Swiss" ss:Size="14" ss:Color="#000080"
    ss:Bold="1"/>
   <Interior/>
  </Style>
  <Style ss:ID="s65">
   <NumberFormat ss:Format="@"/>
  </Style>
  <Style ss:ID="s66">
   <NumberFormat ss:Format="d\-mmm"/>
  </Style>
  <Style ss:ID="s67" ss:Parent="s62">
   <Alignment ss:Horizontal="Center" ss:Vertical="Bottom" ss:Rotate="75"
    ss:WrapText="1"/>
   <Borders>
    <Border ss:Position="Bottom" ss:LineStyle="Continuous" ss:Weight="1"
     ss:Color="#000080"/>
    <Border ss:Position="Left" ss:LineStyle="Continuous" ss:Weight="2"
     ss:Color="#000080"/>
    <Border ss:Position="Right" ss:LineStyle="Continuous" ss:Weight="2"
     ss:Color="#000080"/>
    <Border ss:Position="Top" ss:LineStyle="Continuous" ss:Weight="3"
     ss:Color="#000080"/>
   </Borders>
   <Font ss:FontName="Verdana" x:Family="Swiss" ss:Size="14" ss:Color="#000080"
    ss:Bold="1"/>
   <Interior/>
   <NumberFormat ss:Format="@"/>
  </Style>
 </Styles>',"\n";

 
# Convert each textfile to XML and print to OUTPUT
# Exonic SNPs
if ($goahead[0]) {
	SNP('ExonicSNPs',$snpInExon);
	SNP('NovelSNPs',$NovelSNPs);
}

# Indels
if ($goahead[1]) {
	Indel('Indels',$Indels_filename);
}

# Print last line of XML to OUTPUT
print OUTPUT '</Workbook>',"\n";

##############################################################################
# Clean up
##############################################################################
# unlink ($NovelSNPs) or die $!;



##############################################################################
# Subroutines
##############################################################################

# SNP
#
# Converts SNP data to an XML worksheet
# Input:
#		$wksht_name = name of worksheet
#		$input		   = textfile name
 
sub SNP {

	my ($wksht_name,$input) = @_;
	
	 # Determine number of lines of $input file
	my ($num_lines,$num_cols) = count_lines_col($input);

	# Check if correct number of columns
	my $correct_num_cols = 16;
	if ($num_cols != $correct_num_cols) {
		print "\nWarning: $input does have the correct number of columns ($correct_num_cols) and won't be converted to XML\n";
		return; # break out of subroutine
	}

	# Print headers to OUTPUT
	print OUTPUT 
	 ' <Worksheet ss:Name="',$wksht_name,'">
	  <Table ss:ExpandedColumnCount="16" ss:ExpandedRowCount="',$num_lines,'" x:FullColumns="1"
	   x:FullRows="1" ss:DefaultRowHeight="15">
	   <Column ss:Index="2" ss:Width="52.5"/>
	   <Column ss:Index="5" ss:Width="48.75" ss:Span="4"/>
	   <Column ss:Index="14" ss:StyleID="s65" ss:AutoFitWidth="0"/>
	   <Row ss:AutoFitHeight="0" ss:Height="141">
		<Cell ss:StyleID="s64"><Data ss:Type="String">Chromosome</Data></Cell>
		<Cell ss:StyleID="s64"><Data ss:Type="String">Position</Data></Cell>
		<Cell ss:StyleID="s64"><Data ss:Type="String">Reference Base</Data></Cell>
		<Cell ss:StyleID="s64"><Data ss:Type="String">Alleles</Data></Cell>
		<Cell ss:StyleID="s64"><Data ss:Type="String">Consensus Quality</Data></Cell>
		<Cell ss:StyleID="s64"><Data ss:Type="String">SNP Quality</Data></Cell>
		<Cell ss:StyleID="s64"><Data ss:Type="String">Map Quality</Data></Cell>
		<Cell ss:StyleID="s64"><Data ss:Type="String">Read Depth</Data></Cell>
		<Cell ss:StyleID="s64"><Data ss:Type="String">Gene</Data></Cell>
		<Cell ss:StyleID="s64"><Data ss:Type="String">NCBI_id</Data></Cell>
		<Cell ss:StyleID="s64"><Data ss:Type="String">CCDS_id</Data></Cell>
		<Cell ss:StyleID="s64"><Data ss:Type="String">GVS_func</Data></Cell>
		<Cell ss:StyleID="s64"><Data ss:Type="String">Amino Acid</Data></Cell>
		<Cell ss:StyleID="s67"><Data ss:Type="String">Amino Acid Position</Data></Cell>
		<Cell ss:StyleID="s64"><Data ss:Type="String">polyPhen phenotype</Data></Cell>
		<Cell ss:StyleID="s64"><Data ss:Type="String">SNP Database</Data></Cell>
	   </Row>',"\n";
	 
	 
	 # Print data to OUTPUT
	 open (INPUT, "<", $input) or die $!; 
	 while(<INPUT>) {
		my $line = $_;
		chomp $line;
		

		if ($. == 1) { # Delete header line (starts at 1)
			next;
		}
		
		print OUTPUT '   <Row>',"\n"; 
		
		my @data = split(/\t/,$line);
		
		for (my $i = 0; $i < scalar(@data); $i++) {
			my $entry = $data[$i];
			
			if ($entry =~ m/^(\d+\.?\d*|\.\d+)$/) { # Is numeric
				print OUTPUT '    <Cell><Data ss:Type="Number">',$entry,'</Data></Cell>',"\n";
			
			} elsif ($i == 13) { # Column inported as text (looks like a number (or NA), so Excel would normally interpret as date
				print OUTPUT '    <Cell ss:StyleID="s62"><Data ss:Type="String">',$entry,'</Data></Cell>'."\n";
			
			} else { # Is char string
				# Remove quotes
				$entry =~ s/"//g;
				print OUTPUT '    <Cell><Data ss:Type="String">',$entry,'</Data></Cell>',"\n";
			}		
		}
		
		print OUTPUT '   </Row>',"\n";
		
	}
	
	wksht_end(); # Print rest of xml formatting text to OUTPUT
}


# Indel
#
# Converts indel data to an XML worksheet
# Input:
#		$wksht_name = name of worksheet
#		$input		   = textfile name

sub Indel {
	
	my ($wksht_name,$input) = @_;
	
	 # Determine number of lines of $input file
	my ($num_lines,$num_cols) = count_lines_col($input);
	
	# Check if correct number of columns
	my $correct_num_cols = 11;
	if ($num_cols != $correct_num_cols) {
		print "\nWarning: $input does have the correct number of columns ($correct_num_cols) and won't be converted to XML\n";
		return; # break out of subroutine
	}

	# Print headers to OUTPUT
	 print OUTPUT 
	 ' <Worksheet ss:Name="',$wksht_name,'">
	  <Table ss:ExpandedColumnCount="16" ss:ExpandedRowCount="',$num_lines,'" x:FullColumns="1"
	   x:FullRows="1" ss:DefaultRowHeight="15">
	   <Column ss:Index="2" ss:Width="52.5"/>
	   <Column ss:Index="5" ss:Width="48.75" ss:Span="4"/>
	   <Column ss:Index="14" ss:StyleID="s65" ss:AutoFitWidth="0"/>
	   <Row ss:AutoFitHeight="0" ss:Height="141">
		<Cell ss:StyleID="s64"><Data ss:Type="String">Chromosome</Data></Cell>
		<Cell ss:StyleID="s64"><Data ss:Type="String">Position</Data></Cell>
		<Cell ss:StyleID="s64"><Data ss:Type="String">Reference Base</Data></Cell>
		<Cell ss:StyleID="s64"><Data ss:Type="String">Alleles</Data></Cell>
		<Cell ss:StyleID="s64"><Data ss:Type="String">Consensus Quality</Data></Cell>
		<Cell ss:StyleID="s64"><Data ss:Type="String">SNP Quality</Data></Cell>
		<Cell ss:StyleID="s64"><Data ss:Type="String">Map Quality</Data></Cell>
		<Cell ss:StyleID="s64"><Data ss:Type="String">Coverage</Data></Cell>
		<Cell ss:StyleID="s64"><Data ss:Type="String">Gene</Data></Cell>
		<Cell ss:StyleID="s64"><Data ss:Type="String">Gene_id</Data></Cell>
		<Cell ss:StyleID="s64"><Data ss:Type="String">CCDS_id</Data></Cell>
	   </Row>',"\n";
	 
	 
	 # Print data to OUTPUT
	 open (INPUT, "<", $input) or die $!;
	 while(<INPUT>) {
		my $line = $_;
		chomp $line;

		if ($. == 1) { # Delete header line (starts at 1)
			next;
		}
		
		print OUTPUT '   <Row>',"\n"; #Preserve new lines?
		
		my @data = split(/\t/,$line);
		
		for (my $i = 0; $i < scalar(@data); $i++) {
			my $entry = $data[$i];
			
			if ($entry =~ m/^(\d+\.?\d*|\.\d+)$/) { # Is numeric
				print OUTPUT '    <Cell><Data ss:Type="Number">',$entry,'</Data></Cell>',"\n";
						
			} else { # Is char string
				# Remove quotes
				$entry =~ s/"//g;
				print OUTPUT '    <Cell><Data ss:Type="String">',$entry,'</Data></Cell>',"\n";
			}		
		}
		
		print OUTPUT '   </Row>',"\n";
		
	}
	
	wksht_end(); # Print rest of xml formatting text to OUTPUT
}


# count_lines_and_col
#
# Determines the number of lines and columns of a file
# Input:
#		$input 	   = textfile name
# Output:
#		$num_lines = number of lines of textfile
#		$num_cols = number of columns of textfile

sub count_lines_col {
	
	my($input) = @_;
	
	my $num_lines = 0;
	my $num_cols = 0;
	
	close INPUT; # Must reset line counter
	open (INPUT, "<", $input) or die $!;
	while(<INPUT>) {
		# Count number of columns		
		if ($. == 1) {
			my $line = $_;
			chomp $line;
			my @data = split(/\t/,$line);
			$num_cols = scalar(@data);
		}	
		
		# Count number of lines
		$num_lines++;
	}
	close INPUT;
	return ($num_lines,$num_cols);
}


# wksht_end
#
# Prints the necessary XML lines to OUTPUT to end a worksheet

sub wksht_end {

	# Print rest of XML format to OUTPUT
	print OUTPUT 
	'  </Table>
	  <WorksheetOptions xmlns="urn:schemas-microsoft-com:office:excel">
	   <PageSetup>
		<Header x:Margin="0.3"/>
		<Footer x:Margin="0.3"/>
		<PageMargins x:Bottom="0.75" x:Left="0.7" x:Right="0.7" x:Top="0.75"/>
	   </PageSetup>
	   <Unsynced/>
	   <Print>
		<ValidPrinterInfo/>
		<HorizontalResolution>600</HorizontalResolution>
		<VerticalResolution>600</VerticalResolution>
	   </Print>
	   <ProtectObjects>False</ProtectObjects>
	   <ProtectScenarios>False</ProtectScenarios>
	  </WorksheetOptions>
	 </Worksheet>',"\n";
}










 
