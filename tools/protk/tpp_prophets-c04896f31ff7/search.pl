#!/usr/bin/perl -w

use strict;
use warnings;

open (IN, "<$ARGV[0]");
open (OUT_output ,">$ARGV[1]");
#open (OUT_info, ">$ARGV[2]");

my %options;
#reading all the parameters
while(<IN>){
    #print OUT_output $_;
    chop;
    if ($_=~ /^\*\+\*/){
	#print OUT $_."\n";
	if ($_=~ /searchParams\.([\w]+?)\="([^\"]+)"/){
	    #print OUT $1."+++".$2."\n";
	    $options{$1} = $2;
	}
    }
}

my $queryDir = "query";
my $libDir = "lib";
my $library = "test.splib";
my $outputFormat = "pep.xml";

#1. create essential links to files
if (exists $options{"splibType"}){
    if ($options{"splibType"} eq "single"){
	if (exists $options{"LibraryFile"} ){
	   
	    #my $splib = $options{"LibraryFile"};
	    #my $spidx = $options{"MzIdx"};
	    #my $pepidx = $options{"PepIdx"};
	    #my $link_splib = "temp.splib";
	    #my $link_spidx = "temp.spidx";
	    #my $link_pepidx = "temp.pepidx";
	    my $sptxt = $options{"LibraryFile"};
	    my $link_sptxt = "temp.sptxt";
	    #`mkdir $libDir;ln -s $splib $libDir/$link_splib;ln -s $spidx $libDir/$link_spidx;ln -s $pepidx $libDir/$link_pepidx`;
	    `mkdir $libDir; ln -s $sptxt $libDir/$link_sptxt;`;
	    `cd $libDir;spectrast -cNtest $link_sptxt;cd ..`;
	    #`cp lib/* /extra8/yhu/galaxy/example/test/`;
	    #print OUT_info `ls lib/`;
	}else{
	    die "ERROR: ONE of the file disappeared (lib,mzidx or pepidx)\n";
	}
    }elsif ($options{"splibType"} eq "fileSet"){
	if (exists $options{"LibraryFile"}){
	    my $zipLib = $options{"LibraryFile"};
	    my $link_zipLib = "tempLib.zip";
	    
	    `ln -s $zipLib $link_zipLib;unzip $link_zipLib -d $libDir/`;
	    #HERE we may need to check if the files all have the correct extensions
	    my $libraryName = `grep splib $libDir/fileSet.xml`;
	    if ($libraryName =~ /input=\"(.*?\.splib)\"/){
		$library = $1;
	    }else {
		my $display= `ls $libDir/`;
		die "ERROR: NO legal spectral library file exists in the unzipped fileSet!\n $display \n";
	    }
	}else{
	    die "ERROR: THE fileSet of library files disappeared.\n";
	}
    }else {
	die "ERROR: UNKNOWN file type of spectral libary: ".$options{"splibType"}."\n";
    }
}else{
    die "ERROR: NOT specify the file type of spectal libary (fileSet or single?)!\n";
}

if (exists $options{"queryType"}){
    if ($options{"queryType"} eq "single"){
	$outputFormat = "pep.xml";
	if ((exists $options{"queryFile"})
	    && (exists $options{"queryFileName"})){
	    my $mzxml = $options{"queryFile"};
	    my $link_mzxml = $options{"queryFileName"};
	    
	    `mkdir $queryDir`;
	    `ln -s $mzxml $queryDir/$link_mzxml`;
	    #fix the index of mzXML
            if ($options{'inputType'} eq "mzxml"){
	        `indexmzXML $queryDir/$link_mzxml`;
  	        `mv $queryDir/$link_mzxml.new $queryDir/$link_mzxml`;
            }
	}else{
	    die "ERROR: FILE or FileName of the query dataset disappeared\n";
	}
    }elsif ($options{"queryType"} eq "fileSet"){
	$outputFormat = "zip";
	if (exists $options{"queryFile"}){
	    my $zipQuery = $options{"queryFile"};
	    my $link_zipQuery = "tempQuery.zip";
	    
	    `ln -s $zipQuery $link_zipQuery; unzip $link_zipQuery -d $queryDir`;
	}else{
	    die "ERROR: THE fileSet of query dataset(s) disappeared!\n";
	}
    }else {
	die "ERROR: UNKNOWN file type of query dataset(s):".$options{"queryType"}."\n";
    }
}else {
    die "ERROR: NOT specify the file type of query dataset(s)\n";
}

#make directory for search results
my $database_search = $options{"search_database"};
my $outputDir = "output";
`mkdir $outputDir`;

my $command = "";
if ($options{'inputType'} eq "mzxml"){
    $command = "spectrast -sO$outputDir/ -sD$database_search -sL $libDir/$library $queryDir/*.mzXML";
}else {
    $command = "spectrast -sO$outputDir/ -sD$database_search -sL $libDir/$library $queryDir/*.mzML";
}

#print OUT_info $command;
my $run = `$command`;
#print OUT_output $command;
if(1){
    if ($outputFormat eq "pep.xml"){
	`mv $outputDir/*.pep.xml $outputDir/output.pep.xml`;
	print OUT_output `cat $outputDir/output.pep.xml`;
    }elsif($outputFormat eq "zip"){
    `cd $outputDir/; zip -r output.zip *; mv output.zip ../; cd ..`;
    print OUT_output `cat output.zip`;
    }else{
	die "ERROR: UNKNOWN output file format: $outputFormat!\n";
    }
}


#print OUT $run;


close IN;
#close OUT_info;
close OUT_output;
