#!/usr/bin/perl -w

use strict;
use warnings;

open (IN, "<$ARGV[0]");
open (OUT, ">$ARGV[1]");

my %options;


while(<IN>){
    chop;
    if ($_ =~ /\*\+\*/){
	#print OUT $_."\n";
	if ($_=~ /params\.([\w]+?)="([^\"]+?)"/){
	    $options{$1} = $2;
	}
    }
}
#die "OK";

my $labelSrcFormat = "source";
my $srcFormat;
if (exists $options{$labelSrcFormat}){
    $srcFormat = $options{$labelSrcFormat};
    unless(($srcFormat eq "pepxml")
	   ||($srcFormat eq "sptxt")){
	die "ERROR: UNKNOWN format of source file: $srcFormat!\n";
    }
}else{
    die "ERROR: NOT specify the format of the source file!\n";
}

my $command;
if ($srcFormat eq "pepxml"){
    &cmd4Pepxml();
}elsif ($srcFormat eq "sptxt"){
    &cmd4Sptxt();
}else{
    die "ERROR: UNKNOWN format of source file: $srcFormat!\n";
}


#my $run = `$command`;
#print OUT $run;
#print OUT_splib `cat interact.splib`;
#print OUT_sptxt `cat interact.sptxt`;
#print OUT_spidx `cat interact.spidx`;
#print OUT_pepidx `cat interact.pepidx`;

close IN;
close OUT;
#close OUT_splib;
#close OUT_sptxt;
#close OUT_spidx;
#close OUT_pepidx;

sub cmd4Sptxt(){
    #setting source library file
    my $labelSourceFile = "sourceFile";
    my $srcFile;
    my $link_srcFile="temp.sptxt";
    my $tempLibName = "tmp";
    if (exists $options{$labelSourceFile}){
	$srcFile = $options{$labelSourceFile};
	`ln -s $srcFile $link_srcFile`;#link to temp.sptxt
	`spectrast -cN$tempLibName $link_srcFile`;
    }else{
	die "ERROR: NOT specify the source library file !\n";
    }
    
    #define the type of operation/manipulation
    my $labelOperationType = "operationType";
    my $operationType;
    if (exists $options{$labelOperationType}){
	$operationType = $options{$labelOperationType};
	if ($operationType eq "template"){
	}elsif ($operationType eq "advanced"){
	}elsif ($operationType eq "commonDecoy"){
	    my $labelDecoycc = "decoyConcatenate";
	    my $labelDecoycy = "decoySizeRatio";
	    my $decoycc;
	    my $decoycy;
	    if ((exists $options{$labelDecoycc})
		&&(exists $options{$labelDecoycy})){
		$decoycc = $options{$labelDecoycc};
		$decoycy = $options{$labelDecoycy};
		my $newLibName = $tempLibName."_decoy".$decoycy;
		my $cmd;
		if ($decoycc){
		    $cmd = "spectrast -cAD -cc -cy$decoycy -cN$newLibName $tempLibName.splib";
		}else{
		    $cmd = "spectrast -cAD -cy$decoycy -cN$newLibName $tempLibName.splib";
		}
		#my $run =`$cmd`;
		my $run =`$cmd > runinfo`;
		my $error = `grep error runinfo`;
		#print OUT $cmd;
		#print OUT `cat runinfo`;
		unless ($error =~ /without error/ ){
		    print OUT "ERROR:\n";
		    print OUT `cat runinfo`;
		}
		print OUT `cat $newLibName.sptxt`;
		return;
	    }else{
		die "ERROR: NOT specify decoyConcatenate and/or decoySizeRatio for decoy operation!\n";
	    }
	}else{
	    die "ERROR: UNKNOWN type of operation $operationType\n";
	}
    }else{
	die "ERROR: NOT specify the type of operations: template,common or advanced?\n";
    }

    
}


sub cmd4Pepxml(){

#1. create essential links to files
    my $labelQueryType = "queryType";
    my $labelProb = "prob" ;
    my $labelPepxml = "pepxml";
    my $labelQueryFile = "queryFile";
    my $labelQueryFileName = "queryFileName";
    
    my $prob;
    if (exists $options{$labelProb}){
	$prob = $options{$labelProb};
	eval {
	    if (($prob <= 1) && ($prob > 0)){
		return 1;
	    }else{
		die "ERROR: $prob is out of range [0:1]!\n";
	    }
	} or do{
	    if ($@){warn $@};
	}
    }else{
	die "ERROR: NO prob value is specified!\n";
    }
    
    my $queryFile;
    my $queryFileName;
    #my $queryDir = "query";
    if ((exists $options{$labelQueryFile})
	&& (exists $options{$labelQueryFile})){
	$queryFile = $options{$labelQueryFile};
	$queryFileName = $options{$labelQueryFileName};
	#print OUT `pwd`;
	#print OUT `ls`;
	`ln -s $queryFile $queryFileName`;
	#repair mzXML

	`indexmzXML $queryFileName`;
	`mv $queryFileName.new $queryFileName`;
    }else{
	die "ERROR: Query file can not be found or file name is lost!\n";
    }

    my $pepxml;
    my $link_pepxml = "interact.pep.xml";
    if (exists $options{$labelPepxml}){
	$pepxml = $options{$labelPepxml};
	`ln -s $pepxml $link_pepxml`;
    }else{
	die "ERROR: NOT specify the validated identification pep.xml!\n";
    }
    my $newLibName = "new";
    my $cmd = "spectrast -cN$newLibName -cP".$prob." $link_pepxml";
    my $run =`$cmd > runinfo`;
    my $error = `grep error runinfo`;
    unless ($error =~ /without error/ ){
	print OUT "ERROR:\n";
	print OUT `cat runinfo`;
    }
    print OUT `cat $newLibName.sptxt`;
}
