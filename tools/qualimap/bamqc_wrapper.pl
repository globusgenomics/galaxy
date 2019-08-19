#!/usr/bin/perl
# bamqc_wrapper.pl
# Joachim Jacob - joachim.jacob@gmail.com - 2013

use strict;
use File::Temp 'tempdir';
use File::Basename;
use Log::Log4perl qw(:easy);

# ---------------------- Prepping Logging -----------------------------#
########################################################################
# Log levels: 			$DEBUG	$INFO	$WARN	$ERROR	$FATAL
# ConversionPattern:	%d %-5p %F{1} [%M] (line %L): %m%n%n
my $log_conf = q/ 
    log4perl.category = ERROR, Screen 
    log4perl.appender.Screen.stderr=1
    log4perl.appender.Screen.layout=Log::Log4perl::Layout::PatternLayout
    log4perl.appender.Screen.layout.ConversionPattern = %d %-5p %m%n
    log4perl.appender.Screen        = Log::Log4perl::Appender::Screen 
/;

Log::Log4perl::init( \$log_conf );
my $logger = get_logger();

# ----------------- Getting parameters file ---------------------------#
########################################################################
my ($configfile) = @ARGV;
my (%para);
open(CONFIG,"<$configfile");
while (<CONFIG>) {
	if (/(\S+)==(.+)$/){ $para{ $1 } = $2 ; }
}
close(CONFIG);

=Excerpt Config parameters

    ## first we pass some galaxy environment variables
        galtemp==${__new_file_path__}

        bam==$bam
        c==$c
        hm==$hm
        nr==$nr
        #if $customgtf.upload=="yes"
         gff==$customgtf.gff
         os==$customgtf.os
         p==$customgtf.p
        #end if

=cut

for my $para (keys %para){
	INFO "$para\tset to\t$para{$para}";
}


# ---------------------- Prepping temp dir ----------------------------#
########################################################################
# within the temporary directory of Galaxy, we create a temporary dir

my $galtemp = $para{'galtemp'};
DEBUG "\nReceived Galaxy temporary directory:\n$galtemp";

my $tempdir = File::Temp->tempdir('tmpXXXXX', DIR => $galtemp, CLEANUP => 1); 
mkdir "$tempdir", 077 unless -d "$tempdir";
INFO "\nTemporary directory:\n$tempdir";

# -------------------- Assembling command  ----------------------------#
########################################################################
my $outdir="bamqc_output";
my $qualimap_path="";						# for tool_dependency later on
my $command = $qualimap_path."qualimap bamqc --outdir $outdir"; # this will ultimately be executed

# based on parameters, which are called by $para{'parametername'}
# the command is assembled
$command .= " -bam $para{'bam'} -hm $para{'hm'} -nr $para{'nr'} -nt 32 --java-mem-size=128000M";
if($para{'c'}){$command .= " $para{'c'} "; }
if ( $para{'gff'} ){
	$command .= " -gff $para{'gff'} -p $para{'p'} ";
	if($para{'os'}) {
	    $command .= $para{'os'}; 
	}
} 
	
$command .= " -nt 8 -outformat HTML";


# --------------------- Executing command  ----------------------------#
########################################################################
run_process($command, "qualimap bamqc", ".");


# ----------------- Postprocessing command  ---------------------------#
########################################################################
# Hackish !!

my $librarydir = dirname($para{'bamqc_result'});
my $jobworkdir = dirname($para{'outputdir'});
my $files_path_name = basename($para{'outputdir'});
INFO "\n Librarydir = $librarydir \n jobworkdir = $jobworkdir \n files_path_name = $files_path_name \n";
system("mv $jobworkdir/bamqc_output $librarydir/$files_path_name") == 0 or die "Could not move output data to $librarydir\n";
system("rm -rf $para{'bamqc_result'}") == 0 or die "Could not post-process the results...\n";
system("ln -s $librarydir/$files_path_name/qualimapReport.html $para{'bamqc_result'}") == 0 or die "Could not move output data to $librarydir\n";


# --------------------------- Exiting  --------------------------------#
########################################################################
exit 0;



###	      Subroutines 					     ###
########################################################################
sub run_process {
	my ($command, $name, $tempdir)= @_;
	my $logger = get_logger();
	INFO "\nProcess to launch:\n $command\nIn directory\n$tempdir\n";
	system("cd $tempdir; $command 2>/dev/null") == 0 or die "$name failed\nExit status $?\nCommand: $command";
	if ($? == -1) {
		print "failed to execute: $!\n";
	} elsif ($? & 127) {
		printf "child died with signal %d, %s coredump\n", ($? & 127), ($? & 128) ? 'with' : 'without';
	} else {
		printf "$name executed successfully\n", $? >> 8;
	}	
}



