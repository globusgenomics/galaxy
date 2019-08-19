my ($day, $month, $year) = (localtime)[3,4,5];
$date = sprintf("%04d_%02d_%02d", $year+1900, $month+1, $day);

open (LOG,">log.out");

open (FILE,"/local/ngs/MTB_pipeline_files/Pipeline_shared_tools/Spoligotyper/NY_spolygo_types.txt");
@in=<FILE>;
close FILE;
$join = join ('',@in);

if ($join !~ $ARGV[2]) {
die "NY Code NOT found in the database!";
} 

system ("chmod 775 /local/ngs/MTB_pipeline_files/Pipeline_shared_tools/Spoligotyper/NY_spolygo_types.txt");
system ("cp /local/ngs/MTB_pipeline_files/Pipeline_shared_tools/Spoligotyper/NY_spolygo_types.txt /local/ngs/MTB_pipeline_files/Pipeline_shared_tools/Spoligotyper/NY_spolygo_types\_$date.txt");



@lines = split ('\n',$join);
open (OUT,">/local/ngs/MTB_pipeline_files/Pipeline_shared_tools/Spoligotyper/tmp_spolygo_types.txt");
foreach $in(@lines) {
if ($in =~ $ARGV[2]) {
@split = split ('\t',$in);
if (length($ARGV[1]) > 0) {
$split[0] = $ARGV[1];
}
if (length($ARGV[3]) > 0) {
$split[5] = $ARGV[3];
}
if (length($ARGV[4]) > 0) {
$split[6] = $ARGV[4];
}

print OUT join "\t",@split;
print OUT "\n";


} else {
	print OUT "$in\n";
}

}


close OUT;


print LOG "NY code $ny_code succesfully modify.";
close LOG;
system ("mv log.out $ARGV[0]");
system ("mv /local/ngs/MTB_pipeline_files/Pipeline_shared_tools/Spoligotyper/tmp_spolygo_types.txt /local/ngs/MTB_pipeline_files/Pipeline_shared_tools/Spoligotyper/NY_spolygo_types.txt");
system ("chmod 555 /local/ngs/MTB_pipeline_files/Pipeline_shared_tools/Spoligotyper/NY_spolygo_types.txt");
$ARGV[5] =~ s/\_\_at\_\_/\@/g;
system ("sendemail -t $ARGV[5] -m \"Please do not reply to this email.  For any questions or issues, contact the pipeline administrator at pascal.lapierre\@health.ny.gov\" -u \"MTB pipeline updated Spoligotype list\" -a /local/ngs/MTB_pipeline_files/Pipeline_shared_tools/Spoligotyper/NY_spolygo_types.txt -f bioinfo\@health.ny.gov");

