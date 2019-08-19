my ($day, $month, $year) = (localtime)[3,4,5];
$date = sprintf("%04d_%02d_%02d", $year+1900, $month+1, $day);

open (LOG,">log.out");

open (FILE,"/scratch/galaxy/test/pxl10/Projects/MTB_pipeline_upgrade/Pipeline_shared_tools/Spoligotyper/NY_spolygo_types.txt");
@in=<FILE>;
close FILE;
$join = join ('',@in);
if ($join =~ $ARGV[2]) {
die "NY Code already in the database!";
} 
if ($join =~ $ARGV[3]) {
die "Octal code already in the database!";
} 
if (length($ARGV[3]) != 15) {
die "Octal code not 15 characters long as expected!";
}
if (length($ARGV[2]) != 6) {
die "Type code not 6 characters long as expected!";
}

system ("chmod 775 /scratch/galaxy/test/pxl10/Projects/MTB_pipeline_upgrade/Pipeline_shared_tools/Spoligotyper/NY_spolygo_types.txt");
system ("cp /scratch/galaxy/test/pxl10/Projects/MTB_pipeline_upgrade/Pipeline_shared_tools/Spoligotyper/NY_spolygo_types.txt /local/ngs/MTB_pipeline_files/Pipeline_shared_tools/Spoligotyper/NY_spolygo_types\_$date.txt");



open (OUT,">>/local/ngs/MTB_pipeline_files/Pipeline_shared_tools/Spoligotyper/NY_spolygo_types.txt");
print OUT "$ARGV[1]\t\t$ARGV[2]\t$ARGV[3]\t\t$ARGV[4]\t$ARGV[5]\n";
print LOG "NY code $ny_code succesfully added to the database.";
close LOG;
system ("mv log.out $ARGV[0]");
$ARGV[6] =~ s/\_\_at\_\_/\@/g;
system ("sendemail -t $ARGV[6] -m \"Please do not reply to this email.  For any questions or issues, contact the pipeline administrator at pascal.lapierre\@health.ny.gov\" -u \"MTB pipeline updated Spoligotype list\" -a /local/ngs/MTB_pipeline_files/Pipeline_shared_tools/Spoligotyper/NY_spolygo_types.txt -f bioinfo\@health.ny.gov");

system ("chmod 555 /local/ngs/MTB_pipeline_files/Pipeline_shared_tools/Spoligotyper/NY_spolygo_types.txt");

