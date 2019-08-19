## TYPE DECLARATIONS
type file;

## PROCEDURES
## execute an arbitrary command through the shell
app (file olog, file elog) filter_vcf (string cmd) {
   bash "-c" cmd stdout=@olog stderr=@elog;
}


## PARSE COMMAND LINE ARGUMENTS
string vcfdir=@arg("vcf_dir", "Nonsense");
string outputdir=@arg("outputdir", "Nonsense");
string mapper_exec=@arg("mapper_exec", "Nonsense");

string filter_cmd=@arg("sort_cmd");
#trace(filter_cmd);

file vcffiles[] <ext; exec=mapper_exec, location=vcfdir>;

foreach f,i in vcffiles {

  string final_filtervcfcmd;

  string input_vcf = @strcat("/", @f);
  string input_vcf_parts[] = @strsplit(input_vcf, "/");
  string input_vcf_basename = input_vcf_parts[@length(input_vcf_parts)-1];

  file outlog <single_file_mapper; file=@strcat("logs/", input_vcf_basename, ".stdout.log")>;
  file errlog <single_file_mapper; file=@strcat("logs/", input_vcf_basename, ".stderr.log")>;
  string vcfsort_outfile = @strcat(outputdir, "/", @regexp(input_vcf_basename, ".vcf", ""), ".vcf");

  ## user regular expressions to modify the command line
  string vcf_flag = @strcat(input_vcf, " > ", vcfsort_outfile);
  string cmd1 = @regexp(@arg("sort_cmd"), "(--vcf-sort_file)", vcf_flag);
  trace(cmd1);

  ## send the constructed command through the shell
  (outlog, errlog) = filter_vcf(cmd1);
}

