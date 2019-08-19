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

string filter_cmd=@arg("filter_cmd");
#trace(filter_cmd);

file vcffiles[] <ext; exec=mapper_exec, location=vcfdir>;

foreach f,i in vcffiles {

  string final_filtervcfcmd;

  string input_vcf = @strcat("/", @f);
  string input_vcf_parts[] = @strsplit(input_vcf, "/");
  string input_vcf_basename = input_vcf_parts[@length(input_vcf_parts)-1];

  file outlog <single_file_mapper; file=@strcat("logs/", input_vcf_basename, ".stdout.log")>;
  file errlog <single_file_mapper; file=@strcat("logs/", input_vcf_basename, ".stderr.log")>;
  string vcftools_outfile = @strcat(outputdir, "/", @regexp(input_vcf_basename, ".vcf", ""), ".tmp");
  string vcffixup_outfile = @strcat(outputdir, "/", @regexp(input_vcf_basename, ".vcf", ""), ".tmp2");
  string vcfrecode_outfile = @strcat(outputdir, "/", @regexp(input_vcf_basename, ".vcf", ""), ".recode.vcf");


  ## user regular expressions to modify the command line
  string vcf_flag = @strcat("--vcf ", input_vcf);
  string cmd1 = @regexp(@arg("filter_cmd"), "(--vcf)", vcf_flag);

  string output_flag = @strcat("--out ",  vcftools_outfile);
  string cmd2 = @regexp(cmd1, "(--out)", output_flag);

  string vcffixup_flag = @strcat("vcffixup ",  vcftools_outfile, ".recode.vcf", " > ", vcffixup_outfile);
  string cmd3 = @regexp(cmd2, "(vcffixup)", vcffixup_flag);

  string vcffilter_flag = @strcat( vcffixup_outfile, " > ", vcfrecode_outfile, "; rm -rf ", outputdir, "/", @regexp(input_vcf_basename, ".vcf", ""), "*.tmp*");
  string cmd4 = @regexp(cmd3, "(--vcffilter_file)", vcffilter_flag);


  trace(cmd4); 

  ## send the constructed command through the shell
  (outlog, errlog) = filter_vcf(cmd4);
}

