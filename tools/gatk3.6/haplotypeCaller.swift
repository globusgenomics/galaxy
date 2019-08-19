## TYPE DECLARATIONS
type file;

## PROCEDURES
## execute an arbitrary command through the shell
app (file olog, file elog) gatk (string cmd) {
   bash "-c" cmd stdout=@olog stderr=@elog;
}


## PARSE COMMAND LINE ARGUMENTS
string bamdir=@arg("bam_dir", "Nonsense");
string outputdir=@arg("outputdir", "Nonsense");

string gatk_cmd=@arg("gatk_cmd");
#trace(gatk_cmd);

file bamfiles[] <ext; exec="/opt/galaxy/tools/swift/bam_mapper.sh", location=bamdir>;

foreach f,i in bamfiles {

  string final_atlascmd;

  string input_bam = @strcat("/", @f);
  string input_bam_parts[] = @strsplit(input_bam, "/");
  string input_bam_basename = input_bam_parts[@length(input_bam_parts)-1];

  file outlog <single_file_mapper; file=@strcat("logs/", input_bam_basename, ".stdout.log")>;
  file errlog <single_file_mapper; file=@strcat("logs/", input_bam_basename, ".stderr.log")>;
  string vcf_outfile = @strcat(outputdir, "/", @regexp(input_bam_basename, ".bam", ""), ".vcf");

  string cmd1 = @regexp(@arg("gatk_cmd"), "INPUT_BAM", input_bam);
  string cmd2 = @regexp(cmd1, "OUTPUT_VCF", vcf_outfile);
  trace(cmd2);

  ## send the constructed command through the shell
  (outlog, errlog) = gatk(cmd2);
}

