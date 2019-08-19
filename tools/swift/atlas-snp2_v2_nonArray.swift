## TYPE DECLARATIONS
type file;

## PROCEDURES
## execute an arbitrary command through the shell
app (file olog, file elog) atlas_snp2 (string cmd) {
   bash "-c" cmd stdout=@olog stderr=@elog;
}


## PARSE COMMAND LINE ARGUMENTS
string bamdir=@arg("bam_dir", "Nonsense");
string outputdir=@arg("outputdir", "Nonsense");

string atlas_cmd=@arg("atlas_cmd");
#trace(atlas_cmd);

file bamfiles[] <ext; exec="/opt/galaxy/tools/swift/bam_mapper.sh", location=bamdir>;
file inputsamples[] <ext; exec="/opt/galaxy/tools/swift/bam_mapper_name.sh", bam=bamdir>;

foreach f,i in bamfiles {

  string final_atlascmd;

  string input_bam = @strcat("/", @f);
  string input_bam_parts[] = @strsplit(input_bam, "/");
  string input_bam_basename = input_bam_parts[@length(input_bam_parts)-1];

  file outlog <single_file_mapper; file=@strcat("logs/", input_bam_basename, ".stdout.log")>;
  file errlog <single_file_mapper; file=@strcat("logs/", input_bam_basename, ".stderr.log")>;
  string snp_outfile = @strcat(outputdir, "/", @regexp(input_bam_basename, ".bam", ""), ".snp");
  string complete_outfile = @strcat(outputdir, "/", @regexp(input_bam_basename, ".bam", ""), ".log");
  trace(complete_outfile);

  ## user regular expressions to modify the command line
  string bam_flag = @strcat("-i ", input_bam);
  string cmd1 = @regexp(@arg("atlas_cmd"), "(-i)", bam_flag);

  string sample_flag = @strcat("-n ", @inputsamples[i]);
  string cmd2 = @regexp(cmd1, "(-n)", sample_flag);
  
  string output_flag = @strcat("-o ", snp_outfile);
  string cmd3 = @regexp(cmd2, "(-o)", output_flag);

  string cmd4 = @regexp(cmd3, "OUTPUT", complete_outfile);
  trace(cmd4); 

  ## send the constructed command through the shell
  (outlog, errlog) = atlas_snp2(cmd4);
}

