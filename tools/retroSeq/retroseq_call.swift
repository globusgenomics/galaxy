## TYPE DECLARATIONS
type file;

## PROCEDURES
## execute an arbitrary command through the shell
app (file olog, file elog) job_runner(string cmd) {
   bash "-c" cmd stdout=@olog stderr=@elog;
}


## PARSE COMMAND LINE ARGUMENTS
string filedir=@arg("input_bam", "Nonsense");
string outputdir=@arg("output_dir", "Nonsense");

string app_cmd=@arg("app_cmd");
#trace(filter_cmd);

file inputchrs[] <ext; exec="/nfs/software/galaxy/tools/retroSeq/retroseq_call_file_mapper.sh", bam=filedir>;

foreach f,i in inputchrs {

  file outlog <single_file_mapper; file=@strcat("logs/", @filename(f), ".stdout.log")>;
  file errlog <single_file_mapper; file=@strcat("logs/", @filename(f), ".stderr.log")>;

  string cmd1 = @regexp(@arg("app_cmd"), "(CHR)", @strcat(@filename(f)));
  trace(cmd1); 

  string cmd = @regexp(cmd1, "(OUTPUT_NAME)", @strcat(outputdir, "/", @filename(f), ".retroseq_call"));
  trace(cmd);

  ## send the constructed command through the shell
  (outlog, errlog) = job_runner(cmd);
}

