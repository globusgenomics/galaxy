## TYPE DECLARATIONS
type file;

## PROCEDURES
## execute an arbitrary command through the shell
app (file olog, file elog) job_runner(string cmd) {
   bash "-c" cmd stdout=@olog stderr=@elog;
}


## PARSE COMMAND LINE ARGUMENTS
string filedir=@arg("input_dir", "Nonsense");
string outputdir=@arg("output_dir", "Nonsense");

string app_cmd=@arg("app_cmd");
#trace(filter_cmd);

file inputfiles[] <ext; exec="/nfs/software/galaxy/tools/tangram/tangram_scan_file_mapper.sh", location=filedir>;

foreach f,i in inputfiles {

  file outlog <single_file_mapper; file=@strcat("logs/", @filename(f), ".stdout.log")>;
  file errlog <single_file_mapper; file=@strcat("logs/", @filename(f), ".stderr.log")>;

  ## user regular expressions to modify the command line
  string basename = @strcut(@filename(f), ".*/(.*)$");
  string outputfile = @strcat(outputdir, "/", basename); 

  string cmd1 = @regexp(@arg("app_cmd"), "(INFILE)", @strcat("/",@filename(f)));
  trace(cmd1); 

  string cmd = @regexp(cmd1, "(OUTDIR)", @strcat("/",outputfile, ".scan_out"));
  trace(cmd);

  ## send the constructed command through the shell
  (outlog, errlog) = job_runner(cmd);
}

