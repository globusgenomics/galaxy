## TYPE DECLARATIONS
type file;

## PROCEDURES
## execute an arbitrary command through the shell
app (file olog, file elog) job_runner(string cmd) {
   bash "-c" cmd stdout=@olog stderr=@elog;
}


## PARSE COMMAND LINE ARGUMENTS
string filedir=@arg("input_dir", "Nonsense");

string app_cmd=@arg("app_cmd");
#trace(filter_cmd);

file inputfiles[] <ext; exec="/nfs/software/galaxy/tools/GG_tools/variant_calling/forestSV/file_mapper.sh", location=filedir>;

foreach f,i in inputfiles {

  file outlog <single_file_mapper; file=@strcat("logs/", @filename(f), ".stdout.log")>;
  file errlog <single_file_mapper; file=@strcat("logs/", @filename(f), ".stderr.log")>;

  ## user regular expressions to modify the command line
  string cmd = @regexp(@arg("app_cmd"), "(info_txt)", @strcat("/",@filename(f)));
  trace(cmd); 

  ## send the constructed command through the shell
  (outlog, errlog) = job_runner(cmd);
}
