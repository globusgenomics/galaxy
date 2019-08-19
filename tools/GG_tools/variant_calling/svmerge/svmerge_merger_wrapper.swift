## TYPE DECLARATIONS
type file;

## PROCEDURES
## execute an arbitrary command through the shell
app (file olog, file elog) job_runner(string cmd, file tarred) {
   bash "-c" cmd stdout=@olog stderr=@elog;
}


## PARSE COMMAND LINE ARGUMENTS
string filedir=@arg("input_dir", "Nonsense");

string app_cmd=@arg("app_cmd");
string search_grep=@arg("mapper_grep");

file inputfiles[] <ext; exec="/nfs/software/galaxy/tools/GG_tools/variant_calling/svmerge/file_mapper_svmerger.sh", location=filedir, grep_file=search_grep>;
string tarfile=@arg("input_tar", "Nonsense");

foreach f,i in inputfiles {

  file outlog <single_file_mapper; file=@strcat("logs/", @filename(f), ".stdout.log")>;
  file errlog <single_file_mapper; file=@strcat("logs/", @filename(f), ".stderr.log")>;
  file assemblyresults <single_file_mapper;file=tarfile>;

  ## use regular expressions to modify the command line
  string cmd0 = @regexp(@arg("app_cmd"), "(tarredFile)", @strcat(@filename(assemblyresults)));

  string cmd1 = @regexp(cmd0, "(sub_merge)", @strcat("/",@filename(f)));

  string basename = @strcut(@filename(f), ".*/(.*)$"); 
  #trace(@filename(f));
  #trace(basename);

  string outfile_flag = @strcat("outfile.", basename);
  #trace(outfile_flag);
  string cmd = @regexp(cmd1, "(outfile)", outfile_flag);
  
  trace(cmd); 

  ## send the constructed command through the shell
  (outlog, errlog) = job_runner(cmd, assemblyresults);
}
