## TYPE DECLARATIONS
type file;

## PROCEDURES
## execute an arbitrary command through the shell
app (file olog, file elog) tandem_search (string cmd) {
   bash "-c" cmd stdout=@olog stderr=@elog;
}


## PARSE COMMAND LINE ARGUMENTS
string mzmldir=@arg("mzml_dir", "Nonsense");
string outputdir=@arg("outputdir", "Nonsense");

string tool_cmd=@arg("tool_cmd");
file mzmlfiles[] <ext; exec="/opt/galaxy/tools/protk/tpp_prophets-c04896f31ff7/file_mapper.sh", location=mzmldir>;

foreach f,i in mzmlfiles {
  string final_cmd;

  string input_mzml = @strcat("/", @f);
  string input_mzml_parts[] = @strsplit(input_mzml, "/");
  string input_mzml_basename = input_mzml_parts[@length(input_mzml_parts)-1];

  file outlog <single_file_mapper; file=@strcat("logs/", input_mzml_basename, ".stdout.log")>;
  file errlog <single_file_mapper; file=@strcat("logs/", input_mzml_basename, ".stderr.log")>;
  string mzml_outfile = @strcat(outputdir, "/", @regexp(input_mzml_basename, ".mzML.txt", ""), ".pepXML");

  ## user regular expressions to modify the command line
  string input_flag = @strcat(input_mzml);
  string cmd1 = @regexp(@arg("tool_cmd"), "(INPUT)", input_flag);

  string cmd2 = @regexp(cmd1, "(OUTPUT)", mzml_outfile);

  string cmd3 = @regexp(cmd2, "(OUTPUT)", mzml_outfile);

  trace(cmd3); 

  ## send the constructed command through the shell
  (outlog, errlog) = tandem_search(cmd3);
}

