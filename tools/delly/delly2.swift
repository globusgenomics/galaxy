## TYPE DECLARATIONS
type file;

## PROCEDURES
## execute an arbitrary command through the shell
app (file olog, file elog) run_cmd (string cmd) {
   bash "-c" cmd stdout=@olog stderr=@elog;
}

## PARSE COMMAND LINE ARGUMENTS
string input_list=@arg("inputfiles", "Nonsense");
string items[] = @strsplit(input_list, ",");

string samples_list=@arg("samplenames", "Nonsense");
string itemNames[] = @strsplit(samples_list, ",");

string outputdir=@arg("outputdir", "Nonsense");
string tool_cmd=@arg("tool_cmd", "Nonsense");
string env_cmd=@arg("envcmd", "Nonsense");
tracef("\ntool_cmd=%s\n",tool_cmd);
string tool_cmd_list[]=@strsplit(tool_cmd,";");

foreach tcl, iy in tool_cmd_list {
    string tmp[] =@strsplit(tcl," "); 
    string sv_type=tmp[6];
    foreach f,ix in items {
        string input_file_parts[] = @strsplit(f, "/");
        string input_basename = input_file_parts[@length(input_file_parts)-1];
        string run_output = @strcat(outputdir, "/", itemNames[ix], "_", sv_type, ".bcf");

        file outlog <single_file_mapper; file=@strcat("logs/", input_basename, "_", sv_type, ".stdout.log")>;
        file errlog <single_file_mapper; file=@strcat("logs/", input_basename, "_", sv_type, ".stderr.log")>;

#    string cmd1 = @regexp(@arg("tool_cmd"), "INPUTFILE", f);
        string cmd1 = @regexp(tcl, "INPUTFILE", f);
        string cmd2 = @regexp(cmd1, "OUTPUTDIR", run_output);
        string cmd3 = @strcat(env_cmd, " ", cmd2);
        trace(cmd3);
        (outlog, errlog) = run_cmd(cmd3);
    }
}
