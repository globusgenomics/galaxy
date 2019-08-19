## TYPE DECLARATIONS
type file;

## PROCEDURES
## execute an arbitrary command through the shell
app (file olog, file elog) run_cmd (string cmd) {
   bash "-c" cmd stdout=@olog stderr=@elog;
}

## PARSE COMMAND LINE ARGUMENTS
string samples_list=@arg("samplenames", "Nonsense");
string itemNames[] = @strsplit(samples_list, ",");

string outputdir=@arg("outputdir", "Nonsense");
string tool_cmd=@arg("tool_cmd", "Nonsense");
string env_cmd=@arg("envcmd", "Nonsense");
string tool_cmd_list[]=@strsplit(tool_cmd,";");

foreach tcl, iy in tool_cmd_list {
    string tmp[] =@strsplit(tcl," ");
    string sv_type=tmp[3];

    file outlog <single_file_mapper; file=@strcat("logs/", sv_type, ".stdout.log")>;
    file errlog <single_file_mapper; file=@strcat("logs/", sv_type, ".stderr.log")>;
    string cmd1 = @strcat(env_cmd, " ", tcl);
    trace(cmd1);
    (outlog, errlog) = run_cmd(cmd1);
}
