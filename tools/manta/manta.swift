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
string tool_cmd=@arg("tool_cmd");

foreach f,ix in itemNames {
    file outlog <single_file_mapper; file=@strcat("logs/", f, ".stdout.log")>;
    file errlog <single_file_mapper; file=@strcat("logs/", f, ".stderr.log")>;

    string cmd1 = @regexp(@arg("tool_cmd"), "SAMPLEID", f);
    trace(cmd1);
    (outlog, errlog) = run_cmd(cmd1);
}
