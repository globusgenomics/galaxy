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

string outputdir=@arg("output_dir", "Nonsense");
string tool_cmd=@arg("tool_cmd");
tracef("\n tool_cmd: %s \n ", tool_cmd);

foreach f,ix in items {
    string output = @strcat(outputdir, "/../" , itemNames[ix], ".vcf");
    string tmpfile = @strcat(outputdir, "/" , itemNames[ix], ".tmp");
    string recodefile = @strcat(outputdir, "/" , itemNames[ix], ".tmp.recode.vcf");
    file outlog <single_file_mapper; file=@strcat("logs/", itemNames[ix], ".stdout.log")>;
    file errlog <single_file_mapper; file=@strcat("logs/", itemNames[ix], ".stderr.log")>;
    
    string cmd1 = @regexp(@arg("tool_cmd"), "INPUTFILE", f);
    string cmd2 = @regexp(cmd1, "OUTPUTFILE", output);
    string cmd3 = @regexp(cmd2, "TMPFILE", tmpfile);
    string cmd4 = @regexp(cmd3, "TMPFILE", recodefile);
    trace(cmd4);

    (outlog, errlog) = run_cmd(cmd4);
}

