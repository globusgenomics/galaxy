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
string tool_cmd=@arg("tool_cmd");

foreach f,ix in items {
    string input_file_parts[] = @strsplit(f, "/");
    string input_basename = input_file_parts[@length(input_file_parts)-1];
    string run_output = @strcat(outputdir, "/", itemNames[ix]);
#    string output_vcf = @strcat(run_output, "/", itemNames[ix], ".vcf");
#    string final_vcf = @strcat(outputdir, "/../", itemNames[ix], ".vcf");
    file outlog <single_file_mapper; file=@strcat("logs/", input_basename, ".stdout.log")>;
    file errlog <single_file_mapper; file=@strcat("logs/", input_basename, ".stderr.log")>;

    string cmd1 = @regexp(@arg("tool_cmd"), "INPUTFILE", f);
    string cmd2 = @regexp(cmd1, "SAMPLENAME", itemNames[ix]);
    string cmd3 = @regexp(cmd2, "OUTPUTDIR", run_output);

    trace(cmd3);
    (outlog, errlog) = run_cmd(cmd3);
}
