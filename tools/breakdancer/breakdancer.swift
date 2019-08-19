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
#tracef("\ntool_cmd=%s\n",tool_cmd);

foreach f,ix in items {
    string input_file_parts[] = @strsplit(f, "/");
    string input_basename = input_file_parts[@length(input_file_parts)-1];
#    tracef("\ninput_basename:%s\n",input_basename);
    string config = @strcat(outputdir, "/", itemNames[ix], "/", itemNames[ix], ".config");
    string output_raw = @strcat(outputdir, "/", itemNames[ix], "/", itemNames[ix], ".txt");
    string output_vcf = @strcat(outputdir, "/", itemNames[ix], "/", itemNames[ix], ".vcf");
#    tracef("\noutput_raw:%s\n",output_raw);
#    tracef("\noutput_vcf:%s\n",output_vcf);    
    file outlog <single_file_mapper; file=@strcat("logs/", input_basename, ".stdout.log")>;
    file errlog <single_file_mapper; file=@strcat("logs/", input_basename, ".stderr.log")>;

    string cmd1 = @regexp(@arg("tool_cmd"), "INPUTFILE", f);
    string cmd2 = @regexp(cmd1, "CONFIG1", config);
    string cmd3 = @regexp(cmd2, "CONFIG2", config);
    string cmd4 = @regexp(cmd3, "OUTPUT1", output_raw);
    string cmd5 = @regexp(cmd4, "OUTPUT2", output_raw);
    string cmd6 = @regexp(cmd5, "OUTPUT3", output_vcf);
    
    trace(cmd6);
    (outlog, errlog) = run_cmd(cmd6);
}
