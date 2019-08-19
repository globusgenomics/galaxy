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

string tools_list=@arg("tools", "Nonsense");
string tools[] = @strsplit(tools_list, ",");
string tools_ix[];
foreach tool,ix in tools {
    #tracef("%s\n", tool);
    string parameter_input=@arg(tool, "Nonsense");
    tools_ix[ix] = parameter_input; 
    #tracef("%s\n", tools_ix[ix]);
}

string outputdir=@arg("outputdir", "Nonsense");
string inputbamdir=@arg("inputbamdir", "Nonsense");
string tool_cmd=@arg("tool_cmd");

foreach f,ix in itemNames {
    file outlog <single_file_mapper; file=@strcat("logs/", f, ".stdout.log")>;
    file errlog <single_file_mapper; file=@strcat("logs/", f, ".stderr.log")>;
    string input_basename = f;
    string bam_name = @strcat(inputbamdir,"/", f, ".bam");
    trace(bam_name);
    string output_sample_dir = @strcat(outputdir, "/", f);
    string[auto] tools_section;
    foreach tool,ix_tool in tools {
        #string tool_files[] = @strsplit(tools_ix[ix_tool], ",");
        #tracef("%s\n", tools_ix[ix_tool]);
        #tools_section << @strcat(tool, " ", tool_files[ix]);
        #tracef("%s\n", @strcat(tool, " ", tool_files[ix]));
        if (tool == "--breakdancer_native"){
            string tool_path = tools_ix[ix_tool];
            tools_section << @strcat(tool, " ", @strcat(tool_path, "/", f, ".txt"));
        } else {
            string tool_path = tools_ix[ix_tool];
            tools_section << @strcat(tool, " ", @strcat(tool_path, "/", f, ".vcf"));
        }
    }
    string tooly = @strjoin(tools_section, " ");
    string cmd1 = @regexp(@arg("tool_cmd"), "TOOLS", tooly);
    string cmd2 = @regexp(cmd1, "OUTDIR", output_sample_dir);
    string cmd3 = @regexp(cmd2, "SAMPLE", input_basename);
    string cmd4 = @regexp(cmd3, "BAM", bam_name);
    trace(cmd4);
    (outlog, errlog) = run_cmd(cmd4);
}
