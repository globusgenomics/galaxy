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

string files_list=@arg("inputfiles", "Nonsense");
string filepaths[] = @strsplit(files_list, ",");

string outputdir=@arg("outputdir", "Nonsense");
string tool_cmd=@arg("tool_cmd");

foreach f,ix in itemNames {
    file outlog <single_file_mapper; file=@strcat("logs/", f, ".stdout.log")>;
    file errlog <single_file_mapper; file=@strcat("logs/", f, ".stderr.log")>;
    string input_basename = f;
    string output_sample_dir = @strcat(outputdir, "/", f);
    string input_cns = @strcat(outputdir, "/", f, "/", input_basename, ".cns");
    string output_call_cns = @strcat(outputdir, "/", f, "/", input_basename, ".call.cns");
    string output_vcf = @strcat(outputdir, "/", f, "/", input_basename, ".vcf");
    string final_vcf = @strcat(outputdir, "/", input_basename, ".vcf");
    string mycnn = @strcat(outputdir, "/", f, "/", input_basename, ".cnn");

    string cmd1 = @regexp(@arg("tool_cmd"), "INPUTBAM", filepaths[ix]);
    string cmd2 = @regexp(cmd1, "OUTPUTDIR", output_sample_dir);
    string cmd3 = @regexp(cmd2, "INPUTCNS", input_cns);
    string cmd4 = @regexp(cmd3, "OUTPUTCALLFILE", output_call_cns);
    string cmd5 = @regexp(cmd4, "INPUTCALLFILE", output_call_cns);
    string cmd6 = @regexp(cmd5, "OUTPUTVCF", output_vcf);
    string cmd7 = @regexp(cmd6, "INT_VCF", output_vcf);
    string cmd8 = @regexp(cmd7, "FINAL_VCF", final_vcf);
    trace(cmd8);
    (outlog, errlog) = run_cmd(cmd8);
}
