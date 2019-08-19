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

foreach f,ix in itemNames {
    string tmp[]=@strsplit(f, " ");
    tracef("\n sampleName=%s\n",tmp[0]);
    string output_vcf = @strcat(outputdir, "/", tmp[0], ".vcf");
    tracef("\n output_vcf=%s\n",output_vcf);
    string final_vcf = @strcat(outputdir, "/../", tmp[0], ".vcf");
    tracef("\n final_vcf=%s\n",final_vcf);
    file outlog <single_file_mapper; file=@strcat("logs/", tmp[0], ".delly_to_metasv.stdout.log")>;
    file errlog <single_file_mapper; file=@strcat("logs/", tmp[0], ".delly_to_metasv.stderr.log")>;

    string cmd1 = @regexp(@arg("tool_cmd"), "INPUTFILE", output_vcf);
    string cmd2 = @regexp(cmd1, "OUTPUTFILE", final_vcf);

    trace(cmd2);
    (outlog, errlog) = run_cmd(cmd2);
}
