## TYPE DECLARATIONS
type file;

## PROCEDURES
## execute an arbitrary command through the shell
app (file olog, file elog) run_cmd (string cmd) {
   bash "-c" cmd stdout=@olog stderr=@elog;
}


## PARSE COMMAND LINE ARGUMENTS
string parallel_list=@arg("parallel_list", "Nonsense");
string items[] = @strsplit(parallel_list, ",");

string outputdir=@arg("outputdir", "Nonsense");
string parallel=@arg("parallel", "Nonsense");
string tool_cmd=@arg("tool_cmd");

if (parallel == "chromosome") {
    foreach f in items {
        trace(f);
    }
} else {
    foreach f in items {
        string input_file_parts[] = @strsplit(f, "/");
        string input_basename = input_file_parts[@length(input_file_parts)-1];
        string outfile = @strcat(outputdir, "/", @regexp(input_basename, ".bed", ""), ".vcf");
        string tracefile = @strcat(outputdir, "/", @regexp(input_basename, ".bed", ""), ".trace");
        string failedallelesfile = @strcat(outputdir, "/", @regexp(input_basename, ".bed", ""), ".failed_alleles");
        file outlog <single_file_mapper; file=@strcat("logs/", input_basename, ".stdout.log")>;
        file errlog <single_file_mapper; file=@strcat("logs/", input_basename, ".stderr.log")>;

        string cmd1 = @regexp(@arg("tool_cmd"), "BED", f);
        string cmd2 = @regexp(cmd1, "OUTPUT", outfile);
        string cmd3 = @regexp(cmd2, "TRACE", tracefile);
        string cmd4 = @regexp(cmd3, "FAILED", failedallelesfile);

        trace(cmd4);
        (outlog, errlog) = run_cmd(cmd4);
    }
}
