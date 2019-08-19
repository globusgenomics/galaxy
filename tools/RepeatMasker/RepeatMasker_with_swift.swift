## TYPE DECLARATIONS
type file;

## PROCEDURES
## execute an arbitrary command through the shell
app (file olog, file elog) job_runner(string cmd) {
   bash "-c" cmd stdout=@olog stderr=@elog;
}


## PARSE COMMAND LINE ARGUMENTS
string filedir=@arg("input_dir", "Nonsense");

string app_cmd=@arg("app_cmd");
#trace(filter_cmd);

file inputfiles[] <ext; exec="/nfs/software/galaxy/tools/RepeatMasker/file_mapper.sh", location=filedir>;

foreach f,i in inputfiles {

#  file input_file = @strcat("/", @f);
#  string input_file_parts[] = @strsplit(input_file, "/");
#  string input_file_basename = input_file_parts[@length(input_file_parts)-1];

  file outlog <single_file_mapper; file=@strcat("logs/", @filename(f), ".stdout.log")>;
  file errlog <single_file_mapper; file=@strcat("logs/", @filename(f), ".stderr.log")>;

  ## user regular expressions to modify the command line
  #string app_flag = @strcat( input_file);
  string cmd = @regexp(@arg("app_cmd"), "(--input)", @strcat("/",@filename(f)));
  trace(cmd); 

  # RepeatMasker -parallel 4 -nolow -noint -norna -species human -dir /scratch/uci/galaxy/files/000/dataset_54_files -gccalc -gff -html -s /scratch/uci/galaxy/files/000/dataset_28_files/consensus.chr1_gl000191_random.fasta
  ## send the constructed command through the shell
  (outlog, errlog) = job_runner(cmd);
}
