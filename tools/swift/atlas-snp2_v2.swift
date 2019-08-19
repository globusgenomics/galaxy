## TYPE DECLARATIONS
type file;
type commandline;

## PROCEDURES
app (file olog, file elog) atlas_snp2 (string bam_input, string fasta_input, string output, string samplename) {
   bash "-c" @strcat("/scratch/davidk/atlas-snp2/Atlas-SNP2.sh -i ", bam_input, " -r ", fasta_input, " -o ", output, " -s Illumina -v -n ", samplename) stdout=@olog stderr=@elog;
}

## PARSE COMMAND LINE ARGUMENTS
## Note: these arguments are passed directly from the python wrapper command line
## mandatory args
string fasta=@arg("fasta");
string outputdir=@arg("outputdir");
string bamdir=@arg("bam_dir");
## optional args
string target_region=@arg("target_region", "None");
string out_prefix=@arg("out_prefix", "None");
string platform=@arg("platform", "None");
string sample_name=@arg("sample_name", "None");
string post_cutoff=@arg("post_cutoff", "None");
string min_coverage=@arg("min_coverage", "None");
string prior_e=@arg("prior_e", "None");
string prior_l=@arg("prior_l", "None");
string base_sub_max=@arg("base_sub_max", "None");
string base_indel_max=@arg("base_indel_max", "None");
string insert_size_max=@arg("insert_size_max", "None");
string alignment_max=@arg("alignment_max", "None");
string input_sites=@arg("input_sites", "None");
string site_eval_flag=@arg("site_eval_flag", "None");

## Now we want to build up an Atlas command line argument
## initilize an empty array for command line arguments
## append arguments which are not "None", i.e. are passed in this run
## concatenate array into atlascmd within bam loop below
string[auto] atlascmd;

## add required commands
## TODO: find atlas jar or wrapper location
append(atlascmd, "/scratch/davidk/atlas-snp2/Atlas-SNP2.sh" );
append(atalscmd, " -r "); append(atlascmd, fasta);
if (platform != "None") {
  append(atlascmd,  " --"); append(atlascmd, platform);
}
if (post_cutoff != "None") {
  append(atlascmd,  " -c "); append(atlascmd, post_cutoff);
}
if (min_coverage != "None") {
  append(atlascmd,  " -y "); append(atlascmd, min_coverage);
}
if (prior_e != "None") {
  append(atlascmd,  " -e "); append(atlascmd, prior_e);
}
if (prior_l != "None") {
  append(atlascmd,  " -l "); append(atlascmd, prior_l);
}
if (base_sub_max != "None") {
  append(atlascmd,  " -m "); append(atlascmd, base_sub_max);
}
if (base_indel_max != "None") {
  append(atlascmd,  " -g "); append(atlascmd, base_indel_max);
}
if (insert_size_max != "None") {
  append(atlascmd,  " -p "); append(atlascmd, insert_size_max);
}
if (alignment_max != "None") {
  append(atlascmd,  " -f "); append(atlascmd, alignment_max);
}
if (input_sites != "None") {
  append(atlascmd,  " -a "); append(atlascmd, input_sites);
}
if (site_eval_flag != "None") {
  append(atlascmd,  " -w ");
}



file bamfiles[] <ext; exec="/scratch/davidk/atlas-snp2/bam_mapper.sh", location=bamdir>;

foreach f,i in bamfiles {
  string input_bam = @strcat("/", @f);
  string input_bam_parts[] = @strsplit(input_bam, "/");
  string input_bam_basename = input_bam_parts[@length(input_bam_parts)-1];

  file outlog <single_file_mapper; file=@strcat("logs/", input_bam_basename, ".stdout.log")>;
  file errlog <single_file_mapper; file=@strcat("logs/", input_bam_basename, ".stderr.log")>;
  string snp_outfile = @strcat(outputdir, "/", @regexp(input_bam_basename, ".bam", ""), ".snp");


  append(atlascmd, " -i "); append(input_bam);
  append(atlascmd, " -o "); append(snp_outfile);
  append(atlascmd, " -n "); append(input_bam_basename);
 
  trace(atlascmd);
 
  #(outlog, errlog) = atlas_snp2( input_bam, fasta, snp_outfile, input_bam_basename );
}

