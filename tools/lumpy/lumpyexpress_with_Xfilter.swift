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

#string samples_list=@arg("samplenames", "Nonsense");
#string itemNames[] = @strsplit(samples_list, ",");

string outputdir=@arg("outputdir", "Nonsense");
string tool_cmd=@arg("tool_cmd");

foreach f,ix in items {
    string input_file_parts[] = @strsplit(f, "/");
    string input_basename = input_file_parts[@length(input_file_parts)-1];
    string output_vcf = @strcat(outputdir, "/", input_basename, ".vcf");
    string output_tmp_dir = @strcat(outputdir, "/", input_basename, "_tmpdir/");
    string splitter_bam_unsorted = @strcat(outputdir, "/", input_basename, "_tmpdir/", input_basename, ".vcf.sample1.splitters.unsorted.bam" );
    string splitter_bam_sorted = @strcat(outputdir, "/", input_basename, "_tmpdir/", input_basename, ".vcf.sample1.splitters.sorted.bam" );
    string prefix_splitter_bam_sorted = @strcat(outputdir, "/", input_basename, "_tmpdir/", input_basename, ".vcf.sample1.splitters.sorted" );
    string discordant_bam_unsorted = @strcat(outputdir, "/", input_basename, "_tmpdir/", input_basename, ".vcf.sample1.discordant.unsorted.bam" );
    string discordant_bam_sorted = @strcat(outputdir, "/", input_basename, "_tmpdir/", input_basename, ".vcf.sample1.discordant.sorted.bam" );
    string prefix_discordant_bam_sorted = @strcat(outputdir, "/", input_basename, "_tmpdir/", input_basename, ".vcf.sample1.discordant.sorted" );
    string tmp_bai_file = @strcat(outputdir, "/", input_basename, "_tmpdir/", input_basename, ".vcf.sample1.lib1.splitters.sorted.bam.bai" );
    string bai_file = @strcat(outputdir, "/", input_basename, "_tmpdir/", input_basename, ".vcf.sample1.splitters.sorted.bam.bai" );
    string output_gt_vcf = @strcat(outputdir, "/", input_basename, ".gt.vcf");
    string output_gt1_vcf = @strcat(outputdir, "/", input_basename, ".gt.vcf");
    string output_final_vcf = @strcat(outputdir, "/../", input_basename, ".gt.sorted.vcf");
    file outlog <single_file_mapper; file=@strcat("logs/", input_basename, ".stdout.log")>;
    file errlog <single_file_mapper; file=@strcat("logs/", input_basename, ".stderr.log")>;

    string cmd1 = @regexp(@arg("tool_cmd"), "INPUTFILE_RAW", f);
    string cmd2 = @regexp(cmd1, "TEMPDIR", output_tmp_dir);
    string cmd3 = @regexp(cmd2, "OUTPUT", output_vcf);
    string cmd4 = @regexp(cmd3, "OUTPUT_BAM", f);
    string cmd5 = @regexp(cmd4, "SPLITTER_BAM_SORTED", splitter_bam_sorted);
    string cmd6 = @regexp(cmd5, "OUTPUT_VCF", output_vcf);
    string cmd7 = @regexp(cmd6, "OUTPUT_GT_VCF", output_gt_vcf);
    string cmd8 = @regexp(cmd7, "TMP_BAI_FILE", tmp_bai_file);
    string cmd9 = @regexp(cmd8, "BAI_FILE", bai_file);
    string cmd10 = @regexp(cmd9, "INPUT_GT_VCF1", output_gt1_vcf);
    string cmd11 = @regexp(cmd10, "INPUT_GT_VCF2", output_gt1_vcf);
    string cmd12 = @regexp(cmd11, "OUTPUT_FINAL_VCF", output_final_vcf);
    string cmd13 = @regexp(cmd12, "SPLITTER_BAM_UNSORTED", splitter_bam_unsorted);
    string cmd14 = @regexp(cmd13, "DISCORDANT_BAM_SORTED", discordant_bam_sorted);
    string cmd15 = @regexp(cmd14, "DISCORDANT_BAM_UNSORTED", discordant_bam_unsorted);
    string cmd16 = @regexp(cmd15, "INPUTFILE_RAW", f);
    string cmd17 = @regexp(cmd16, "INPUTFILE_RAW", f);
    string cmd18 = @regexp(cmd17, "DISCORDANT_BAM_UNSORTED", discordant_bam_unsorted);
    string cmd19 = @regexp(cmd18, "PREFIX_DISCORDANT", prefix_discordant_bam_sorted);
    string cmd20 = @regexp(cmd19, "SPLITTER_BAM_UNSORTED", splitter_bam_unsorted);
    string cmd21 = @regexp(cmd20, "PREFIX_SPLITTER", prefix_splitter_bam_sorted);
    string cmd22 = @regexp(cmd21, "SPLITTER_BAM_SORTED", splitter_bam_sorted);
    string cmd23 = @regexp(cmd22, "TEMPDIR", output_tmp_dir);
    string cmd24 = @regexp(cmd23, "SPLITTER_BAM_SORTED", splitter_bam_sorted);

    trace(cmd24);
    (outlog, errlog) = run_cmd(cmd24);
}
