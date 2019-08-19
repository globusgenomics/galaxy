#!/usr/bin/python
import os, time, glob, optparse

def _vcf_header(sample_name):
    header = """##fileformat=VCFv4.1
##fileDate=%s
##source=SeqPilotV4.1.2
##INFO=<ID=TI,Number=.,Type=String,Description="Transcript ID">
##INFO=<ID=GI,Number=.,Type=String,Description="Gene ID">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership: ">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth at this position for this sample">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depth, unfiltered count of all reads">
##FORMAT=<ID=VF,Number=A,Type=Float,Description="Allele frequency for each ALT allele in the same order as listed">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	%s
""" % (str(time.strftime("%Y-%m-%d %H:%M")), sample_name)
    return header

# return the info field with given columns
# assume field is only: "GT:AD:VF:DP"
def _info_field(stored_values, unique_id, sample_id):
    # assumption: columns are: ["Zygosity", "Count", "Frequency", "Coverage"]
    columns = ["Zygosity", "Count", "Frequency", "Coverage"]
    zygozity = None
    count = stored_values[columns[1]][unique_id][sample_id]
    frequency = stored_values[columns[2]][unique_id][sample_id]
    coverage = stored_values[columns[3]][unique_id][sample_id]
    count_field = "%s,%s" % (str(int(coverage)-int(count)), str(count))
    if stored_values[columns[0]][unique_id][sample_id] == "Heterozygous":
        zygozity = "0/1"
    elif stored_values[columns[0]][unique_id][sample_id] == "Homozygous":
        zygozity = "1/1"
    field = "%s:%s:%s:%s" % (zygozity, count_field, str(frequency), str(coverage))
    return field

# convert annotated CLC VCF to regular VCF
def _convert_CLC_to_vcf(output_dir, annotations_dir):
    # get the necessary annotations from the raw annotated VCF files
    columns = ["Zygosity", "Count", "Frequency", "Coverage"]
    stored_values, sample_list = get_variant_values_from_raw_input(annotations_dir, columns)

    for anno_file in glob.glob("%s/*.txt" % annotations_dir):
        sample_name = os.path.basename(anno_file)
        modified_sample_name = sample_name 
        if " (paired)" in sample_name:
            values = sample_name.split(" ")
            modified_sample_name = values[0]
        vcf_fh = open(anno_file, "r")
        header = vcf_fh.readline()
        vcfout_fh = open("%s/%s.vcf" % (output_dir, modified_sample_name), "w")
        vcfout_fh.write(_vcf_header(modified_sample_name))
        info_field = "GT:AD:VF:DP"
        for line in vcf_fh:
            if line.startswith("#"):
                continue
            values = line.split("\t")
            chrm = values[0]
            ref = values[3]
            allele = values[4]
            pos = values[1]
            start = pos
            if values[2] == "Insertion":
                start = pos.split("^")[0]
                pos = start
            elif values[2] == "Deletion" or values[2] == "Replacement" or values[2] == "MNV":
                start = pos.split("..")[0]
                pos = start
            uniq_id = "%s_%s_%s" % (values[0], values[1], values[4])
            info_content = _info_field(stored_values, uniq_id, sample_name)
            new_line = "\t".join([chrm, pos, ".", ref, allele, ".", ".", ".", info_field, info_content])
            vcfout_fh.write(new_line + "\n")
        vcfout_fh.close()

# loop through anno files to get the rdf and avf values in a dictionary
def get_variant_values_from_raw_input(variant_annotation_dir, columns):
    anno_files = glob.glob("%s/*.txt" % variant_annotation_dir)
    # loop through anno files to get the rdf and avf values
    sample_list = []
    stored_values = {}
    for col in columns:
        stored_values[col] = {}
    #rdf_values = {}
    #avf_values = {}
    for anno in anno_files:
        sample_name = os.path.basename(anno)
        sample_list.append(sample_name)
        fh_anno = open(anno, "r")
        # assumption: each file's first line is a header
        # assumption: each row's first 5 columns are [chromosome, region, type, reference, allele]
        header = fh_anno.readline()
        column_indexes = match_column_names(header, columns)
        for line in fh_anno:
            line = line.rstrip("\n")
            values = line.split("\t")
            if ".." in values[1]:
                split_values = values[1].split("..")
                start = split_values[0]
            elif "^" in values[1]:
                split_values = values[1].split("^")
                start = split_values[0]
            else:
                start = values[1]
            #uniq_id = "%s_%s_%s_%s" % (values[0], start, values[3], values[4])
            uniq_id = "%s_%s_%s" % (values[0], values[1], values[4])
            for col in columns:
                if stored_values[col].get(uniq_id, None) is None:
                    stored_values[col][uniq_id] = {sample_name : values[column_indexes[col]]}
                else:
                    stored_values[col][uniq_id][sample_name] = values[column_indexes[col]]

    return [stored_values, sample_list]

# remove quote characters from string
def remove_quotes(value):
    value = value.rstrip("\"")
    value = value.lstrip("\"")
    return value

# return column index that matches the column names
def match_column_names(header, sample_list):
    header = header.rstrip("\n")
    header = header.rstrip("\r")
    values = [remove_quotes(x) for x in header.split("\t")]
    indexes = {}
    for i in sample_list:
        indexes[i] = values.index(i)
    return indexes

# return a header index matrix
def get_header_index(header):
    header_index = {}
    header = header.rstrip("\n")
    values = header.split("\t")
    count = 0
    for i in values:
        header_index[i] = count
        count += 1
    return header_index

def get_filepaths(directory):
    """
    This function will generate the file names in a directory
    tree by walking the tree either top-down or bottom-up. For each
    directory in the tree rooted at directory top (including top itself),
    it yields a 3-tuple (dirpath, dirnames, filenames).
    """
    file_paths = []  # List which will store all of the full filepaths.

    # Walk the tree.
    for root, directories, files in os.walk(directory):
        for filename in files:
            # Join the two strings in order to form the full filepath.
            filepath = os.path.join(".", filename)
            file_paths.append(filepath)  # Add it to the list.

    return file_paths  # Self-explanatory.

def __main__():
    descr = ""
    parser = optparse.OptionParser(description=descr)
    parser.add_option( '-d', '--dir', dest="output_dir", help="output_dir" )
    parser.add_option( '-o', '--output', dest="output", help="output" )
    parser.add_option( '-a', '--variant-annotations-directory', dest="variant_annotation_dir", help="if provided, the vcf annotations will be taken from here" )
    (options, args) = parser.parse_args()

    _convert_CLC_to_vcf(options.output_dir, options.variant_annotation_dir)
    directory_contents = get_filepaths(options.output_dir)
    #fh = open(options.output, "w")
    #fh.write("\n".join(directory_contents))
    #fh.close()
    for i in directory_contents:
        print "<p><a href='%s'>%s</a></p>" % (i, i)

if __name__ == "__main__":
    __main__()

