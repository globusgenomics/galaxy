#!/usr/bin/python
import pybedtools, glob
import os, sys, tempfile, optparse, re
from numpy import sum
import time

### Written by Alex Rodriguez
### Navipoint

# column mapping in all_variants file
all_var_cols = {"chrom_col" : 0, "start" : 1, "stop" : 2, "strand" : 3, "thousandG_matrix" : 4,
                "unique_id" : 5, "exac_matrix" : 6, "coding_matrix" : 7, "gene_matrix" : 8, "location_matrix" : 9,
                "avf_matrix" : 10, "min_avf_index" : 11, "min_avf_value" : 12, "avf_above_zero" : 13, "rdf_matrix" : 14,
                "dist_matrix" : 15}

# get a value from a matrix
def valid_value_of(matrix):
    for i in matrix:
        if i != ".":
            return i
    return "."

# get all column values for an ID
def get_all_values_for_variant(unique_id, stored, header_index, sample_list, exclude_columns=[]):
    values = []
    for key in sorted(header_index, key=lambda i: int(header_index[i])):
        if key in exclude_columns:
            continue
        column_values = []
        for sample in sample_list:
            if stored[key][unique_id].get(sample, None) is not None:
                column_values.append(stored[key][unique_id].get(sample))
        values.append(";".join(set(column_values)))
    return values

# print out final format of output
def _create_final_output(final_variants, final_file, sample_list, stored, header):
    header_index = get_header_index(header)
    exclude_columns = ["Allele Frequency Allele 1", "Allele Frequency Allele 2", "Read depth (infoDP)", "Allelic Depth Allele 1", "Allelic Depth Allele 2"]

    header = header.rstrip("\n")
    new_header = []
    for i in header.split("\t"):
        if i not in exclude_columns:
            new_header.append(i)
    header = "\t".join(new_header)

    final_fh = open(final_file, "w")
    final_fh.write("%s\t%s\t%s\n" % ("Uniq_variants", header, "\t".join(sample_list)))
    for variant in final_variants:
        unique_id = variant[all_var_cols["unique_id"]]
        avf_values = variant[all_var_cols["avf_matrix"]].split(",")
        rdf_values = variant[all_var_cols["rdf_matrix"]].split(",")
        paired_values = []
        index = 0
        for avf in avf_values:
            paired_values.append("%s/%d" % (str(float(avf)*100), float(rdf_values[index])))
            index += 1

        values = get_all_values_for_variant(unique_id, stored, header_index, sample_list, exclude_columns=exclude_columns)
        final_fh.write(unique_id + "\t" + "\t".join(values) + "\t" + "\t".join(paired_values) + "\n")
    final_fh.close()

# write final filtered VCF files given the variants filtered by samples
def _create_final_vcf_output(sample_list, sample_variants, output_dir, stored, vcf_dir): 
    final_variants = []
    for variant in sample_variants:
        final_variants.append(variant[all_var_cols["unique_id"]])
        # break the unique id to sections 
        # and add a new variant id that takes out the Reference variant value
        values = variant[all_var_cols["unique_id"]].split("_")
        final_variants.append("_".join([values[0], values[1], values[3]]))

    # go through each VCF in the vcf_dir and remove variants that are not in sample_variants
    for vcf in glob.glob("%s/*" % vcf_dir):
        vcf_fh = open(vcf, "r")
        vcfout_fh = open("%s/%s" % (output_dir, os.path.basename(vcf)), "w")
        for line in vcf_fh:
            if line.startswith("#"):
                vcfout_fh.write(line)
            else:
                values = line.split("\t")
                unique_id = "%s_%s_%s_%s" % (values[0], values[1], values[3], values[4])
                unique_id = unique_id.replace("-", ".")
                unique_id_homo = "%s_%s_%s" % (values[0], values[1], values[4])
                unique_id_homo = unique_id_homo.replace("-", ".")
                if unique_id in final_variants or unique_id_homo in final_variants:
                    vcfout_fh.write(line)
        vcfout_fh.close()
    return 0

# loop through anno files to get the rdf and avf values in a dictionary
def get_variant_values_from_raw_input(variant_annotation_dir, columns):
    anno_files = glob.glob("%s/*.txt" % variant_annotation_dir)
    # loop through anno files to get the rdf and avf values
    sample_list = []
    stored_values = {}
    variant_dict = {}
    for col in columns:
        stored_values[col] = {}
    for anno in anno_files:
        sample_name = os.path.basename(anno)
        sample_list.append(sample_name)
        fh_anno = open(anno, "r")
        # assumption: each file's first line is a header
        # assumption: each row's first 5 columns are [chromosome, region, type, reference, allele]
        for line in fh_anno:
            if not line.startswith("##"):
                header = line
                break
        header_index = get_header_index(header)
        column_indexes = match_column_names(header, columns)
        for line in fh_anno:
            line = line.rstrip("\n")
            values = line.split("\t")
            uniq_id = "%s_%s_%s_%s" % (values[header_index["Chromosome"]], values[header_index["Start"]], values[header_index["Allele 1"]], values[header_index["Allele 2"]])
            #uniq_id = "%s_%s_%s_%s" % (values[header_index["Chromosome"]], values[header_index["Start"]], values[header_index["Reference"]], values[header_index["Allele 2"]])
            variant_dict[uniq_id] = variant_dict.get(uniq_id,0) + 1
            for col in columns:
                if values[column_indexes[col]] == "":
                    values[column_indexes[col]] = "."
                if stored_values[col].get(uniq_id, None) is None:
                    stored_values[col][uniq_id] = {sample_name : values[column_indexes[col]]}
                else:
                    stored_values[col][uniq_id][sample_name] = values[column_indexes[col]]
    return [stored_values, sample_list, variant_dict]

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

# return matrix of sample values
def get_value_matrix(anno_values, sample_list, var_id):
    matrix = []
    for i in sample_list:
        if anno_values[var_id].get(i) == ".":
            anno_values[var_id][i] = 0
        matrix.append(float(anno_values[var_id].get(i,0)))
    return matrix

# return matrix of sample string values
def get_string_matrix(anno_values, sample_list, var_id):
    matrix = []
    for i in sample_list:
        matrix.append(str(anno_values[var_id].get(i,".")))
    return matrix

# gdt the indexes in an AVF string where the avf value is greater than threshold
def get_index_list_above_theshold(string, threshold):
    index_list = []
    values = [float(i) for i in string.split(",")]
    index = 0
    for value in values:
        if value >= threshold and value != 0:
            index_list.append(index)
        index = index+1
    return index_list

def subtract_variants(bed_A, bed_B):
    variants_to_remove = []
    for variant in bed_B:
        variants_to_remove.append("\t".join(variant.fields))
    final_variants = []
    for variant in bed_A:
        if "\t".join(variant.fields) not in variants_to_remove:
            final_variants.append("\t".join(variant.fields))
    fd, path = tempfile.mkstemp()
    tmp_bed = open(path, "w")
    tmp_bed.write("\n".join(final_variants))
    tmp_bed.close()
    final_bed = pybedtools.BedTool(path).sort()
    return final_bed

def get_max_sample_avf_cluster(cluster, feature, index):
    max_value = None
    for neighbor in cluster:
        if feature != neighbor:
            max_sample_value = get_index_value_from_sample(neighbor[all_var_cols["avf_matrix"]], index)
            if max_value is None:
                max_value = [max_sample_value, index, neighbor]
            elif max_sample_value >= max_value[0]:
                max_value = [max_sample_value, index, neighbor]
    return max_value

def get_max_cluster_avf(cluster, feature):
    max_value = [0, 0]
    for neighbor in cluster:
        if feature != neighbor:
            max_sample_value, index = get_max_from_numerical_string(neighbor[all_var_cols["avf_matrix"]])
            if max_sample_value >= max_value[0]:
                max_value = [max_sample_value, index, neighbor]
    return max_value

def get_max_from_numerical_string(string):
    values = [float(i) for i in string.split(",")]
    max_value = max(values)
    index = values.index(max_value)
    return [max_value, index] 

def get_min_above_threshold_from_numerical_string(string, threshold):
    values = [float(i) for i in string.split(",")]
    min_above_threshold = min(i for i in values if i >= threshold)
    if min_above_threshold is None:
        min_above_threshold = 0
    if min_above_threshold == 0:
        index = 0
    else:
        index = values.index(min_above_threshold)
    return [min_above_threshold, index]

def get_index_value_from_sample(string, index):
    values = [float(i) for i in string.split(",")]
    return values[index]

def clusters_to_list(clusters):
    # get cluster dictionary
    cluster_dict = {}
    id_to_cluster = {}
    for feature in clusters:
        id_to_cluster[feature.fields[all_var_cols["unique_id"]]] = feature.fields[-1]
        if cluster_dict.get(feature.fields[-1]):
            cluster_dict[feature.fields[-1]].append(feature.fields[0:-1])
        else:
            cluster_dict[feature.fields[-1]] = [feature.fields[0:-1]]
    return cluster_dict, id_to_cluster

def remove_avf_abundant_variants(bed, window=10, occurrence_threshold=1.0):
    clusters = bed.cluster(d=window).sort()
    # get cluster dictionary
    cluster_dict, id_to_cluster = clusters_to_list(clusters)

    variants_to_remove = []
    for cluster_id in cluster_dict:
        cluster = cluster_dict[cluster_id]
        sample_count = len(cluster[0][all_var_cols["avf_matrix"]].split(","))
        occurrence_count = 0
        cluster_sample_values = []
        for feature in cluster:
            cluster_sample_values.append([float(i) for i in feature[all_var_cols["avf_matrix"]].split(",")])

        # get the complete count of all samples in a cluster which have more than a x% precense
        sums = sum(cluster_sample_values, axis=0)
        occurrence_avf_value = sum(x > 0 for x in sums)
        if float(occurrence_avf_value)/float(sample_count) > float(occurrence_threshold):
            variants_to_remove.append("\n".join("\t".join(map(str,i)) for i in cluster))

    fd, path = tempfile.mkstemp()
    tmp_bed = open(path, "w")
    tmp_bed.write("\n".join(variants_to_remove))
    tmp_bed.close()
    remove = pybedtools.BedTool(path).sort()
    #final_variants = bed.subtract(remove, A=True).sort()
    final_variants = subtract_variants(bed, remove)
    #os.remove(path)
    return final_variants

def get_multiple_clusters(bed, window):
    clusters = bed.cluster(d=window).sort()
    tmp_clusters_bed = tempfile.NamedTemporaryFile()
    cluster_counts_bed = clusters.groupby(c=len(all_var_cols)+1, g=[len(all_var_cols)+1], o=['count'], output=tmp_clusters_bed.name)
    cluster_counts = get_dict(tmp_clusters_bed.name)
    multiple_clusters = subset_clusterCount(clusters, (1), '<', cluster_counts, len(clusters[0].fields)-1).cut(range(len(all_var_cols))).sort()
    return multiple_clusters

def remove_low_variance_in_window(bed, occurrence_threshold, avf_threshold, cluster_window):
    multiple_clusters = get_multiple_clusters(bed, cluster_window)
    clusters = multiple_clusters.cluster(d=cluster_window).sort()
    cluster_dict, id_to_cluster = clusters_to_list(clusters)

    low_occurrence_bed = subset_featuretypes(multiple_clusters, (str(occurrence_threshold[0])), occurrence_threshold[1], operator="==" ).sort()
    low_occurrence_avf_bed = subset_featureMath(low_occurrence_bed, (str(avf_threshold[0])), '<=', avf_threshold[1]).sort()
   
    # search through low_occurrence_avf_bed for variants that are less than 5% AVF in matching sample cluster 
    variants_to_remove = []
    for variant in low_occurrence_avf_bed:
        unique_id = variant.fields[all_var_cols["unique_id"]]
        cluster_id = id_to_cluster[unique_id]
        avf_value = variant.fields[avf_threshold[1]]
        avf_sampleID = variant.fields[avf_threshold[1]-1]
        flag = None
        for feature in cluster_dict[cluster_id]:
            feature_unique_id = feature[all_var_cols["unique_id"]]
            if feature_unique_id != unique_id:
                neighbor_avf = get_index_value_from_sample(feature[all_var_cols["avf_matrix"]], int(avf_sampleID))
                if float(neighbor_avf) >= avf_threshold[0]:
                    flag = 1
                    break
        if flag is None:
            variants_to_remove.append("\t".join(variant.fields))
    if len(variants_to_remove) > 0:
        fd, path = tempfile.mkstemp()
        tmp_bed = open(path, "w")
        tmp_bed.write("\n".join(variants_to_remove))
        tmp_bed.close()
        remove = pybedtools.BedTool(path).sort()
        #final_variants = bed.subtract(remove).sort()
        final_variants = subtract_variants(bed, remove).sort()
        #os.remove(path)
        return final_variants
    else:
        return bed

def get_dict(bed_file):
    bed_dict = {}
    fh_bed = open(bed_file, "r")
    for line in fh_bed:
        line = line.rstrip("\n")
        values = line.split("\t")
        bed_dict[values[0]] = int(values[1])
    fh_bed.close()
    return bed_dict

def clusterCount_filter(feature, countThreshold, operator, cluster_dict, field_index):
    if operator == "==":
        if countThreshold == cluster_dict[feature[field_index]]:
            return True
    elif operator == "<=":
        if countThreshold <= cluster_dict[feature[field_index]]:
            return True
    elif operator == "<":
        if countThreshold < cluster_dict[feature[field_index]]:
            return True
    elif operator == ">=":
        if countThreshold >= cluster_dict[feature[field_index]]:
            return True
    elif operator == ">":
        if countThreshold > cluster_dict[feature[field_index]]:
            return True
    return False

def subset_clusterCount(bed, countThreshold, operator, cluster_dict, field_index):
    result = bed.filter(clusterCount_filter, countThreshold, operator, cluster_dict, field_index).saveas()
    return pybedtools.BedTool(result.fn)

def featuretype_filter(feature, featuretype, field_index, operator=None):
    if operator is None:
        if featuretype in feature[field_index]:
            return True
    elif operator == "in":
        if featuretype in feature[field_index]:
            return True
    elif operator == "==":
        if featuretype == feature[field_index]:
            return True
    return False

def subset_featuretypes(gff, featuretype, field_index, operator=None):
    result = gff.filter(featuretype_filter, featuretype, field_index, operator=operator).saveas()
    return pybedtools.BedTool(result.fn)

def featureMath_filter(feature, featureValue, operator, field_index):
    if operator == "<=":
        if float(feature[field_index]) <= float(featureValue):
            return True
    elif operator == "<":
        if float(feature[field_index]) < float(featureValue):
            return True
    elif operator == ">=":
        if float(feature[field_index]) >= float(featureValue):
            return True
    elif operator == ">":
        if float(feature[field_index]) > float(featureValue):
            return True
    elif operator == "==":
        if float(feature[field_index]) == float(featureValue):
            return True
    return False

def subset_featureMath(bed, featureValue, operator, field_index):
    result = bed.filter(featureMath_filter, featureValue, operator, field_index)
    return pybedtools.BedTool(result.fn) 

def featureMax_filter(feature, featureValue, field_index):
    max_sample_avf, sample_index = get_max_from_numerical_string(feature[field_index])
    if float(max_sample_avf) <= float(featureValue):
        return True
    return False

def subset_featuresall_under5(bed, featureValue, field_index):
    result = bed.filter(featureMax_filter, featureValue, field_index)
    return pybedtools.BedTool(result.fn)

def thresholdvalue_filter(feature, percent_threshold, min_rdf_threshold, max_rdf_threshold):
    if float(feature[all_var_cols["avf_above_zero"]]) == 1:
        sample_avf = float(feature[all_var_cols["min_avf_value"]])
        sample_index = int(feature[all_var_cols["min_avf_index"]])
        sample_rdf = get_index_value_from_sample(feature[13], sample_index)
    else:
        sample_avf, sample_index = get_min_above_threshold_from_numerical_string(feature[all_var_cols["avf_matrix"]], 5)
        sample_rdf = get_index_value_from_sample(feature[all_var_cols["rdf_matrix"]], sample_index)
    if sample_rdf < float(max_rdf_threshold) and sample_rdf >= float(min_rdf_threshold) and sample_avf < float(percent_threshold):
        return True
    return False

def subset_variant_percent(bed, percent_threshold, min_rdf_threshold, max_rdf_threshold):
    result = bed.filter(thresholdvalue_filter, percent_threshold, min_rdf_threshold, max_rdf_threshold).saveas()
    return pybedtools.BedTool(result.fn)

# test a set of avf and rdf conditions
def _check_all_samples_satisfy_conditions(conditions, feature):
    sample_index_list = get_index_list_above_theshold(feature[all_var_cols["avf_matrix"]], 0)
    flag_list = []
    
    # test that for each sample it satisfies each condition
    # if a sample does not satisfy all of the conditions, return False
    # if all samples satisfy one of the conditions, return True
    for sample_index in sample_index_list:
        sample_avf = get_index_value_from_sample(feature[all_var_cols["avf_matrix"]], sample_index)
        sample_rdf = get_index_value_from_sample(feature[all_var_cols["rdf_matrix"]], sample_index)

        flag_list.append(_check_sample_satisfy_conditions(conditions, [sample_avf, sample_rdf]))
    
    # loop through flag_list. If all elements are True, then return True
    # If any element is False, then return False
    for i in flag_list:
        if i is False:
            return False
    return True

# given a sample_rdf and sample_avf, check that it satisfies a condition
def _check_sample_satisfy_conditions(conditions, sample_values):
    flag = None
    for condition in conditions:
        max_rdf_threshold = condition[2]
        min_rdf_threshold = condition[1]
        percent_threshold = condition[0]
        if sample_values[1] < max_rdf_threshold and sample_values[1] >= min_rdf_threshold and sample_values[0] < percent_threshold:
            return True
    return False

def _get_sample_neighbor_values(variant, cluster, index):
    neighbor_values = []
    for feature in cluster:
        if feature[all_var_cols["unique_id"]] != variant[all_var_cols["unique_id"]]:
            neighbor_values.append([get_index_value_from_sample(feature[all_var_cols["avf_matrix"]], index), get_index_value_from_sample(feature[all_var_cols["rdf_matrix"]], index)])
    return neighbor_values

#def subset_by_rdf_avf(bed, percent_threshold, min_rdf_threshold, max_rdf_threshold, cluster_window):
def subset_by_rdf_avf(bed, conditions, cluster_window):
    clusters = bed.cluster(d=cluster_window).sort()
    cluster_dict, id_to_cluster = clusters_to_list(clusters)
    variants_to_remove = []
    for cluster_id in cluster_dict:
        cluster = cluster_dict[cluster_id]
        if len(cluster) == 1:
            if _check_all_samples_satisfy_conditions(conditions, cluster[0]):
                variants_to_remove.append("\t".join(cluster[0]))
        else:
            for variant in cluster:
                flag = None
                # first check if each sample in variant satisfies a condition
                if _check_all_samples_satisfy_conditions(conditions, variant):
                    # for each feature in cluster check that it satisfies one of the conditions in the same sample
                    sample_index_list = get_index_list_above_theshold(variant[all_var_cols["avf_matrix"]], 0)
                    for index in sample_index_list:
                        sample_neighbors = _get_sample_neighbor_values(variant, cluster, index)
                        for neighbor in sample_neighbors:
                            if not _check_sample_satisfy_conditions(conditions, neighbor):
                                flag = 1
                                break
                    if flag is None:
                        variants_to_remove.append("\t".join(variant))

    if len(variants_to_remove) > 0:
        fd, path = tempfile.mkstemp()
        tmp_bed = open(path, "w")
        tmp_bed.write("\n".join(variants_to_remove))
        tmp_bed.close()
        remove = pybedtools.BedTool(path)
        #os.remove(path)
        return remove
    else:
        return variants_to_remove

def distance_value_filter(feature, window):
    if float(feature[-1]) <= window:
        return True
    return False

def subset_closest(bed_with_distance_in_last_column, window):
    result = bed_with_distance_in_last_column.filter(distance_value_filter, window).saveas()
    return pybedtools.BedTool(result.fn)

def __main__():
    descr = ""
    parser = optparse.OptionParser(description=descr)
    parser.add_option( '-v', '--vcf', dest="vcf_input", help='vcf' )
    parser.add_option( '-d', '--dir', dest="output_dir", help="output_dir" )
    parser.add_option( '', '--vcf-output-list', dest="output_dir_contents_file", help="output_dir_contets_file" )
    parser.add_option( '-g', '--gtf', dest="gtf_input", help='gtf' )
    parser.add_option( '-r', '--rdf', dest="rdf_input", help='read_depth_input' )
    parser.add_option( '-o', '--output', dest="output_bed", help="bed_output" )
    parser.add_option( '-e', '--exclude', dest="exclude_bed", help="exclude bed" )
    parser.add_option( '-w', '--window', dest="window", help="window size for clusters" )
    parser.add_option( '-t', '--population-threshold', dest="population_threshold", help="Population variant occurrence percentage threshold for a given cluster (0,1)" )
    parser.add_option( '-a', '--variant-annotations-directory', dest="variant_annotation_dir", help="if provided, the vcf annotations will be taken from here" )
    parser.add_option( '', '--vcf-directory', dest="vcf_dir", help="if provided, the raw vcf files will be printed out" )
    parser.add_option( '-s', '--sample-list', dest="sample_list", help="text file with exact name of sample name as it appears on input matrix annotated file" )
    parser.add_option( '-m', '--min-avf', dest="min_avf", help="Minimum AVF threshold for a given variant (0,1)" )
    parser.add_option( '', '--debug', dest="debug", action="store_true", help="By running in DEBUG mode, it will print out all intermediate files in the outputs directory" )
    (options, args) = parser.parse_args()

    if len(args) == 1:
        function = args[0]
    else:
        sys.exit("You will need to specify if you are converting (convert) from CLC raw to VCF or filtering variant output (filter)")
    if options.debug:
        DEBUG = True
    else:
        DEBUG = False

    if not options.window:
        cluster_window = 20
    else:
        cluster_window = options.window
    if not options.population_threshold:
        population_threshold = 0.75
    else:
        population_threshold = options.population_threshold
    if not options.min_avf:
        min_avf = 0.05
    else:
        min_avf = options.minimum_avf


    #conditions = [[1000, 0, 60], [20, 60, 100], [10, 100, 200], [5, 200, 100000]]
    conditions = [[1, 0, 60], [.10, 60, 200], [.05, 200, 1000000]]

    # is variant annotation directory provided
    if options.variant_annotation_dir:
        # get the column names from one of the VCF in the directory
        anno_files = glob.glob("%s/*" % options.variant_annotation_dir)
        anno_file = open(anno_files[0], "r")
        for line in anno_file:
            if not line.startswith("##"):
                header = line
                break
        header_index = get_header_index(header)
        anno_file.close()
        columns = header_index.keys()
        #columns = ["Allele Frequency Allele 2", "Read depth (infoDP)"]
        stored_values, sample_list, variant_dict = get_variant_values_from_raw_input(options.variant_annotation_dir, columns)
        rdf_values = stored_values["Read depth (infoDP)"]
        avf_values = stored_values["Allele Frequency Allele 2"]

    # create matrix with info from all samples
    # create a tmp BED file
    # Need to put the chr, start and end locations at the beggnining of the line
    vcf_tmp, vcf_tmp_filename = tempfile.mkstemp()
    vcf_fh = open(vcf_tmp_filename, "w")
    for unique_id in sorted(variant_dict):
        avf_matrix = get_value_matrix(avf_values, sample_list, unique_id)
        rdf_matrix = get_value_matrix(rdf_values, sample_list, unique_id)
        chrm = valid_value_of(get_string_matrix(stored_values["Chromosome"], sample_list, unique_id))
        start = valid_value_of(get_string_matrix(stored_values["Start"], sample_list, unique_id))
        stop = valid_value_of(get_string_matrix(stored_values["Stop"], sample_list, unique_id))
        thousandG_matrix = valid_value_of(get_string_matrix(stored_values["1000 genomes phase 3 allele frequency all (1000Gp3_FreqAll)"], sample_list, unique_id))
        if "=" in str(thousandG_matrix):
            ratio, thousandG_matrix = thousandG_matrix.split("=")
        if thousandG_matrix == ".":
            thousandG_matrix = 0
        exac_matrix = valid_value_of(get_string_matrix(stored_values["ExAC allele frequency all (ExAC_alleleFrequencyAll)"], sample_list, unique_id))
        if "=" in str(exac_matrix):
            ratio, exac_matrix = exac_matrix.split("=")
        if exac_matrix == ".":
            exac_matrix = 0
        coding_matrix = valid_value_of(get_string_matrix(stored_values["Effect (codingEffect)"], sample_list, unique_id))
        gene_matrix = valid_value_of(get_string_matrix(stored_values["Gene (gene)"], sample_list, unique_id))
        location_matrix = valid_value_of(get_string_matrix(stored_values["Location (varLocation)"], sample_list, unique_id))
        dist_matrix = valid_value_of(get_string_matrix(stored_values["Positions from nearest splice site (distNearestSS)"], sample_list, unique_id))

        if max(avf_matrix) > 0:
            m = min(i for i in avf_matrix if i > 0)
            count = sum(x > 0 for x in avf_matrix)
            selected_sampleId = avf_matrix.index(m)
            sample_name = sample_list[selected_sampleId]
            rdv = rdf_values[unique_id][sample_name]
        else:
            m = 0.0
            count = 0
            selected_sampleId = "0"
            rdv = "0"

        vcf_fh.write("%s\n" % ("\t".join([chrm, str(start), str(stop), ".", str(thousandG_matrix), unique_id, str(exac_matrix), coding_matrix, gene_matrix, location_matrix, ",".join(map(str, avf_matrix)), str(selected_sampleId), str(m), str(count), ",".join(map(str, rdf_matrix)), dist_matrix])))
    vcf_fh.close()

    # create a BEDtools object
    if DEBUG is True:
        variants = pybedtools.BedTool(vcf_tmp_filename).sort().saveas('%s/all_variants.bed' % options.output_dir)
    else:
        variants = pybedtools.BedTool(vcf_tmp_filename).sort()

    # step 1.1
    # get exonic
    if DEBUG is True:
        variants_exonic = subset_featuretypes(variants, ('exonic'), all_var_cols["location_matrix"], operator="in" ).saveas('%s/exonic.bed' % options.output_dir)
    else:
        variants_exonic = subset_featuretypes(variants, ('exonic'), all_var_cols["location_matrix"], operator="in" )
    # step 1.2
    # get splcing
    if DEBUG is True:
        variants_intronic = subset_featuretypes(variants, ('intronic'), all_var_cols["location_matrix"], operator="in" ).saveas('%s/intronic.bed' % options.output_dir)
        variants_splicing_tmp = subset_featureMath(variants_intronic, (2), '<=', all_var_cols["dist_matrix"]).sort().saveas('%s/splicing2.bed' % options.output_dir)
        variants_splicing = subset_featureMath(variants_splicing_tmp, (-2), '>=', all_var_cols["dist_matrix"]).sort().saveas('%s/splicing_2.bed' % options.output_dir)
    else:
        variants_intronic = subset_featuretypes(variants, ('intronic'), all_var_cols["location_matrix"], operator="in" )
        variants_splicing_tmp = subset_featureMath(variants_intronic, (2), '<=', all_var_cols["dist_matrix"]).sort()
        variants_splicing = subset_featureMath(variants_splicing_tmp, (-2), '>=', all_var_cols["dist_matrix"]).sort()

    # step 2
    # merge 1.1, 1.2
    if DEBUG is True:
        step2_variants = variants_exonic.cat(variants_splicing, postmerge=False).sort().saveas('%s/step2.bed' % options.output_dir)
    else:
        step2_variants = variants_exonic.cat(variants_splicing, postmerge=False).sort()

    # step 3
    # remove variants where occurrence count is 1 and AVF is less than 5% and there are no variants within 20nt
    clusters = step2_variants.cluster(d=cluster_window).sort()
    tmp_clusters_bed = tempfile.NamedTemporaryFile()
    cluster_counts_bed = clusters.groupby(c=len(all_var_cols)+1, g=[len(all_var_cols)+1], o=['count'], output=tmp_clusters_bed.name)
    cluster_counts = get_dict(tmp_clusters_bed.name)
    single_clusters = subset_clusterCount(clusters, (1), '==', cluster_counts, len(clusters[0].fields)-1).cut(range(len(all_var_cols))).sort()
    step3_occurrence1 = subset_featuretypes(single_clusters, ('1'), all_var_cols["avf_above_zero"], operator="==" ).sort()
    if DEBUG is True:
        step3_occurrence1_avfunder5 = subset_featureMath(step3_occurrence1, (min_avf), '<', all_var_cols["min_avf_value"]).sort().saveas("%s/step3_avfunder5_occurrence1.bed" % options.output_dir)
        step3_variants = subtract_variants(step2_variants, step3_occurrence1_avfunder5).sort().saveas("%s/step3.bed" % options.output_dir)
    else:
        step3_occurrence1_avfunder5 = subset_featureMath(step3_occurrence1, (min_avf), '<', all_var_cols["min_avf_value"]).sort()
        step3_variants = subtract_variants(step2_variants, step3_occurrence1_avfunder5).sort()

    # step 4
    # remove variants where occurrence is 1 and AVF <5% and if for variants within 20nt it is also <5% AVF or 0
    occurrence_threshold = [1, all_var_cols["avf_above_zero"]] # tuple [occurence_threhsold, index_field_in_bed]
    avf_threshold = [min_avf, all_var_cols["min_avf_value"]] # tuple [avf_threhsold, index_field_in_bed]
    if DEBUG is True:
        step4_variants = remove_low_variance_in_window(step3_variants, occurrence_threshold, avf_threshold, cluster_window).sort().saveas("%s/step4.bed" % options.output_dir)
    else:
        step4_variants = remove_low_variance_in_window(step3_variants, occurrence_threshold, avf_threshold, cluster_window).sort()

    # step 5
    # Remove if >1 sample with <5% AVF (but none above 5%) and none of the same samples have a second variant within 20nt >5% AVF
    clusters = step4_variants.cluster(d=cluster_window).sort()
    tmp_clusters_bed = tempfile.NamedTemporaryFile()
    cluster_counts_bed = clusters.groupby(c=len(all_var_cols)+1, g=[len(all_var_cols)+1], o=['count'], output=tmp_clusters_bed.name)
    cluster_counts = get_dict(tmp_clusters_bed.name)
    cluster_dict, id_to_cluster = clusters_to_list(clusters)
    occurrence_gt1 = subset_featureMath(clusters, (1), '>', all_var_cols["avf_above_zero"]).cut(range(len(all_var_cols))).sort()
    if DEBUG is True:
        occurrence_gt1_all_avfunder5 = subset_featuresall_under5(occurrence_gt1, (min_avf), all_var_cols["avf_matrix"]).sort().saveas('%s/step5_occurrence_gt1_all_avfunder5.bed' % options.output_dir)
    else:
        occurrence_gt1_all_avfunder5 = subset_featuresall_under5(occurrence_gt1, (min_avf), all_var_cols["avf_matrix"]).sort()

    # search through low_occurrence_avf_bed for variants that are less than 5% AVF in matching sample cluster 
    variants_to_remove = []
    for variant in occurrence_gt1_all_avfunder5:
        unique_id = variant.fields[all_var_cols["unique_id"]]
        cluster_id = id_to_cluster[unique_id]
        sample_index_list = get_index_list_above_theshold(variant[all_var_cols["avf_matrix"]], 0)
        flag = None
        for feature in cluster_dict[cluster_id]:
            feature_unique_id = feature[all_var_cols["unique_id"]]
            if feature_unique_id != unique_id:
                for index in sample_index_list:
                    feature_avf_value = get_index_value_from_sample(feature[all_var_cols["avf_matrix"]], index) 
                    if float(feature_avf_value) >= avf_threshold[0]:
                        flag = 1
                        break
                if flag == 1:
                    break
        if flag is None:
            variants_to_remove.append("\t".join(variant.fields))
    if len(variants_to_remove) > 0:
        fd, path = tempfile.mkstemp()
        tmp_bed = open(path, "w")
        tmp_bed.write("\n".join(variants_to_remove))
        tmp_bed.close()
        if DEBUG is True:
            remove = pybedtools.BedTool(path).sort().saveas("%s/step5_remove.bed" % options.output_dir)
            step5_variants = subtract_variants(step4_variants, remove).sort().saveas("%s/step5.bed" % options.output_dir)
        else:
            remove = pybedtools.BedTool(path).sort()
            step5_variants = subtract_variants(step4_variants, remove).sort()
    else:
        if DEBUG is True:
            step5_variants = step4_variants.sort().saveas("%s/step5.bed" % options.output_dir)
        else:
            step5_variants = step4_variants.sort()

    # step 6 and 7
    # If read depth is <60x, variants are filtered out.
    # each conditions is a triple ==> [percent_threshold, min_rdf_threshold, max_rdf_threshold]
    #variants_to_remove = subset_by_rdf_avf(step5_variants,conditions, cluster_window).sort().saveas("%s/remove_step7.bed" % options.output_dir)
    if DEBUG is True:
        variants_to_remove = subset_by_rdf_avf(step5_variants,conditions, cluster_window).sort().saveas("%s/step7_remove.bed" % options.output_dir)
    else:
        variants_to_remove = subset_by_rdf_avf(step5_variants,conditions, cluster_window).sort()
    if len(variants_to_remove) > 0:
        if DEBUG is True:
            step7_variants = subtract_variants(step5_variants, variants_to_remove.sort()).sort().saveas("%s/step7.bed" % options.output_dir)
        else:
            step7_variants = subtract_variants(step5_variants, variants_to_remove.sort()).sort()
    else:
        if DEBUG is True:
            step7_variants = step5_variants.sort().saveas("%s/step7.bed" % options.output_dir)
        else:
            step7_variants = step5_variants.sort()

    # step 8
    # Remove CTBP2, KDM6B, DCP1B, PDGFRA, KMT2C
    # exclude BED file
    # BEDTOOLS exclude object
    if DEBUG is True:
        exclude = pybedtools.BedTool(options.exclude_bed).sort().saveas("%s/step8_exclude.bed" % options.output_dir)
        intersect_s8 = step7_variants.window(exclude, w=1)
        if len(intersect_s8) > 0:
            step8_variants = step7_variants.subtract(intersect_s8.sort(), A=True).sort().saveas("%s/step8.bed" % options.output_dir)
        else:
            step8_variants = step7_variants.saveas("%s/step8.bed" % options.output_dir)
    else:
        exclude = pybedtools.BedTool(options.exclude_bed).sort()
        intersect_s8 = step7_variants.window(exclude, w=1)
        if len(intersect_s8) > 0:
            step8_variants = step7_variants.subtract(intersect_s8.sort(), A=True).sort()
        else:
            step8_variants = step7_variants

    # step 9
    # remove variants where in a cluster within window of 20nt more than 75% of samples have AVF present
    if DEBUG is True:
        step9_variants = remove_avf_abundant_variants(step8_variants, window=cluster_window, occurrence_threshold=population_threshold).sort().saveas("%s/step9.bed" % options.output_dir)
    else:
        step9_variants = remove_avf_abundant_variants(step8_variants, window=cluster_window, occurrence_threshold=population_threshold).sort()

    # step 10
    # Remove synonymous unless in a 20nt cluster within the same sample
    clusters = step9_variants.cluster(d=cluster_window).sort()
    cluster_dict, id_to_cluster = clusters_to_list(clusters)
    variants_to_remove = []
    for cluster_id in cluster_dict:
        cluster = cluster_dict[cluster_id]
        if len(cluster) == 1 and re.search(r'\bsynonymous\b', cluster[0][all_var_cols["coding_matrix"]]):
            variants_to_remove.append("\t".join(cluster[0]))
        elif len(cluster) >= 2:
            # remove only if neighboring variants have same sample present
            for variant in cluster:
                 if re.search(r'\bsynonymous\b', variant[all_var_cols["coding_matrix"]]):
                     flag = None
                     if int(variant[all_var_cols["avf_above_zero"]]) == 1:
                         variant_avf_index = int(variant[all_var_cols["min_avf_index"]])
                     else:
                         variant_avf, variant_avf_index = get_max_from_numerical_string(variant[all_var_cols["avf_matrix"]])
                     for neighbor in cluster:
                         if neighbor != variant:
                             neighbor_avf_value = get_index_value_from_sample(neighbor[all_var_cols["avf_matrix"]], variant_avf_index)
                             if neighbor_avf_value >= min_avf:
                                 flag = 1
                                 break
                     if flag is None:
                         variants_to_remove.append("\t".join(variant))

    fd, path = tempfile.mkstemp()
    tmp_bed = open(path, "w")
    tmp_bed.write("\n".join(variants_to_remove))
    tmp_bed.close()
    if DEBUG is True:
        remove = pybedtools.BedTool(path).sort().saveas("%s/step10_remove.bed" % options.output_dir)
        step10_variants = subtract_variants(step9_variants,remove).sort().saveas("%s/step10.bed" % options.output_dir)
    else:
        remove = pybedtools.BedTool(path).sort()
        step10_variants = subtract_variants(step9_variants,remove).sort()

    # step 11
    # Remove Exac (BP) or 1000G (BV) if MAF is >0.002 unless in a 20nt cluster
    clusters = step10_variants.cluster(d=cluster_window).sort()
    cluster_dict, id_to_cluster = clusters_to_list(clusters)
    variants_to_remove = []
    maf_threshold = 0.002
    for cluster_id in cluster_dict:
        cluster = cluster_dict[cluster_id]
        if len(cluster) == 1 and ( float(cluster[0][all_var_cols["thousandG_matrix"]]) >= maf_threshold or float(cluster[0][all_var_cols["exac_matrix"]]) >= maf_threshold ):
            variants_to_remove.append("\t".join(cluster[0]))
    if len(variants_to_remove) > 0:
        fd, path = tempfile.mkstemp()
        tmp_bed = open(path, "w")
        tmp_bed.write("\n".join(variants_to_remove))
        tmp_bed.close()
        if DEBUG is True:
            remove = pybedtools.BedTool(path).sort().saveas("%s/step11_remove.bed" % options.output_dir)
            step11_variants = subtract_variants(step10_variants, remove).sort().saveas("%s/step11.bed" % options.output_dir)
        else:
            remove = pybedtools.BedTool(path).sort()
            step11_variants = subtract_variants(step10_variants, remove).sort()
    else:
        if DEBUG is True:
            step11_variants = step10_variants.saveas("%s/step11.bed" % options.output_dir)
        else:
            step11_variants = step10_variants

    # step 12
    # repeat step 7
    variants_to_remove = subset_by_rdf_avf(step11_variants,conditions, cluster_window)
    if len(variants_to_remove) > 0:
        if DEBUG is True:
            step12_variants = subtract_variants(step11_variants, variants_to_remove.sort()).sort().saveas("%s/step12.bed" % options.output_dir)
        else:
            step12_variants = subtract_variants(step11_variants, variants_to_remove.sort()).sort()
    else:
        if DEBUG is True:
            step12_variants = step11_variants.sort().saveas("%s/step12.bed" % options.output_dir)
        else:
            step12_variants = step11_variants.sort()

    final_variants = step12_variants
    _create_final_output(final_variants, options.output_bed, sample_list, stored_values, header)

    # write final vcfs if an input vcf directory is given
    _create_final_vcf_output(sample_list, final_variants, options.output_dir, stored_values, options.vcf_dir)
    

if __name__ == "__main__":
    __main__()
