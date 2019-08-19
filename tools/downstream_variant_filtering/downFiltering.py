#!/usr/bin/python
import pybedtools, glob
import os, sys, tempfile, optparse, re
from numpy import sum
import time

# print out final format of output
def _create_final_output(final_variants, raw_file, final_file, sample_list, uniq_id_column):
    raw_fh = open(raw_file, "r")
    raw_dict = {}
    header = raw_fh.readline()
    for line in raw_fh:
        values = line.split("\t")
        raw_dict[values[uniq_id_column]] = line
    raw_fh.close()

    final_fh = open(final_file, "w")
    final_fh.write(header)
    for variant in final_variants:
        unique_id = variant[4]
        rdf_values = variant[13].split(",")
        values = raw_dict[unique_id].split("\t")
        final_fh.write(raw_dict[unique_id])
        final_fh.write("\t%s\t\t%s\t\t\t\t\t%s\t%s\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t%s\n" % (values[1], values[3], values[8], values[9],  "\t".join(rdf_values)))
    final_fh.close()

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
        vcf_fh = open(anno_file, "r")
        header = vcf_fh.readline()
        vcfout_fh = open("%s/%s.vcf" % (output_dir, sample_name), "w")
        vcfout_fh.write(_vcf_header(sample_name))
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

# write final filtered VCF files given the variants filtered by samples
def _create_final_vcf_output(sample_list, sample_variants, vcf_dir, output_dir, annotations_dir): 
    # get the necessary annotations from the raw annotated VCF files
    columns = ["Zygosity", "Count", "Frequency", "Coverage"]
    stored_values, sample_list = get_variant_values_from_raw_input(annotations_dir, columns)
    #print stored_values

    for vcf_file in glob.glob("%s/*.vcf" % vcf_dir):
        sample_name = os.path.basename(vcf_file)
        if sample_name.endswith(".vcf"):
            sample_id = sample_name.replace(".vcf", ".txt")
        vcf_fh = open(vcf_file, "r")
        vcfout_fh = open("%s/%s" % (output_dir, sample_name), "w")
        vcfout_fh.write(_vcf_header())
        info_field = "GT:AD:VF:DP"
        for line in vcf_fh:
            if not line.startswith("#"):
                values = line.split("\t")
                ref = values[3]
                allele = values[4]
                start = values[1]
                chrm = values[0]
                if len(ref) > 1:
                    start = str(int(start) + 1)
                    allele = "-"
                    ref = ref[1:]
                    unique_id = "%s_%s_%s_%s" % (chrm, start, ref, allele)
                    if unique_id in sample_variants[sample_id]:
                        info_content = _info_field(stored_values, unique_id, sample_id)
                        new_line = "\t".join(values[0:8]) + "\t" + "\t".join([info_field, info_content])
                        vcfout_fh.write(new_line + "\n")
                elif len(allele) > 1:
                    if "," not in allele:
                        ref = "-"
                        allele = allele[1:]
                        unique_id = "%s_%s_%s_%s" % (chrm, start, ref, allele)
                        if unique_id in sample_variants[sample_id]:
                            info_content = _info_field(stored_values, unique_id, sample_id)
                            new_line = "\t".join(values[0:8]) + "\t" + "\t".join([info_field, info_content])
                            vcfout_fh.write(new_line + "\n")
                    else:
                        ref = "-"
                        alleles = allele.split(",")
                        keep_alleles = []
                        index = 0
                        for allele_version in alleles:
                            allele_modified = allele_version[1:]
                            unique_id = "%s_%s_%s_%s" % (chrm, start, ref, allele_modified)
                            if unique_id in sample_variants[sample_id]:
                                keep_alleles.append([allele_version, index])
                            index += 1
                        # create new line in VCF for remaining alleles to keep
                        if len(keep_alleles) > 0:
                            allele_field = ",".join([x[0] for x in keep_alleles])
                            gt_field = values[9]
                            gt, info, dp = gt_field.split(":")
                            new_gt = '/'.join(str(x) for x in xrange(1,len(keep_alleles)+1))
                            info_values = info.split(",")
                            new_info = ",".join([info_values[x[1]] for x in keep_alleles]) + "," + info_values[-1]
                            new_gt_field = ":".join([new_gt, new_info, dp])
                            values[4] = allele_field
                            values[9] = new_gt_field
                            vcfout_fh.write("\t".join(values))
                else:
                    unique_id = "%s_%s_%s_%s" % (chrm, start, ref, allele)
                    if unique_id in sample_variants[sample_id]:
                        info_content = _info_field(stored_values, unique_id, sample_id)
                        new_line = "\t".join(values[0:8]) + "\t" + "\t".join([info_field, info_content])
                        vcfout_fh.write(new_line + "\n")
        vcfout_fh.close()
    return 0

# get each samples variants
def get_sample_variants(sample_list, final_variants):
    storage = {}
    for sample in sample_list:
        storage[sample] = {}

    for variant in final_variants:
        unique_id = variant[4]
        avf_values = variant[9].split(",")
        index = 0
        variant_dict = {}
        for avf in avf_values:
            sample_name = sample_list[index]
            if float(avf) > 0:
                storage[sample_name][unique_id] = float(avf)
            index += 1
    return storage

# get list of sample names as they appear in the raw annotated file given a list of sample names
def order_list_from_input(bed_file, sample_list):
    raw_fh = open(bed_file, "r")
    header = raw_fh.readline()
    header_values = header.split("\t")
    ordered_list = []
    for col in header_values:
        for sample in sample_list:
            if sample in col:
                ordered_list.append(sample)
    return ordered_list

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

# get the rdf values from a pure rdf matrix
def get_values_from_annotation_matrix(file_input, sample_list, uniq_id_index):
    fh_matrix = open(file_input, "r")
    matrix_values = {}
    header = fh_matrix.readline()
    sample_indexes = match_column_names(header, sample_list)
    for line in fh_matrix:
        values = line.split("\t")
        matrix_values[values[uniq_id_index]] = {}
        for i in sample_list:
            matrix_values[values[uniq_id_index]][i] = values[sample_indexes[i]]
    fh_matrix.close()
    return matrix_values

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
        matrix.append(float(anno_values[var_id].get(i,0)))
    return matrix

# return list of sample names from a file
# each line represents a sample name
def get_sample_list(sample_list_file):
    fh_file = open(sample_list_file, "r")
    sample_list = []
    for line in fh_file:
        line = line.rstrip("\n")
        sample_list.append(line)
    return sample_list

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
            max_sample_value = get_index_value_from_sample(neighbor[9], index)
            if max_value is None:
                max_value = [max_sample_value, index, neighbor]
            elif max_sample_value >= max_value[0]:
                max_value = [max_sample_value, index, neighbor]
    return max_value

def get_max_cluster_avf(cluster, feature):
    max_value = [0, 0]
    for neighbor in cluster:
        if feature != neighbor:
            max_sample_value, index = get_max_from_numerical_string(neighbor[9])
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
        id_to_cluster[feature.fields[4]] = feature.fields[-1]
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
        sample_count = len(cluster[0][9].split(","))
        occurrence_count = 0
        cluster_sample_values = []
        for feature in cluster:
            cluster_sample_values.append([float(i) for i in feature[9].split(",")])

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
    cluster_counts_bed = clusters.groupby(c=16, g=[16], ops=['count'], output=tmp_clusters_bed.name)
    cluster_counts = get_dict(tmp_clusters_bed.name)
    multiple_clusters = subset_clusterCount(clusters, (1), '<', cluster_counts, len(clusters[0].fields)-1).cut(range(15)).sort()
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
        unique_id = variant.fields[4]
        cluster_id = id_to_cluster[unique_id]
        avf_value = variant.fields[avf_threshold[1]]
        avf_sampleID = variant.fields[avf_threshold[1]-1]
        flag = None
        for feature in cluster_dict[cluster_id]:
            feature_unique_id = feature[4]
            if feature_unique_id != unique_id:
                neighbor_avf = get_index_value_from_sample(feature[9], int(avf_sampleID))
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
    #if feature[12] == 1 and feature[11] < 5:
    if float(feature[12]) == 1:
        sample_avf = float(feature[11])
        sample_index = int(feature[10])
        sample_rdf = get_index_value_from_sample(feature[13], sample_index)
    else:
        sample_avf, sample_index = get_min_above_threshold_from_numerical_string(feature[9], 5)
        sample_rdf = get_index_value_from_sample(feature[13], sample_index)
    if sample_rdf < float(max_rdf_threshold) and sample_rdf >= float(min_rdf_threshold) and sample_avf < float(percent_threshold):
        return True
    return False

def subset_variant_percent(bed, percent_threshold, min_rdf_threshold, max_rdf_threshold):
    result = bed.filter(thresholdvalue_filter, percent_threshold, min_rdf_threshold, max_rdf_threshold).saveas()
    return pybedtools.BedTool(result.fn)

# test a set of avf and rdf conditions
def _check_all_samples_satisfy_conditions(conditions, feature):
    sample_index_list = get_index_list_above_theshold(feature[9], 0)
    flag_list = []
    
    # test that for each sample it satisfies each condition
    # if a sample does not satisfy all of the conditions, return False
    # if all samples satisfy one of the conditions, return True
    for sample_index in sample_index_list:
        sample_avf = get_index_value_from_sample(feature[9], sample_index)
        sample_rdf = get_index_value_from_sample(feature[13], sample_index)

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
        if feature[4] != variant[4]:
            neighbor_values.append([get_index_value_from_sample(feature[9], index), get_index_value_from_sample(feature[13], index)])
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
                    sample_index_list = get_index_list_above_theshold(variant[9], 0)
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
    parser.add_option( '-g', '--gtf', dest="gtf_input", help='gtf' )
    parser.add_option( '-r', '--rdf', dest="rdf_input", help='read_depth_input' )
    parser.add_option( '-o', '--output', dest="output_bed", help="bed_output" )
    parser.add_option( '-e', '--exclude', dest="exclude_bed", help="exclude bed" )
    parser.add_option( '-w', '--window', dest="window", help="window size for clusters" )
    parser.add_option( '-t', '--population-threshold', dest="population_threshold", help="Population variant occurrence percentage threshold for a given cluster (0,1)" )
    parser.add_option( '-a', '--variant-annotations-directory', dest="variant_annotation_dir", help="if provided, the vcf annotations will be taken from here" )
    parser.add_option( '', '--vcf-directory', dest="vcf_dir", help="if provided, the raw vcf files will be printed out" )
    parser.add_option( '-s', '--sample-list', dest="sample_list", help="text file with exact name of sample name as it appears on input matrix annotated file" )
    (options, args) = parser.parse_args()

    if len(args) == 1:
        function = args[0]
    else:
        sys.exit("You will need to specify if you are converting (convert) from CLC raw to VCF or filtering variant output (filter)")

    if function == "convert":
        _convert_CLC_to_vcf(options.output_dir, options.variant_annotation_dir)
        return 0

    id_column = "Unique_Variants"
    #id_column = "Chromosome_Start_Reference_Allele"

    if not options.window:
        cluster_window = 20
    else:
        cluster_window = options.window
    if not options.population_threshold:
        population_threshold = 0.75
    else:
        population_threshold = options.population_threshold

    #conditions = [[1000, 0, 60], [20, 60, 100], [10, 100, 200], [5, 200, 100000]]
    conditions = [[1000, 0, 60], [10, 60, 200], [5, 200, 1000000]]
    # parse the gtf file
    ###gtf = pybedtools.BedTool(options.gtf_input)
    exons = None
    ###exons = subset_featuretypes(gtf, ('exon'), 2, operator="in").sort()

    # is variant annotation directory provided
    if options.variant_annotation_dir:
        columns = ["Count", "Coverage"]
        stored_values, sample_list = get_variant_values_from_raw_input(options.variant_annotation_dir, columns)
        rdf_values = stored_values["Coverage"]
        avf_values = stored_values["Count"]
        sample_list = order_list_from_input(options.vcf_input, sample_list)
    else:
        # parse read_depth_file
        sample_list = get_sample_list(options.sample_list)
        rdf_values = get_values_from_annotation_matrix(options.rdf_input, sample_list, 0)
      
        # parse avf values
        avf_values = get_values_from_annotation_matrix(options.vcf_input, sample_list, 2)

    # create a tmp BED file
    # Need to put the chr, start and end locations at the beggnining of the line
    vcf = open(options.vcf_input, "r")
    header = vcf.readline()
    header_index = get_header_index(header)
    vcf_tmp, vcf_tmp_filename = tempfile.mkstemp()
    vcf_fh = open(vcf_tmp_filename, "w")
    for line in vcf:
        values = line.split("\t")
        new_values = list()
        for value in values:
            if value=="":
                value = "."
            new_values.append(value)
        line = "\t".join(new_values)
        matrix = get_value_matrix(avf_values, sample_list, values[header_index[id_column]])
        rdf_matrix = get_value_matrix(rdf_values, sample_list, values[header_index[id_column]])
        if max(matrix) > 0:
            m = min(i for i in matrix if i > 0)
            count = sum(x > 0 for x in matrix)
            selected_sampleId = matrix.index(m)
            sample_name = sample_list[selected_sampleId]
            rdv = rdf_values[new_values[header_index[id_column]]][sample_name]
        else:
            m = 0.0
            count = 0
            selected_sampleId = "0"
            rdv = "0"
        if new_values[header_index["1000G_ALL"]] == ".":
            new_values[header_index["1000G_ALL"]] = "0"
        if new_values[header_index["ExAC_Freq"]] == ".":
            new_values[header_index["ExAC_Freq"]] = "0"

        vcf_fh.write("%s\n" % ("\t".join([new_values[header_index["Chromosome"]], new_values[header_index["Start"]], new_values[header_index["End"]], new_values[header_index["1000G_ALL"]], new_values[header_index[id_column]], new_values[header_index["ExAC_Freq"]], new_values[header_index["ExonicFuncrefgene"]], new_values[header_index["gene_name (Homo_sapiens_ensembl_v74_mRNA)"]], new_values[header_index["Funcrefgene"]], ",".join(map(str, matrix)), str(selected_sampleId), str(m), str(count), ",".join(map(str, rdf_matrix)), new_values[header_index["Non-synonymous"]] ])))
    vcf_fh.close()

    # create a BEDtools object
    #variants = pybedtools.BedTool(vcf_tmp_filename).sort().saveas('%s/all_variants.bed' % options.output_dir)
    variants = pybedtools.BedTool(vcf_tmp_filename).sort()

    # step 1.1
    # get exonic
    variants_exonic = subset_featuretypes(variants, ('exonic'), 8, operator="in" )
    # step 1.2
    # get splcing
    variants_splicing = subset_featuretypes(variants, ('splicing'), 8, operator="in" )
    # step 1.3
    # get non-annotated regions
    variants_none = subset_featuretypes(variants, ('.'), 8, operator="in" )
    # step 1.3.1
    # keep variants_none with non_synonymous yes, yes
    variants_none_yes = subset_featuretypes(variants_none, ('Yes, Yes'), 14, operator="in" )
    # step 1.3.2
    # keep variants_none with non_synonymous no, no
    variants_none_no = subset_featuretypes(variants_none, ('No, No'), 14, operator="in" )
    # step 1.3.3
    # keep variants_none with non_synonymous -, -
    variants_none_empty = subset_featuretypes(variants_none, ('-, -'), 14, operator="in" )

    # step 1.3.3.1
    # keep variants_none with non_synonymous -, - within +/- 2nt of an exon
    window = 2
    if len(variants_none_empty) > 0 and exons is not None:
        closest_exon_variants = variants_none_empty.closest(exons, io=True, d=True)
        variants_within_2nt_of_exon = subset_closest(closest_exon_variants, window)
    else:
        variants_within_2nt_of_exon = None
    variants_within_2nt_of_exon_final = None
    if variants_within_2nt_of_exon is not None:
        variants_within_2nt_of_exon_final = variants_within_2nt_of_exon.cut([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14])

    # step 2
    # merge 1.1, 1.2, 1.3.1, 1.3.2, 1.3.3.1
    step2_variants_annotated = variants_exonic.cat(variants_splicing, postmerge=False).sort()
    if len(variants_none_no) > 0 and len(variants_none_yes) > 0:
        step2_variants_none = variants_none_yes.cat(variants_none_no, postmerge=False).sort()
    elif len(variants_none_no) > 0:
        step2_variants_none = variants_none_no.sort()
    elif len(variants_none_yes) > 0:
        step2_variants_none = variants_none_yes.sort()
    else:
        step2_variants_none = []

    if len(step2_variants_none) > 0:
        step2_variants_annotated_none = step2_variants_annotated.cat(step2_variants_none, postmerge=False).sort()
    else:
        step2_variants_annotated_none = step2_variants_annotated.sort()

    if variants_within_2nt_of_exon_final is not None:
        #step2_variants = step2_variants_annotated_none.cat(variants_within_2nt_of_exon_final, postmerge=False).sort().saveas("%s/variants_step2.bed" % options.output_dir)
        step2_variants = step2_variants_annotated_none.cat(variants_within_2nt_of_exon_final, postmerge=False).sort()
    else:
        #step2_variants = step2_variants_annotated_none.saveas("%s/variants_step2.bed" % options.output_dir)
        step2_variants = step2_variants_annotated_none

    # step 3
    # remove variants where occurrence count is 1 and AVF is less than 5% and there are no variants within 20nt
    clusters = step2_variants.cluster(d=cluster_window).sort()
    tmp_clusters_bed = tempfile.NamedTemporaryFile()
    cluster_counts_bed = clusters.groupby(c=16, g=[16], ops=['count'], output=tmp_clusters_bed.name)
    cluster_counts = get_dict(tmp_clusters_bed.name)
    single_clusters = subset_clusterCount(clusters, (1), '==', cluster_counts, len(clusters[0].fields)-1).cut(range(15)).sort()
    step3_occurrence1 = subset_featuretypes(single_clusters, ('1'), 12, operator="==" ).sort()
    step3_occurrence1_avfunder5 = subset_featureMath(step3_occurrence1, (5), '<', 11).sort()
    #step3_variants = subtract_variants(step2_variants, step3_occurrence1_avfunder5).sort().saveas("%s/variants_step3.bed" % options.output_dir)
    step3_variants = subtract_variants(step2_variants, step3_occurrence1_avfunder5).sort()

    # step 4
    # remove variants where occurrence is 1 and AVF <5% and if for variants within 20nt it is also <5% AVF or 0
    occurrence_threshold = [1, 12] # tuple [occurence_threhsold, index_field_in_bed]
    avf_threshold = [5, 11] # tuple [avf_threhsold, index_field_in_bed]
    #step4_variants = remove_low_variance_in_window(step3_variants, occurrence_threshold, avf_threshold, cluster_window).sort().saveas("%s/variants_step4.bed" % options.output_dir)
    step4_variants = remove_low_variance_in_window(step3_variants, occurrence_threshold, avf_threshold, cluster_window).sort()

    # step 5
    # Remove if >1 sample with <5% AVF (but none above 5%) and none of the same samples have a second variant within 20nt >5% AVF
    clusters = step4_variants.cluster(d=cluster_window).sort()
    tmp_clusters_bed = tempfile.NamedTemporaryFile()
    cluster_counts_bed = clusters.groupby(c=16, g=[16], ops=['count'], output=tmp_clusters_bed.name)
    cluster_counts = get_dict(tmp_clusters_bed.name)
    cluster_dict, id_to_cluster = clusters_to_list(clusters)
    occurrence_gt1 = subset_featureMath(clusters, (1), '>', 12).cut(range(15)).sort()
    #occurrence_gt1_all_avfunder5 = subset_featuresall_under5(occurrence_gt1, (5), 9).sort().saveas('%s/occurrence_gt1_all_avfunder5.bed' % options.output_dir)
    occurrence_gt1_all_avfunder5 = subset_featuresall_under5(occurrence_gt1, (5), 9).sort()

    # search through low_occurrence_avf_bed for variants that are less than 5% AVF in matching sample cluster 
    variants_to_remove = []
    for variant in occurrence_gt1_all_avfunder5:
        unique_id = variant.fields[4]
        cluster_id = id_to_cluster[unique_id]
        sample_index_list = get_index_list_above_theshold(variant[9], 0)
        flag = None
        for feature in cluster_dict[cluster_id]:
            feature_unique_id = feature[4]
            if feature_unique_id != unique_id:
                for index in sample_index_list:
                    feature_avf_value = get_index_value_from_sample(feature[9], index) 
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
        remove = pybedtools.BedTool(path).sort()
        #step5_variants = subtract_variants(step4_variants, remove).sort().saveas("%s/variants_step5.bed" % options.output_dir)
        step5_variants = subtract_variants(step4_variants, remove).sort()
    else:
        #step5_variants = step4_variants.sort().saveas("%s/variants_step5.bed" % options.output_dir)
        step5_variants = step4_variants.sort()

    # step 6 and 7
    # If read depth is <60x, variants are filtered out.
    # each conditions is a triple ==> [percent_threshold, min_rdf_threshold, max_rdf_threshold]
    #variants_to_remove = subset_by_rdf_avf(step5_variants,conditions, cluster_window).sort().saveas("%s/remove_step7.bed" % options.output_dir)
    variants_to_remove = subset_by_rdf_avf(step5_variants,conditions, cluster_window).sort()
    if len(variants_to_remove) > 0:
        #step7_variants = subtract_variants(step5_variants, variants_to_remove.sort()).sort().saveas("%s/variants_step7.bed" % options.output_dir)
        step7_variants = subtract_variants(step5_variants, variants_to_remove.sort()).sort()
    else:
        #step7_variants = step5_variants.sort().saveas("%s/variants_step7.bed" % options.output_dir)
        step7_variants = step5_variants.sort()

    # step 8
    # Remove CTBP2, KDM6B, DCP1B, PDGFRA, KMT2C
    # exclude BED file
    # BEDTOOLS exclude object
    exclude = pybedtools.BedTool(options.exclude_bed).sort()
    #intersect_s8 = step7_variants.window(exclude, w=1).sort().saveas("%s/intersect_step8.bed" % options.output_dir)
    intersect_s8 = step7_variants.window(exclude, w=1).sort()
    #step8_variants = step7_variants.subtract(intersect_s8, A=True).sort().saveas("%s/variants_step8.bed" % options.output_dir)
    step8_variants = step7_variants.subtract(intersect_s8, A=True).sort()

    # step 9
    # remove variants where in a cluster within window of 20nt more than 75% of samples have AVF present
    #step9_variants = remove_avf_abundant_variants(step8_variants, window=cluster_window, occurrence_threshold=0.75).sort().saveas("%s/variants_step9.bed" % options.output_dir)
    step9_variants = remove_avf_abundant_variants(step8_variants, window=cluster_window, occurrence_threshold=population_threshold).sort()

    # step 10
    # Remove synonymous unless in a 20nt cluster within the same sample
    clusters = step9_variants.cluster(d=cluster_window).sort()
    cluster_dict, id_to_cluster = clusters_to_list(clusters)
    variants_to_remove = []
    for cluster_id in cluster_dict:
        cluster = cluster_dict[cluster_id]
        if len(cluster) == 1 and re.search(r'\bsynonymous SNV\b', cluster[0][6]):
            variants_to_remove.append("\t".join(cluster[0]))
        elif len(cluster) >= 2:
            # remove only if neighboring variants have same sample present
            for variant in cluster:
                 if re.search(r'\bsynonymous SNV\b', variant[6]):
                     flag = None
                     if int(variant[12]) == 1:
                         variant_avf_index = int(variant[10])
                     else:
                         variant_avf, variant_avf_index = get_max_from_numerical_string(variant[9])
                     for neighbor in cluster:
                         if neighbor != variant:
                             neighbor_avf_value = get_index_value_from_sample(neighbor[9], variant_avf_index)
                             if neighbor_avf_value >= 5:
                                 flag = 1
                                 break
                     if flag is None:
                         variants_to_remove.append("\t".join(variant))

    fd, path = tempfile.mkstemp()
    tmp_bed = open(path, "w")
    tmp_bed.write("\n".join(variants_to_remove))
    tmp_bed.close()
    #remove = pybedtools.BedTool(path).sort().saveas("%s/remove_step10.bed" % options.output_dir)
    remove = pybedtools.BedTool(path).sort()
    #step10_variants = subtract_variants(step9_variants,remove).sort().saveas("%s/variants_step10.bed" % options.output_dir)
    step10_variants = subtract_variants(step9_variants,remove).sort()

    # step 11
    # Remove Exac (BP) or 1000G (BV) if MAF is >0.002 unless in a 20nt cluster
    clusters = step10_variants.cluster(d=cluster_window).sort()
    cluster_dict, id_to_cluster = clusters_to_list(clusters)
    variants_to_remove = []
    maf_threshold = 0.002
    for cluster_id in cluster_dict:
        cluster = cluster_dict[cluster_id]
        if len(cluster) == 1 and ( float(cluster[0][3]) >= maf_threshold or float(cluster[0][5]) >= maf_threshold ):
            variants_to_remove.append("\t".join(cluster[0]))
    fd, path = tempfile.mkstemp()
    tmp_bed = open(path, "w")
    tmp_bed.write("\n".join(variants_to_remove))
    tmp_bed.close()
    remove = pybedtools.BedTool(path).sort()
    #step11_variants = subtract_variants(step10_variants, remove).sort().saveas("%s/variants_step11.bed" % options.output_dir)
    step11_variants = subtract_variants(step10_variants, remove).sort()

    # step 12
    # repeat step 7
    variants_to_remove = subset_by_rdf_avf(step11_variants,conditions, cluster_window)
    if len(variants_to_remove) > 0:
        #step12_variants = subtract_variants(step11_variants, variants_to_remove.sort()).sort().saveas("%s/variants_step12.bed" % options.output_dir)
        step12_variants = subtract_variants(step11_variants, variants_to_remove.sort()).sort()
    else:
        #step12_variants = step11_variants.sort().saveas("%s/variants_step12.bed" % options.output_dir)
        step12_variants = step11_variants.sort()

    final_variants = step12_variants
    _create_final_output(final_variants, options.vcf_input, options.output_bed, sample_list, header_index[id_column])

    # write final vcfs if an input vcf directory is given
    if options.vcf_dir:
        sample_variants = get_sample_variants(sample_list, final_variants)
        _create_final_vcf_output(sample_list, sample_variants, options.vcf_dir, options.output_dir, options.variant_annotation_dir)
    

if __name__ == "__main__":
    __main__()

