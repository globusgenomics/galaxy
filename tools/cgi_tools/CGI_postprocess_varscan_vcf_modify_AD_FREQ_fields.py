#!/usr/bin/python
import sys, os

vcf_file = sys.argv[1]
sample_name = sys.argv[2]
variant_caller = sys.argv[3]

v_fh = open(vcf_file, "r")
for line in v_fh:
    line = line.rstrip("\n")
    if "Sample1" in line:
        line = line.replace("Sample1", sample_name)
        print line
    elif line.startswith("#") and "ID=FREQ" in line:
        print "##FORMAT=<ID=AVF,Number=1,Type=Float,Description=\"Variant allele frequency\">"
    elif line.startswith("#"):
        print line
    elif not line.startswith("#"):
        values = line.split("\t")

        # re-format the FORMAT column with AVF instead of FREQ
        format_keys_col = values[8]
        keys = format_keys_col.split(":")
        new_keys = []
        index_values = {}
        index = 0
        if variant_caller == "varscan":
            avf_key = "FREQ"
        elif variant_caller == "freebayes":
            avf_key = "NONE"
            #new_keys.append("AVF")
        elif variant_caller == "mutect2":
            avf_key = "AF"
 
        for i in keys:
            if i == avf_key:
                new_keys.append("AVF")
                index_values["AVF"] = index
            else:
                new_keys.append(i)
                index_values[i] = index
            index += 1
        if variant_caller == "freebayes":
            new_keys.append("AVF")
            index_values['AVF'] = index
            new_keys.append("AD")
            index_values['AD'] = index+1

        #if variant_caller == "mutect2":
        #    new_keys.append("AVF")
        #    index_values['AVF'] = index
        #    #new_keys.append("AD")
        #    #index_values['AD'] = index+1

        values[8] = ":".join(new_keys)
        if variant_caller == "varscan":
            values[8] = values[8].replace("DP:RD:AD", "DP:AD")
        #print values[8]
 
        # re-format the sample column with the modified AVF column and AD values
        if variant_caller == "varscan":
            format_values_col = values[9]
            key_values = format_values_col.split(":")
            key_values[index_values["AD"]] = "%s,%s" % (key_values[index_values["RD"]],key_values[index_values["AD"]])
            #print index_values["AD"]
            #print key_values[index_values["AD"]]
            del key_values[index_values["RD"]]
            key_values[index_values["AVF"]-1] = key_values[index_values["AVF"]-1].replace("%", "")
            key_values[index_values["AVF"]-1] = str(float(key_values[index_values["AVF"]-1])/100)
            #print index_values["AVF"]
            values[9] = ":".join(key_values)
            print "%s" % ("\t".join(values))
        elif variant_caller == "freebayes":
            format_values_col = values[9]
            key_values = format_values_col.split(":")
            #print key_values
            #print format_values_col
            AVF = str(float(key_values[index_values["AO"]]) / (float(key_values[index_values["DP"]])))
            AD = "%s,%s" % (key_values[index_values["RO"]], key_values[index_values["AO"]])
            key_values.append(AVF)
            key_values.append(AD)
            values[9] = ":".join(key_values)
            print "%s" % ("\t".join(values))
        elif variant_caller == "mutect2":
            format_values_col = values[9]
            key_values = format_values_col.split(":")

            print "%s" % ("\t".join(values))

