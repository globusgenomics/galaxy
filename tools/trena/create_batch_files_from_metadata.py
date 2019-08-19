#!/usr/bin/python
import sys, os

metadata_file = sys.argv[1]
fh = open(metadata_file, "r")
header = fh.readline()
#print header
groups = {}
files = {}
paired = {}
for line in fh:
    values = line.split("\t")
    acc_id = values[3]
    files[values[0]] = values[39]
    if values[32] == "paired-ended":
        paired[values[0]] = [values[33], values[34]]
    else:
        paired[values[0]] = None
    if acc_id in groups:
        groups[acc_id].append(values[0])
    else:
        groups[acc_id] = [values[0]]

for acc_id in groups:
    fw = open("%s/%s.batch_submit.txt" % (sys.argv[2], acc_id), "w")
    fw.write("### METADATA\n#######################################\nWorkflow Name\tsnap_from_bag\n")
    fw.write("Workflow id\t27ff0cc6ccf7f82f\n")
    fw.write("Project Name\t%s\n#######################################\n\n" % acc_id)
    fw.write("###TABLE DATA\n#######################################\n")
    fw.write("SampleName\t##Param::3339::get_from_encode_url::input\n")

    for file_id in groups[acc_id]:
        if paired[file_id] is not None and paired[file_id][0] == "1":
            fw.write("%s-%s\t%s,%s\n" % (acc_id, file_id, files[file_id], files[paired[file_id][1]]))
        elif paired[file_id] is None:
            fw.write("%s-%s\t%s\n" % (acc_id, file_id, files[file_id]))
    fw.close()
