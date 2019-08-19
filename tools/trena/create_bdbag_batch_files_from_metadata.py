#!/usr/bin/python
import sys, os
from bioblend.galaxy import GalaxyInstance
import requests
requests.packages.urllib3.disable_warnings()

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
    tissue = values[6]
    files[values[0]] = values[39]
    if values[32] == "paired-ended":
        paired[values[0]] = [values[33], values[34]]
    else:
        paired[values[0]] = None
    if acc_id in groups:
        groups[acc_id].append(values[0])
    else:
        groups[acc_id] = [values[0]]

# get history name from history id
gi = GalaxyInstance(url=sys.argv[5], key=sys.argv[6])
history = gi.histories.show_history(sys.argv[4])

for acc_id in groups:
    fw = open("%s/%s.batch_submit.txt" % (sys.argv[2], acc_id), "w")
    fw.write("### METADATA\n#######################################\nWorkflow Name\tsnap_from_minid_bag\n")
    fw.write("Workflow id\tc065e7bd4a1680ca\n")
    fw.write("Project Name\t%s-%s\n#######################################\n\n" % (tissue, acc_id))
    fw.write("###TABLE DATA\n#######################################\n")
    fw.write("SampleName\t0##SourceType::SourceName::BAG\t1##Param::1::SNAP_Alignment_id_bag::input_source::sample_id\n")

    for file_id in groups[acc_id]:
        fw.write("%s-%s\thistory::%s::%s\t%s\n" % (acc_id, file_id, history['name'], sys.argv[3], file_id))
    fw.close()
