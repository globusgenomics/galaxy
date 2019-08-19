#!/usr/bin/python
import sys, os, glob
from bioblend.galaxy import GalaxyInstance
import requests
requests.packages.urllib3.disable_warnings()

output_dir = sys.argv[1]
output_file = sys.argv[2]
output_name = sys.argv[3]

# get history name from history id
gi = GalaxyInstance(url=sys.argv[5], key=sys.argv[6])
history = gi.histories.show_history(sys.argv[4])

fw = open("%s" % (output_file), "w")
fw.write("### METADATA\n#######################################\nWorkflow Name\tSNAP_BAG_MASTER_v.1.0.0\n")
fw.write("Workflow id\t04685446e5198eef\n")
fw.write("Project Name\t%s-%s\n#######################################\n\n" % (history['name'], "MASTER"))
fw.write("###TABLE DATA\n#######################################\n")
fw.write("SampleName\t0##SourceType::SourceName::BATCH_FILE\t1##Param::1::batch_submit_trena::sample_id\n")

contents = glob.glob("%s/*.txt" % output_dir)
files = []
for f in contents:
    base = os.path.basename(f)
    files.append(os.path.basename(f))
    fw.write("%s\thistory::%s::%s\t%s\n" % (base.split(".")[0], history['name'], output_name, base))
    fw.close()
