#!/usr/bin/python
import sys
import requests
import os
from bdbag import bdbag_api
import urllib
import time
import json
import requests.packages.urllib3
requests.packages.urllib3.disable_warnings()

from pprint import pprint

from optparse import OptionParser


TMP_DOWNLOAD_DIR = "/scratch/galaxy/tmp"


args=None
parser = OptionParser()

parser.add_option("--inputs", dest="inputs", help="Inputs JSON")
parser.add_option("--json-inputs-output", dest="json_inputs_output", help="JSON inputs output")

options, args = parser.parse_args(args)

inputs_tmp = options.inputs
json_acceptable_string = inputs_tmp.replace("'", "\"")
inputs = json.loads(json_acceptable_string)
json_inputs_output = options.json_inputs_output


# work on the json inputs file

def decide_path_type(path):
    if path.startswith("ark:"):
        return "minid"
    else:
        msg = "Path type not supported: {0}".format(path)
        sys.exit(msg)

def download_file_from_minid(minid, local_dir):
    QUERY_BASE = None
    if len(minid) > 17:
        QUERY_BASE = "https://identifiers.globus.org"
    else:
        QUERY_BASE = "http://minid.bd2k.org/minid/landingpage"
    minid = minid
    BASE_DOWNLOAD_PATH = os.path.join(local_dir, str(int(time.time()*10000)))

    if not os.path.exists(BASE_DOWNLOAD_PATH):
        os.makedirs(BASE_DOWNLOAD_PATH)

    query = "%s/%s" % (QUERY_BASE, minid)
    #print "Executing query: %s" % query

    r = requests.get(query,  headers = {"Accept" : "application/json"})
    #print r.json()

    location = None
    if len(minid) > 17:
        location = r.json()["location"][0]
    else:
        location = r.json()["locations"][0]['link']

    filename = location.split("/")[-1]
    path = "%s/%s" % (BASE_DOWNLOAD_PATH, filename)

    #print "Downloading result: %s" % location

    testfile = urllib.URLopener()
    testfile.retrieve(location, path)

    # BDBag tooling doesnt let us unzip into a single dir.. so
    # we use a /base/uuid/uuid
    extract_path = ".".join(path.split(".")[0:-1])
    output_path = "%s/%s" %(extract_path, filename.split(".")[0])

    #print "Extracting bag and resolving fetch: %s" % output_path
    bdbag_api.extract_bag(path, extract_path)
    bdbag_api.resolve_fetch(output_path, True)

    tmp_path = os.path.join(output_path, "data")
    file_name = os.listdir(tmp_path)[0]
    file_path = os.path.join(tmp_path, file_name)
    return file_path


def translate_path(path):
    path_type = decide_path_type(path)
    # minid
    if path_type == "minid":
        file_path = download_file_from_minid(path, TMP_DOWNLOAD_DIR)
        return file_path


inputs_to_write = inputs
if isinstance(inputs, dict):
    for k, v in inputs.iteritems():
        if isinstance(v, dict):
            for key, value in v.iteritems():
                if key == "path":
                    local_path = translate_path(value)
                    inputs_to_write[k][key] = local_path
else:
    sys.exit("inputs is not a dictionary.")

with open(json_inputs_output, "w") as write_file:
    json.dump(inputs_to_write, write_file)



"""
{
  "output_name": "output1",
  "output": {
    "path": "/home/ubuntu/test/output1",
    "class": "File"
  },
    "bwa_index": {
        "class": "File",
        "path": "ark:/57799/b91414"
    },
    "dbsnp": {
        "class": "File",
        "path": "ark:/57799/b9wb0k"
    },
    "input_file": {
        "class": "File",
        "path": "ark:/57799/b9rm50"
    },
    "reference_genome": {
        "class": "File",
        "path": "ark:/57799/b9mt4f"
    }
}

hs38DH.fa.tar    hs38DH_fa_tar_bdbag.zip    ark:/57799/b91414
Homo_sapiens_assembly38.dbsnp138.vcf.gz    Homo_sapiens_assembly38_dbsnp138_bdbag.zip    ark:/57799/b9wb0k
NWD176325.0005.recab.cram    NWD176325_0005_recab_cram_bdbag.zip    ark:/57799/b9rm50
hs38DH.fa    hs38DH_fa_bdbag.zip    ark:/57799/b9mt4f

{"bwa_index": {"class": "File", "path": "ark:/57799/b91414" }, "dbsnp": {"class": "File", "path": "ark:/57799/b9wb0k"}, "input_file": {"class": "File", "path": "ark:/57799/b9rm50"}, "reference_genome": {"class": "File", "path": "ark:/57799/b9mt4f"}}

"""