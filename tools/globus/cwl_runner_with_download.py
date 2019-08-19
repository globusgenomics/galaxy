#!/usr/bin/python
import sys
import requests
import os
from bdbag import bdbag_api
import urllib
import time
import json
from shutil import copyfile, rmtree
import subprocess
import boto3
import requests.packages.urllib3
requests.packages.urllib3.disable_warnings()

from pprint import pprint

from optparse import OptionParser


TMP_DIR = "/ephemeral/0/tmp/"

args=None
parser = OptionParser()

parser.add_option("--cwl-file", dest="cwl_file", help="CWL File")
parser.add_option("--inputs", dest="inputs", help="Inputs JSON")
parser.add_option("--output", dest="output", help="Output")
parser.add_option("--output-dir", dest="output_dir", help="Output_dir")
parser.add_option("--stdout", dest="output_stdout", help="Output log stdout")
parser.add_option("--stderr", dest="output_stderr", help="Output log stderr")
options, args = parser.parse_args(args)

#cwl_file = options.cwl_file
#cwl_file = "/home/galaxy/sbg_dockstore_tools/topmed-workflows/alignment/topmed-alignment.cwl"
cwl_file = "/home/galaxy/topmed-workflows/aligner/sbg-alignment-cwl/topmed-alignment.cwl"
inputs_tmp = options.inputs
json_acceptable_string = inputs_tmp.replace("'", "\"")
inputs = json.loads(json_acceptable_string)
output = options.output
output_dir = options.output_dir

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# work on the cwl-worflow


# work on the json inputs file

def decide_path_type(path):
    if path.startswith("ark:"):
        return "minid"
    elif path.startswith("s3:"):
        return "s3"
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
    file_list = os.listdir(tmp_path)

    def find_cram_file(f_l):
        for item in f_l:
            if item.endswith('.cram'):
                return item
        return None

    if len(file_list) > 1:
        file_name = find_cram_file(file_list)
        if file_name == None:
            file_name = file_list[0]
    else:
        file_name = file_list[0]

    file_path = os.path.join(tmp_path, file_name)
    return file_path

def download_file_from_s3(s3_path, local_dir):
    tmp1 = s3_path.split('//')
    tmp2 = tmp1[1].split('/', 1)
    bucket = tmp2[0]
    file_path = tmp2[1]
    file_name = tmp1[1].split('/')[-1]

    local_dir = os.path.join(local_dir, str(int(time.time()*10000)))

    if not os.path.exists(local_dir):
        os.makedirs(local_dir)

    local_file_path = os.path.join(local_dir, file_name)

    s3 = boto3.resource('s3')
    bucket = s3.Bucket(bucket)
    bucket.download_file(file_path, local_file_path)

    return local_file_path
    
def translate_path(path):
    path_type = decide_path_type(path)
    # minid
    if path_type == "minid":
        file_path = download_file_from_minid(path, TMP_DIR)
        return file_path
    elif path_type == "s3":
        file_path = download_file_from_s3(path, TMP_DIR)
        return file_path

if "workflow_params" in inputs:
    inputs = inputs["workflow_params"]
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

json_inputs_dir = os.path.join(TMP_DIR, str(int(time.time()*10000)))
if not os.path.exists(json_inputs_dir):
    os.makedirs(json_inputs_dir)
json_inputs_file = os.path.join(json_inputs_dir, 'inputs.json')

with open(json_inputs_file, "w") as write_file:
    json.dump(inputs_to_write, write_file)

# run the workflow
outdir = os.path.join(TMP_DIR, str(int(time.time()*10000)))
if not os.path.exists(outdir):
    os.makedirs(outdir)

command = 'cwl-runner --no-match-user --rm-tmpdir  --outdir {0} --tmpdir-prefix {3} --tmp-outdir-prefix {3} {1} {2} > {4} 2> {5}'.format(outdir, cwl_file, json_inputs_file, TMP_DIR, options.output_stdout, options.output_stderr)
print "CWL Command: ", command

subprocess.call(command, shell=True)


output_list = os.listdir(outdir)
#print output_list
if output_list == []:
    sys.exit("No output file.")
else:
    output_file = output_list[0]

source_path = os.path.join(outdir, output_file)

#if os.path.exists(output):
#    os.remove(output)

out_final = os.path.join(output_dir, output_file)
copyfile(source_path, out_final)
#os.system("ls -l %s > %s" % (output_dir, output))


# delete all content in the tmp dir
for the_file in os.listdir(TMP_DIR):
    file_path = os.path.join(TMP_DIR, the_file)
    try:
        if os.path.isfile(file_path):
            os.remove(file_path)
        elif os.path.isdir(file_path):
            rmtree(file_path)
    except Exception as e:
        print(e)

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
