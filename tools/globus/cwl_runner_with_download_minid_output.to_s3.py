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
from botocore.client import Config
import requests.packages.urllib3
requests.packages.urllib3.disable_warnings()
sys.path.append("/mnt/galaxyTools/tools/pymodules/python2.7/lib/python/globus_sdk-1.5.0-py2.7.egg")
import globus_sdk
from datetime import datetime
import hashlib
import uuid 
from bdbag import bdbag_ro as ro
import getpass

from identifier_client.identifier_api import IdentifierClient
from minid_client import minid_client_api

from pprint import pprint

from optparse import OptionParser



TMP_DIR = "/ephemeral/0/tmp/"


def file_as_bytes(file):
    with file:
        return file.read()


args=None
parser = OptionParser()

parser.add_option("--cwl-file", dest="cwl_file", help="CWL File")
parser.add_option("--inputs", dest="inputs", help="Inputs JSON")
parser.add_option("--output", dest="output", help="Output")
parser.add_option("--output-dir", dest="output_dir", help="Output_dir")
parser.add_option( '--history', dest="history_id", help="the history api id" )
parser.add_option("--stdout", dest="output_stdout", help="Output log stdout")
parser.add_option("--stderr", dest="output_stderr", help="Output log stderr")
parser.add_option( '-t', '--token-auth', dest="goauth_token", help='Globus auth token' )
parser.add_option( '-s', '--token-service', dest="goauth_service", help='Globus auth token' )

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
    os.mkdir("%s/data" % output_path)
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

    #s3 = boto3.resource('s3', config=Config(signature_version='s3v4'))
    #s3.meta.client.download_file(bucket, file_path, local_file_path)

    # Downloading from a public S3 bucket
    #s3 = boto3.resource('s3')
    #bucket = s3.Bucket(bucket)
    #bucket.download_file(file_path, local_file_path)

    # Downloading from NC
    session = boto3.session.Session(profile_name='ncs3conn')
    client = session.client('sts')
    rolearn = 'arn:aws:iam::600168050588:role/developer_access_gtex'
    #rolearn = 'arn:aws:iam::600168050588:role/developer_access_topmed'
    assumeRoleObject = response = client.assume_role(RoleArn=rolearn, RoleSessionName ='NIH-Test', DurationSeconds=3600 )
    credentials = assumeRoleObject['Credentials']
    s3 = boto3.client('s3',aws_access_key_id = credentials['AccessKeyId'],
            aws_secret_access_key = credentials['SecretAccessKey'],
            aws_session_token = credentials['SessionToken'])
    s3.download_file(Bucket=bucket, Key=file_path, Filename=local_file_path)

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
#copyfile(source_path, out_final)

#s3 transfer


#if os.path.exists(output):
#    os.remove(output)

local_ep = "f78c134c-9295-11e8-9698-0a6d4e044368"
local_root_path = "galaxy/files/output_bdbags_data"
#S3_BUCKET = 'gg-commons-tmp'
S3_BUCKET = 'gg-commons'

results_url = "https://results.fair-research.org"
results_path = "/results"
bdbags_url = "https://bags.fair-research.org"
bdbags_path = "/bags"

identifier = str(int(time.time()*10000))
s3_dir_to_use = "{0}_{1}".format(os.path.basename(outdir).strip()+"_1", identifier)
out_final = os.path.join(output_dir, output_file)
files_to_bdbag = []
# [path, md5, length, filename]
#os.makedirs("%s/%s" % (results_path, s3_dir_to_use))
#os.mkdir("%s/%s" % (results_path, s3_dir_to_use), 0755)
to_loc = "{0}/{1}".format(s3_dir_to_use, os.path.basename(source_path))
to_path = "%s/%s" % (results_path, to_loc)
real_url_path = "%s/%s" % (results_url, to_loc)
#gc_url_path = "globus://%s/%s/%s" % (local_ep, local_root_path, to_loc)
file_size = os.path.getsize(source_path)
md5sum = hashlib.md5(file_as_bytes(open(source_path, 'rb'))).hexdigest()
filename = os.path.basename(source_path)

#s3 transfer

#s3 = boto3.resource('s3')
session = boto3.session.Session(profile_name='gs3conn')
s3_client = session.client("s3", config=Config(signature_version="s3v4"))
#bucket = s3.Bucket(S3_BUCKET)
url_path = "https://s3.amazonaws.com/%s/results/%s" % (S3_BUCKET,  to_loc)
#bucket.upload_file(source_path, "results/%s" % (to_loc), ExtraArgs={'ACL':'public-read'})
#bucket.upload_file(source_path, "results/%s" % (to_loc))
s3_client.upload_file(source_path, S3_BUCKET, "results/%s" % (to_loc))

remote_manifest = []
remote_manifest.append({'url' : real_url_path, "length": file_size, "filename" : filename, "md5" : md5sum})
#remote_gc_manifest.append({'url' : gc_url_path, "length": file_size, "filename" : filename, "md5" : md5sum})

#files_to_bdbag.append([{'source' : source_path, 'md5' : md5sum, 'length' : file_size, 'filename' : os.path.basename(source_path), 'url_path' : url_path, 'gc_url_path' : gc_url_path, 'to_path' : to_path])

#print remote_manifest
tmp_file_path = os.path.join(TMP_DIR, "dbbag_history_info_{0}".format(identifier))
#print tmp_file_path
with open(tmp_file_path, 'w') as fp:
    json.dump(remote_manifest, fp)
bdbag_name = "%s.%s" % (s3_dir_to_use, "outputs.bdbag")
#bdbag_path = os.path.join(bdbags_path, bdbag_name)
bdbag_path = os.path.join(TMP_DIR, bdbag_name)
if not os.path.isdir(bdbag_path):
    os.mkdir(bdbag_path, 0755)
bdbag = bdbag_api.make_bag(bdbag_path, remote_file_manifest=tmp_file_path)

bdbag_fetch = "%s/fetch.txt" % bdbag_path

ro_metadata = dict()
ro_author_name = "Globus Genomics"
ro_author_orcid = "https://orcid.org/0000-0003-0723-665X"
ro_manifest = ro.init_ro_manifest(author_name=ro_author_name, author_orcid=ro_author_orcid)
ro_annotation_about = list()
ro_annotation_filename = "outputs.metadata.json"
ro_annotation_content_path = "annotations/%s" % ro_annotation_filename

fh = open(bdbag_fetch, "r")
for line in fh:
    line = line.rstrip("\n")
    values = line.split("\t")
    url = values[0]
    filename = values[2]
    if filename.startswith("data/"):
        filename = filename.replace("data/", "")
    ro_annotation_about.append(ro.ensure_payload_path_prefix(filename))
    ro.add_file_metadata(ro_manifest,
                         source_url=url, # url parameter from the RFM
                         bundled_as=ro.make_bundled_as(
                                    folder=os.path.dirname(filename), # filename parameter from the RFM
                                    filename=os.path.basename(filename))) # filename parameter from the RFM
#ro.add_annotation(ro_manifest, ro_annotation_about, uri="urn:uuid:%s" % str(uuid.uuid4()), content=ro_annotation_content_path)
ro_metadata["manifest.json"] = ro_manifest

#ro_metadata[ro_annotation_content_path] = history_info_to_write # ALEX: this is the metadata from Galaxy in dict form

# save the BDBAG with the RO data and archive
ro.serialize_bag_ro_metadata(ro_metadata, bdbag_path)
bdbag.save()
# archive the bag
archived_path = bdbag_api.archive_bag(bdbag_path, "zip")

# upload output bdbag to s3
remote_path = "bags/%s" % ( os.path.basename(archived_path))
#bucket.upload_file(archived_path, remote_path, ExtraArgs={'ACL':'public-read'})
#bucket.upload_file(archived_path, remote_path)
#s3_client.upload_file(source_path, S3_BUCKET, remote_path)
s3_client.upload_file(archived_path, S3_BUCKET, remote_path)

s3_url = "%s/%s" % (bdbags_url, os.path.basename(archived_path))

# get user email and name from auth token
#globus_auth_token = '%s \'' % options.goauth_token
service_token = options.goauth_token

minid_title = "TOPMED workflow output results for %s" % os.path.basename(source_path)

ic = IdentifierClient('Identifier', base_url='https://identifiers.globus.org/', app_name='WES Service', authorizer=globus_sdk.AccessTokenAuthorizer(service_token))
#checksum = minid_client_api.compute_checksum(archived_path)
md5sum = hashlib.md5(file_as_bytes(open(archived_path, 'rb'))).hexdigest()
# HHxPIZaVDh9u - namespace for test minids
# kHAAfCby2zdn - namespace for production once ready
kwargs = {
        'namespace': 'kHAAfCby2zdn',
        'visible_to': json.dumps(['public']),
        'location': json.dumps([s3_url]),
        'checksums': json.dumps([{
              'function': 'md5',
              'value': md5sum
          }]),
        'metadata': json.dumps({
                'Title': minid_title
        })
}
response_with_minid = ic.create_identifier(**kwargs)
minid = response_with_minid['identifier']

out_dict = { 'minid' : minid, 'BDbag_files' : remote_manifest, 'BDbag_path' : s3_url, 'title': minid_title}
fh_out = open(output, "w")
fh_out.write( json.dumps(out_dict))
fh_out.close()
#copyfile(source_path, out_final)

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
