#!/usr/bin/python

import globus_sdk
from bioblend.galaxy import GalaxyInstance
import optparse
import sys
from datetime import datetime
import os
import time
import boto3
from pprint import pprint
import hashlib
import uuid, json
from bdbag import bdbag_api
from bdbag import bdbag_ro as ro
import shutil
import subprocess
import getpass
#print getpass.getuser()
sql_path = '/opt/galaxy/.venv/lib/python2.7/site-packages'
sys.path.append(sql_path)
import sqlalchemy as sa
from sqlalchemy import Table
from identifier_client.identifier_api import IdentifierClient
from minid_client import minid_client_api

def get_dataset_ids (obj, gi):
    datasets = []
    if obj['type'] == 'file':
        datasets.append(obj['id'])
    elif obj['type'] == 'collection':
        collection_info = gi.histories.show_dataset_collection(obj['history_id'], obj['id'])
        #print collection_info
        col_type = collection_info['collection_type']
        #print col_type
        if col_type == "list:paired":
            for element in collection_info['elements']:
                for element2 in element['object']['elements']:
                    #pprint(element)
                    datasets.append([element2['object']['id'], element2['element_identifier']])
        elif col_type == "list":
            #print "Obj: %s" % obj
            for element in collection_info['elements']:
                #pprint(element)
                datasets.append([element['object']['id'],element['element_identifier']])
    return datasets

def get_filepath(UUID):
    GALAXY_DATABASE_CONN = "postgresql://galaxy:globus_genomics_pass@rds.ops.globusgenomics.org:5432/galaxy_nihcommons"
    galaxy_db_conn = sa.create_engine(GALAXY_DATABASE_CONN).connect()
    galaxy_db_meta = sa.MetaData(bind=galaxy_db_conn)

    dataset_table = sa.Table("dataset", galaxy_db_meta, autoload=True)

    dataset_id = galaxy_db_conn.execute(sa.sql.select([dataset_table.c.id]).where(dataset_table.c.uuid == UUID)).fetchone()[0]
    dir_num =  get_dir_num(dataset_id)
    file_path = "/scratch/galaxy/files/{0}/dataset_{1}.dat".format(get_dir_num(dataset_id), dataset_id)
    #print file_path
    return file_path

def get_dir_num(dataset_id):
  if dataset_id < 1000:
    return "000"
  tmp = str(dataset_id)[:-3]
  if len(tmp) == 1:
    return "00{0}".format(tmp)
  else:
    return "0{0}".format(tmp)


def file_as_bytes(file):
    with file:
        return file.read()

os.environ['AWS_SHARED_CREDENTIALS_FILE'] = '/home/galaxy/.aws/config'
parser = optparse.OptionParser()
parser.add_option( '-k', '--api-key', dest='api_key', help='Mandatory. Need to get your galaxy API key by visiting the Galaxy webpage http://dev.globusgenomics.edu' )
parser.add_option( '-u', '--url', dest='url_name', help='Galaxy URL' )
parser.add_option( '--history', dest="history_id", help="the history api id" )
parser.add_option( '-t', '--token-auth', dest="goauth_token", help='Globus auth token' )
parser.add_option( '-s', '--token-service', dest="goauth_service", help='Globus auth token' )
(options, args) = parser.parse_args()


# skip these tool ids for bdbad 
SKIP_TOOL_IDS = ['get_bdbag_from_minid', 'create_output_bdbag_and_minid', 'create_output_bdbag_and_minid_from_collections', 'create_output_bdbag_and_minid2', 'cwl_prepare_inputs']
SKIP_COLLECTION_IDS = ["forward", "reverse", "Collections on ark:", "MarkDups_Dupes Marked.html", "STAR", "rsem_sample.isoform_abundances", "rsem_sample.rsem_log", "MarkDups_Dupes Marked.bam"]

TMP_DIR = '/scratch/galaxy/tmp'
if not os.path.isdir(TMP_DIR):
    os.mkdir(TMP_DIR, 0755)

S3_BUCKET = 'gg-commons-tmp'

history_info_to_write = {}
dataset_to_transfer = []

# check file name duplicate
file_name_pool = []
def set_file_name(name):
    name = name.replace(' ', '_').strip("/").strip("\\").strip()
    if name in file_name_pool:
        name_in_use = True
        num = 0
        new_name = None
        while name_in_use == True:
            new_name = name + "_{0}".format(num)
            if new_name not in file_name_pool:
                name_in_use = False
            else:
                num = num + 1
        name = new_name

    file_name_pool.append(name)
    return name

url = options.url_name
if "http:" in options.url_name:
   url = options.url_name.replace("http", "https")
if "dev1" in url:
   url = url.replace("dev1", "dev")

key = options.api_key
#key = "hhhhhhmpcvthfopjdthmpcvthfopjdt"
if len(key) == 0:
    sys.exit("missing api key")

gi = GalaxyInstance(url=url, key=key)

history_info = gi.histories.show_history(options.history_id, contents=True)
##print(history_info)
history_name = gi.histories.get_histories(history_id=options.history_id)[0]['name']
##print history_name
counter = 0
d_seen = []
for dataset in history_info:
    if dataset['id'] in d_seen:
        continue
    #d_in_object = get_dataset_ids(dataset, gi)
    if dataset['type'] == 'collection':
        #collection_info = gi.histories.show_dataset_collection(options.history_id, dataset['id'])
        #pprint(collection_info)
        d_in_object = get_dataset_ids(dataset, gi)
        #print "COLLECTION DATASETS: %s" % d_in_object

        for d,d_name in d_in_object:
            d_seen.append(d)
            #print "D: %s" % d
            dataset_info = gi.datasets.show_dataset(d)
            #pprint(dataset_info)
            job_info = gi.jobs.show_job(dataset_info['creating_job'])
            #pprint(job_info)
            duration = (datetime.strptime(job_info['update_time'], '%Y-%m-%dT%H:%M:%S.%f') - \
                datetime.strptime(job_info['create_time'], '%Y-%m-%dT%H:%M:%S.%f')).total_seconds()
            if job_info['tool_id'] not in SKIP_TOOL_IDS and \
                    dataset_info['state'] == 'ok'  and \
                    dataset_info['deleted'] is False:
                dataset_info_to_append = {
                    "job_name": dataset_info['name'],
                    "tool_id": job_info['tool_id'],
                    "job_create_time": job_info['create_time'],
                    "job_end_time": job_info['update_time'],
                    "duration_in_sec": duration
                }
                history_info_to_write[dataset_info['name']] = dataset_info_to_append
                flag_col = 0
                for ID in SKIP_COLLECTION_IDS:
                    if ID in dataset_info['name']:
                        flag_col = 1
                        break
                
                if dataset['visible'] is True and flag_col == 0:
                    #print dataset_info['name']
                    #pprint(collection_info)
                    name_to_use = d_name + "_" + set_file_name(dataset_info['name'])
                    #print "D: %s" % d
                    #print name_to_use
                    #pprint(dataset_info)
                    file_path = None
                    if 'file_name' in dataset_info:
                        file_path = dataset_info['file_name']
                    else:
                        UUID = dataset_info['uuid'].replace("-", "")
                        #print UUID
                        file_path = get_filepath(UUID)
                    dataset_to_transfer_to_append = {
                        "dataset_location": file_path,
                        "name_to_use": name_to_use
                    }
                    dataset_to_transfer.append(dataset_to_transfer_to_append)
                    # check if there are other directories associated with output file
                    dir_path = file_path.replace(".dat", "_files")
                    if os.path.exists(dir_path):
                        for root, dirs, files in os.walk(dir_path):
                            for file in files:
                                dataset_location = "%s/%s" % (root, file)
                                name_to_use = d_name + "/" + dataset_location.replace(dir_path + "/", "")
                                dataset_to_transfer_to_append = {
                                    "dataset_location": dataset_location,
                                    "name_to_use": name_to_use
                                }
                                dataset_to_transfer.append(dataset_to_transfer_to_append)
    else:
        d_seen.append(dataset['id'])
        dataset_info = gi.datasets.show_dataset(dataset['id'])
        job_info = gi.jobs.show_job(dataset_info['creating_job'])
        duration = (datetime.strptime(job_info['update_time'], '%Y-%m-%dT%H:%M:%S.%f') - \
            datetime.strptime(job_info['create_time'], '%Y-%m-%dT%H:%M:%S.%f')).total_seconds()
#        dataset_info_to_append = {
#            "job_name": dataset_info['name'],
#            "tool_id": job_info['tool_id'],
#            "command_line": job_info['command_line'],
#            "job_create_time": job_info['create_time'],
#            "job_end_time": job_info['update_time'],
#            "duration_in_sec": duration
#        }
        dataset_info_to_append = {
            "job_name": dataset_info['name'],
            "tool_id": job_info['tool_id'],
            "job_create_time": job_info['create_time'],
            "job_end_time": job_info['update_time'],
            "duration_in_sec": duration
        }
        history_info_to_write[dataset_info['name']] = dataset_info_to_append
        if job_info['tool_id'] not in SKIP_TOOL_IDS and \
            dataset_info['visible'] is True and \
            dataset_info['state'] == 'ok'  and \
            dataset_info['deleted'] is False:
            #print dataset_info
            #counter += 1
            name_to_use = set_file_name(dataset_info['name'])
            file_path = None
            if 'file_name' in dataset_info:
                file_path = dataset_info['file_name']
            else:
                UUID = dataset_info['uuid'].replace("-", "")
                file_path = get_filepath(UUID)
            dataset_to_transfer_to_append = {
                "dataset_location": file_path,
                "name_to_use": name_to_use
            }
            dataset_to_transfer.append(dataset_to_transfer_to_append)
            # check if there are other directories associated with output file
            dir_path = file_path.replace(".dat", "_files")
            if os.path.exists(dir_path):
                for root, dirs, files in os.walk(dir_path):
                    for file in files:
                        dataset_location = "%s/%s" % (root, file)
                        name_to_use = dataset_location.replace(dir_path + "/", "")
                        dataset_to_transfer_to_append = {
                            "dataset_location": dataset_location,
                            "name_to_use": name_to_use
                        }
                        dataset_to_transfer.append(dataset_to_transfer_to_append)
    
#pprint(history_info_to_write)
#print(dataset_to_transfer)
#sys.exit()
# write file
identifier = str(int(time.time()*10000))

tmp_file_path = os.path.join(TMP_DIR, "dbbag_history_info_{0}".format(identifier))
with open(tmp_file_path, "w") as write_file:
    json.dump(history_info_to_write, write_file)


# transfer to S3
s3_dir_to_use = "{0}_{1}".format(history_name.replace(' ', '_').strip("/").strip("\\").strip()+"_1", identifier)
s3 = boto3.resource('s3')
bucket = s3.Bucket(S3_BUCKET)

history_info_file_s3_path = "{0}/{1}".format(s3_dir_to_use, "history_info_{0}".format(identifier))
bucket.upload_file(tmp_file_path, history_info_file_s3_path, ExtraArgs={'ACL':'public-read'})

os.remove(tmp_file_path)
remote_manifest = []
remote_gc_manifest = []
#local_ep = "b475f4de-0536-11e8-a6a1-0a448319c2f8"
local_ep = "f78c134c-9295-11e8-9698-0a6d4e044368"
local_root_path = "galaxy/files/output_bdbags_data"

os.makedirs("%s/%s/%s" % ("/scratch", local_root_path, s3_dir_to_use))

for item in dataset_to_transfer:
    from_loc = item["dataset_location"]
    to_loc = "{0}/{1}".format(s3_dir_to_use, item["name_to_use"])
    #print "FROM: %s" % from_loc
    #print "TO: %s" % to_loc
    file_size = os.path.getsize(from_loc)
    md5sum = hashlib.md5(file_as_bytes(open(from_loc, 'rb'))).hexdigest()
    filename = os.path.basename(item["name_to_use"])
    to_path = "%s/%s/%s" % ("/scratch", local_root_path, to_loc)
    url_path = "https://s3.amazonaws.com/%s/%s" % (S3_BUCKET ,to_loc)
    gc_url_path = "globus://%s/%s/%s" % (local_ep, local_root_path, to_loc)
    bucket.upload_file(from_loc, to_loc, ExtraArgs={'ACL':'public-read'})
    if os.path.exists(os.path.dirname(to_path)):
        shutil.copyfile(from_loc,to_path)
    else:
        os.makedirs(os.path.dirname(to_path))
        shutil.copyfile(from_loc,to_path)
    remote_manifest.append({'url' : url_path, "length": file_size, "filename" : filename, "md5" : md5sum})
    remote_gc_manifest.append({'url' : gc_url_path, "length": file_size, "filename" : filename, "md5" : md5sum}) 

#print remote_manifest
tmp_file_path = os.path.join(TMP_DIR, "dbbag_history_info_{0}".format(identifier))
#print tmp_file_path
with open(tmp_file_path, 'w') as fp:
    json.dump(remote_gc_manifest, fp)
bdbag_path = os.path.join(TMP_DIR, "%s.%s" % (s3_dir_to_use, "outputs.bdbag"))
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
ro_metadata[ro_annotation_content_path] = history_info_to_write # ALEX: this is the metadata from Galaxy in dict form

# save the BDBAG with the RO data and archive
ro.serialize_bag_ro_metadata(ro_metadata, bdbag_path)
bdbag.save()
# archive the bag
archived_path = bdbag_api.archive_bag(bdbag_path, "zip")

# upload output bdbag to s3
remote_path = "bdbags/output_bags/%s" % os.path.basename(archived_path)
bucket.upload_file(archived_path, remote_path, ExtraArgs={'ACL':'public-read'})

url_path = "https://s3.amazonaws.com/%s" % S3_BUCKET
s3_url = "%s/%s" % (url_path, remote_path)
config_path = "/home/galaxy/.minid/minid-config.cfg"

# get user email and name from auth token
#globus_auth_token = '%s \'' % options.goauth_token
service_token = options.goauth_token
#cmd_globus = 'curl -X GET  https://auth.globus.org/v2/oauth2/userinfo  -H \'authorization: Bearer %s' %(globus_token)
##print cmd_globus
#p = subprocess.Popen(cmd_globus, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
#line =  p.stdout.readlines()[-1]
##print line
#info = json.loads(line)
#name = info['name']
#email = info['email']
minid_title = "TOPMED workflow output results for %s" % options.history_id

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

#cmd = "minid --config %s --register --location %s --title \"TOPMED RNA-seq output results for %s\" %s" % (config_path, s3_url, options.history_id, archived_path)
#cmd = "minid --register --location %s --globus_auth_token \"%s\" --name \"%s\" --email \"%s\" --title \"TOPMED RNA-seq output results for %s\" %s" % (s3_url, options.goauth, name, email, options.history_id, archived_path)
##print cmd
#p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
#minid = None
#for line in p.stdout.readlines():
#    #print line
#    if "Created/updated minid:" in line:
#        line = line.rstrip("\n")
#        minid = line.split(" ")[-1]
#retval = p.wait()
out_dict = { 'minid' : minid, 'BDbag_files' : remote_manifest, 'BDbag_path' : s3_url, 'title': minid_title}
print json.dumps(out_dict)
os.remove(tmp_file_path)
os.remove("%s.zip" % bdbag_path)
shutil.rmtree(bdbag_path)
