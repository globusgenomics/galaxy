import requests
import sys, os
from bdbag import bdbag_api
import time
import urllib
import uuid
import shutil

QUERY_BASE = "http://minid.bd2k.org/minid/landingpage/"
minid = sys.argv[1]
USER_DOWNLOAD_PATH = sys.argv[2]
user_id = sys.argv[3].split("@")[0]
token = sys.argv[4]
#local_endpoint = "d1c150fa-b43d-11e7-b0a7-22000a92523b"
local_endpoint = "b475f4de-0536-11e8-a6a1-0a448319c2f8"

random_name = uuid.uuid4()
BASE_DOWNLOAD_PATH = "/scratch/shared/%s/%s" % (user_id, random_name)
print BASE_DOWNLOAD_PATH
# build the keychain.json file

kc_obj = """
[
    {
        "uri":"globus://",
        "auth_type": "token",
        "auth_params": {
            "local_endpoint": "%s",
            "transfer_token": "%s"
        }
    }
]
""" % (local_endpoint, token)
kf_name = "keychain.json"
kf = open(kf_name, "w")
kf.write(kc_obj)
kf.close()

if not os.path.exists(BASE_DOWNLOAD_PATH):
    os.makedirs(BASE_DOWNLOAD_PATH)

query = "%s%s" % (QUERY_BASE, minid)

print "Executing query: %s" % query

r = requests.get(query,  headers = {"Accept" : "application/json"})
print r

location = r.json()["locations"][0]['link']

filename = location.split("/")[-1]
path = "%s/%s" % (BASE_DOWNLOAD_PATH, filename)

print "Downloading result: %s" % location

testfile = urllib.URLopener()
testfile.retrieve(location, path)

# BDBag tooling doesnt let us unzip into a single dir.. so 
# we use a /base/uuid/uuid 
extract_path = ".".join(path.split(".")[0:-1])
#output_path = "%s/%s" %(extract_path, filename.split(".")[0])
output_path = "%s/%s" %(extract_path, ".".join(filename.split(".")[0:-1]))

print "Extracting bag and resolving fetch: %s" % output_path
bdbag_api.extract_bag(path, extract_path)
time.sleep(10)
bdbag_api.resolve_fetch(output_path, True, keychain_file=kf_name)

print output_path

# check that transfer is complete
# get fetch file for file size
fetch_file = "%s/fetch.txt" % output_path
print fetch_file
ff = open(fetch_file, "r")
sizes = []
for line in ff:
    values = line.rstrip("\n").split("\t")
    sizes.append({"name" : values[2], "length" : values[1]})

flag = True
while flag:
    for size in sizes:
        file_path = "%s/%s" % (output_path, size['name'])
        if os.path.exists(file_path) and os.path.getsize(file_path) == int(size['length']):
            flag = False
            continue
        else:
            flag = True
            break

# move BDBag with data to the USER_DOWNLOAD_PATH
print "move dir from %s to %s" % (extract_path, USER_DOWNLOAD_PATH)
shutil.move(extract_path, USER_DOWNLOAD_PATH)
os.remove(kf_name)
# Delete input download
shutil.rmtree(BASE_DOWNLOAD_PATH)
