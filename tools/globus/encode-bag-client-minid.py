import requests
import sys, os
from bdbag import bdbag_api
import urllib

BAG_SERVICE = "http://encode.bdbag.org/encode"
QUERY_BASE = "http://minid.bd2k.org/minid/landingpage/"
minid = sys.argv[1]
BASE_DOWNLOAD_PATH = sys.argv[2]

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
output_path = "%s/%s" %(extract_path, filename.split(".")[0])

print "Extracting bag and resolving fetch: %s" % output_path 
bdbag_api.extract_bag(path, extract_path)
bdbag_api.resolve_fetch(output_path, True)

print output_path
