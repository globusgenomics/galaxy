import requests
import sys
from bdbag import bdbag_api
import urllib

BAG_SERVICE = "http://encode.bdbag.org/encode"
#QUERY_BASE = "https://www.encodeproject.org/search/?type=Experiment&assay_slims=DNA+accessibility&assay_title=DNase-seq&"
QUERY_BASE = "https://www.encodeproject.org/search/?type=Experiment&assay_term_name=DNase-seq&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&biosample_term_name="
BASE_DOWNLOAD_PATH = sys.argv[2]

query = "%s%s" % (QUERY_BASE, sys.argv[1])

print "Executing query: %s" % query

r = requests.post(BAG_SERVICE, json = {'q': query}, headers = {"Content-Type" : "application/json", "Accepts" : "application/json"})

print r

url = r.json()["uri"]

filename = url.split("/")[-1]
path = "%s%s" % (BASE_DOWNLOAD_PATH, filename)

print "Downloading result: %s" % url

testfile = urllib.URLopener()
testfile.retrieve(url, path)

# BDBag tooling doesnt let us unzip into a single dir.. so 
# we use a /base/uuid/uuid 
extract_path = path.split(".")[0]
output_path = "%s/%s" %(extract_path, filename.split(".")[0])

print "Extracting bag and resolving fetch: %s" % output_path 
bdbag_api.extract_bag(path, extract_path)
bdbag_api.resolve_fetch(output_path, True)

print output_path
