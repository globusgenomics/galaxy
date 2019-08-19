#!/usr/bin/python
#!/usr/bin/python
# -*- coding: utf-8 -*-

import time, optparse, os, shutil, subprocess, sys, glob, tempfile, psycopg2, socket, re
from bioblend.galaxy import GalaxyInstance
#from galaxy.web import url_for
import requests
requests.packages.urllib3.disable_warnings()
import json, gzip
import pandas as pd
from bdbag import bdbag_api
import urllib
from shutil import copyfile, rmtree
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
import fnmatch
import ConfigParser

def check_globus_token(token, config_file, globus_dir, username, token_type):
   token_type_to_server = {
       'transfer': 'transfer.api.globus.org',
       'auth': 'auth.globus.org',
       'identifier': 'identifiers.globus.org'
   }
   if token_type == 'identifier':
       Config = ConfigParser.ConfigParser()
       Config.read(config_file)
       client_id = Config.get('commons', 'client_id')
       client_secret = Config.get('commons', 'client_secret')
       client = globus_sdk.ConfidentialAppAuthClient(client_id, client_secret)
       token_file = os.path.join(globus_dir, 'tokens-identifiers', username)
   else:
       Config = ConfigParser.ConfigParser()
       Config.read(config_file)
       client_id = Config.get('globus', 'client_id')
       client_secret = Config.get('globus', 'client_secret')
       client = globus_sdk.ConfidentialAppAuthClient(client_id, client_secret)
       if token_type == 'transfer':
           token_file = os.path.join(globus_dir, 'tokens', username)
       elif token_type == 'auth':
           token_file = os.path.join(globus_dir, 'tokens-auth', username)
       else:
           sys.exit('token_type unsupported')
   if not client.oauth2_validate_token(token)['active']:
       refresh_token_file = token_file + '_refresh'
       if os.path.isfile(refresh_token_file):
           with open(refresh_token_file, 'r') as f:
               refresh_token = f.readline().strip()
           token_response = client.oauth2_refresh_token(refresh_token)
           token = token_response.by_resource_server[token_type_to_server[token_type]]['access_token']
           with open(token_file, 'w') as write_file:
               write_file.write(token)
   return token

def file_as_bytes(file):
    with file:
        return file.read()


def download_file_from_minid(minid_set, local_dir):
    minid,location = minid_set
    #QUERY_BASE = None
    #if len(minid) > 17:
    #    QUERY_BASE = "https://identifiers.globus.org"
    #else:
    #    QUERY_BASE = "http://minid.bd2k.org/minid/landingpage"
    #minid = minid
    BASE_DOWNLOAD_PATH = os.path.join(local_dir, str(int(time.time()*10000)))

    if not os.path.exists(BASE_DOWNLOAD_PATH):
        os.makedirs(BASE_DOWNLOAD_PATH)

    #query = "%s/%s" % (QUERY_BASE, minid)
    #print "Executing query: %s" % query

    #r = requests.get(query,  headers = {"Accept" : "application/json"})
    #print r.json()

    #location = None
    #if len(minid) > 17:
    #    location = r.json()["location"][0]
    #else:
    #    location = r.json()["locations"][0]['link']
    filename = location.split("/")[-1]
    path = "%s/%s" % (BASE_DOWNLOAD_PATH, filename)
    sample_name = filename.split("_")[0]

    print "Downloading result: %s" % location
    print "Sample name: %s" % sample_name

    #https://bags.fair-research.org/ with s3://gg-commons/bags/
    location = location.replace("https://bags.fair-research.org/", "bags/")

    #testfile = urllib.URLopener()
    #testfile.retrieve(location, path)

######PULL DATA USING S3 format
    S3_BUCKET = 'gg-commons'
    session = boto3.session.Session(profile_name='gs3conn')
    s3_client = session.client("s3", config=Config(signature_version="s3v4"))
    s3_client.download_file(S3_BUCKET, location, path)

    # BDBag tooling doesnt let us unzip into a single dir.. so
    # we use a /base/uuid/uuid
    extract_path = ".".join(path.split(".")[0:-1])
    #output_path = "%s/%s" %(extract_path, filename.split(".")[0])
    output_path = "%s/%s" %(extract_path, os.path.basename(extract_path))
    uncomp_path = "%s/%s" %(extract_path, "real_data")

    print "Extracting bag and resolving fetch: %s" % output_path
    bdbag_api.extract_bag(path, extract_path)

    # if url starts with https, modify to s3 url
    fetch_file = "%s/fetch.txt" % output_path
    os.mkdir("%s/data" % output_path)
    fh = open(fetch_file, "r")
    get_objs = []
    for line in fh:
        line = line.replace("https://results.fair-research.org/","results/")
        results_location = line.split("\t")[0]
        filename = results_location.split("/")[-1]
        path = "%s/data/%s" % (output_path, filename)
        s3_client.download_file(S3_BUCKET, results_location, path)
    fh.close()

    #print "Resolving fetch: %s" % output_path
    #bdbag_api.resolve_fetch(output_path, True)
    print "Resolved fetch: %s" % output_path

    tmp_path = os.path.join(output_path, "data")
    for f in os.listdir(tmp_path):
        if f.endswith('.gz'):
            return ["%s/%s" % (tmp_path, f) , sample_name]

def sort_file(inputf):
    outf = "%s.sorted" % inputf
    os.system("head -n1 %s > %s; tail -n +2 %s |sort -k1,1 -k2,2h -k3,3 -k4,4  >> %s" % (inputf, outf, inputf, outf))
    return outf

def pull_minids(minid_set, output_dir):
    paths = []
    sample_name = None
    for minid in minid_set:
        print "PULLIN MINID: %s" % minid
        gzip_file, sample_name = download_file_from_minid(minid, output_dir)
        inF = gzip.GzipFile(gzip_file, 'rb')
        s = inF.read()
        inF.close()

        gunzip_file = gzip_file.replace(".gz", "")
        outF = file(gunzip_file, 'wb')
        outF.write(s)
        outF.close()
        os.remove(gzip_file)

        #sort file
        paths.append(sort_file(gunzip_file))
        os.remove(gunzip_file)
    return [paths,sample_name]

def merge_vcfs(anno_vcfs):
#    df1 = []
#    for vcf in anno_vcfs:
#        fh_header = open(vcf, "r")
#        header = fh_header.readline().rstrip("\n").split("\t")
#        dtype = {header[0] : 'string', header[1] : 'int64', header[2] : 'string', header[3] : 'string' }
#        for h in header[4:]:
#            dtype[h] = 'string'
#        fh_header.close()
#        df = pd.read_csv(vcf, sep="\t",dtype=dtype)
#        dfs.append(df)
#
#    new_df = pd.merge(df[0], df[1], on=['chr', 'pos', 'ref', 'alt'],how='outer')
#    for df in dfs[2:]:
#        new_df = new_df.merge(df, on=['chr', 'pos', 'ref', 'alt'],how='outer')
#
#    fd, path = tempfile.mkstemp()
#    df.to_csv(path, index=False, sep="\t")
#
#    return path

    fh_header1 = open(anno_vcfs[0], "r")
    header1 = fh_header1.readline().rstrip("\n").split("\t")
    dtype1 = {header1[0] : 'string', header1[1] : 'int64', header1[2] : 'string', header1[3] : 'string' }
    for h in header1[4:]:
        dtype1[h] = 'string'
    fh_header1.close()
    #di = dd.read_csv(anno_vcfs[0], sep="\t",dtype=dtype1)
    di = pd.read_csv(anno_vcfs[0], sep="\t",dtype=dtype1)
    path = None
    path_old = None
    for i in anno_vcfs[1:]:
        if path_old is not None:
            os.remove(path_old)
        if path is not None:
            path_old = path
        print "Reading paired file %s" % i

        fh_header2 = open(i, "r")
        header2 = fh_header2.readline().rstrip("\n").split("\t")
        dtype2 = {header2[0] : 'string', header2[1] : 'int64', header2[2] : 'string', header2[3] : 'string' }
        for h in header2[4:]:
            dtype2[h] = 'string'
        fh_header2.close()

        #dl = dd.read_csv(i, sep="\t",dtype=dtype2)
        dl = pd.read_csv(i, sep="\t",dtype=dtype2)
        print "mergin with %s" % i
        #df = dd.merge(di, dl, on=['chr', 'pos', 'ref', 'alt'],how='outer').compute()
        df = pd.merge(di, dl, on=['chr', 'pos', 'ref', 'alt'],how='outer')
        fd, path = tempfile.mkstemp()
        print "Writing merge to %s" % path
        df.to_csv(path, index=False, sep="\t")
        print "Opening merged file %s" % path
        fh_header1 = open(path, "r")
        header1 = fh_header1.readline().rstrip("\n").split("\t")
        dtype1 = {header1[0] : 'string', header1[1] : 'int64', header1[2] : 'string', header1[3] : 'string' }
        for h in header1[4:]:
            dtype1[h] = 'string'
        fh_header1.close()

        di = pd.read_csv(path, sep="\t",dtype=dtype1)

    return path

def stop_err( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit()

def is_complete (historyID, gi):
    status = gi.histories.get_status(historyID)
    if status['percent_complete'] == "100":
        return True
    elif status['state'] == 'ok':
        return True
    else:
        return False

parser = optparse.OptionParser()
parser.add_option( '-k', '--api-key', dest='api_key', help='Mandatory. Need to get your galaxy API key by visiting the Galaxy webpage http://dev.globusgenomics.edu' )
parser.add_option( '-u', '--url', dest='url_name', help='Galaxy URL' )
parser.add_option( '--input', dest='table_file', help='The table file containing the parameters for the workflow to submit' )
parser.add_option( '--out', dest='output', help='Log file' )
parser.add_option( '--out-dir', dest='output_dir', help='Output directory to store BAM, FASTQ, VCF and more' )
parser.add_option( '-t', '--token-auth', dest="goauth_token", help='Globus auth token' )
(options, args) = parser.parse_args()

url = options.url_name
if "http:" in options.url_name:
   url = options.url_name.replace("http", "https")

# get the database, user, host, password, url
#print os.getcwd()

output_dir = tempfile.mkdtemp(prefix="merge-")
tmp_dir = tempfile.mkdtemp(prefix="merge-tmp-")

#if not os.path.exists(options.output_dir):
#    os.mkdir(options.output_dir)


key = options.api_key
if len(key) > 0:
        # get an API handle
        gi = GalaxyInstance(url=url, key=key)

        # monitor the jobs
        # get the list of histories to monitor
        monitor_meta = {}
        fh = open(options.table_file, "r")
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("SUBMITTED"):
                status,sampleName,wfName,wfID,historyName,historyID = line.split("\t")
                monitor_meta[historyID] = { 'sampleName':sampleName, 'wfName':wfName, 'wfID':wfID, 'historyName':historyName }
               
        completed_meta = {}
        #summary_output = "%s/summary_output.txt" % options.output_dir
        #fh_out = open(summary_output, "w")
        snps_minids = []
        indels_minids = []
        missing_minids = []
        while len(monitor_meta) != len(completed_meta):
            for hKey in monitor_meta:
                print "Analyzing HKEY %s: %s" % (hKey, monitor_meta[hKey]['sampleName']) 
                if hKey in completed_meta.keys():
                    continue
                elif is_complete(hKey, gi) == True:
                    #fh_out.write("%s\n" % hKey)
                    contents = gi.histories.show_history(hKey, contents=True)
                    dataset_id = None
                    for content in contents:
                        if content['deleted'] is False and content['name'] == "Minid for history":
                            dataset_id = content['id']
                            break
                    #dataset_id = gi.histories.show_history(hKey, contents=True)[0]['id']
                    dataset = gi.histories.show_dataset(hKey, dataset_id)
                    dataset_file = gi.datasets.show_dataset(dataset['id'])['file_name']
                    data = None
                    if os.path.getsize(dataset_file) > 0:
                        with open(dataset_file) as data_file:
                            data = json.loads(data_file.read())
                            bdbag_url = data['BDbag_path']
                            minid = data['minid']

                    if data is None or 'minid' not in data:
                        job = gi.jobs.show_job(dataset['creating_job'], full_details=True)
                        bdbag_name = "NONE"
                        for line in job['stderr'].split("\n"):
                            if line.startswith("INFO:bdbag.bdbag_api:Created bag archive:"):
                                bdbag_name = line.split("/")[-1]
                                bdbag_url = "https://bags.fair-research.org/%s" % bdbag_name
                                minid = "MISSING"
                                completed_meta[hKey] = {'status' : 'complete', 'sampleName': monitor_meta[hKey]['sampleName'], 'minid' : "NONE"}
                                break
                        #missing_minids.append([hKey, monitor_meta[hKey]['sampleName'], bdbag_name])
                        if bdbag_name != "NONE":
                            completed_meta[hKey] = {'status' : 'complete', 'sampleName': monitor_meta[hKey]['sampleName'], 'minid' : "NONE"}
                            vtype = None
                            if '.snps.' in monitor_meta[hKey]['sampleName']:
                                snps_minids.append(["NA", bdbag_url])
                            elif '.indels.' in monitor_meta[hKey]['sampleName']:
                                indels_minids.append(["NA", bdbag_url])
                        else:
                            missing_minids.append([hKey, monitor_meta[hKey]['sampleName'], bdbag_name])
                    else:
                        completed_meta[hKey] = {'status' : 'complete', 'sampleName': monitor_meta[hKey]['sampleName'], 'minid' : minid}
                        vtype = None
                        if '.snps.' in monitor_meta[hKey]['sampleName']:
                            snps_minids.append([minid, bdbag_url])
                        elif '.indels.' in monitor_meta[hKey]['sampleName']:
                            indels_minids.append([minid, bdbag_url])

            if len(monitor_meta) != len(completed_meta):
                print "MISSING_MINIDS: %s" % missing_minids
                time.sleep(600)
        print "INDELS: %s " % indels_minids
        print "SNPS: %s " % snps_minids
        #if len(missing_minids) > 0:
        #    print "MISSING_MINIDS:"
        #    for i in missing_minds:
        #        print i
        #    #sys.exit("CAN'T KEEP GOING. FIX YOUR MINIDS")

        # pull minid data from wherever they are
        final_vcfs = []
        index = 0
        sample_name = None
        print "PULLING IN MINIDS:"
        for minid_set in [snps_minids, indels_minids]:
            anno_vcfs, sample_name = pull_minids(minid_set, output_dir)
            merged_vcf = merge_vcfs(anno_vcfs)
            #inF = file(merged_vcf, 'rb')
            #s = inF.read()
            #inF.close()

            if index == 0:
                vtype = "snp"
            else:
                vtype = "indel"
            final_merged_unzipped =  "%s/%s.out.vcf.annotatedresources.dir.%s" % (output_dir, sample_name, vtype)
            final_merged = "%s/%s.out.vcf.annotatedresources.dir.%s.gz" % (output_dir, sample_name, vtype)
            #outF = gzip.GzipFile(final_merged, 'wb')
            #outF.write(s)
            #outF.close()
            copyfile(merged_vcf, final_merged_unzipped)
            os.system("gzip %s" % final_merged_unzipped)
            final_vcfs.append(final_merged)
            index = 1

        # create BDbag/minid for final_vcfs
        ### s3 transfer of result files
        S3_BUCKET = 'gg-commons'
        results_url = "https://results.fair-research.org"
        results_path = "/results"
        bdbags_url = "https://bags.fair-research.org"
        bdbags_path = "/bags"

        identifier = str(int(time.time()*10000))
        s3_dir_to_use = "{0}_{1}".format( sample_name+"_1", identifier)
        files_to_bdbag = []
        remote_manifest = []
        session = boto3.session.Session(profile_name='gs3conn')
        s3_client = session.client("s3", config=Config(signature_version="s3v4"))
        for outf in final_vcfs:
            to_loc = "{0}/{1}".format(s3_dir_to_use, os.path.basename(outf))
            #to_path = "%s/%s" % (results_path, to_loc)
            real_url_path = "%s/%s" % (results_url, to_loc)
            file_size = os.path.getsize(outf)
            md5sum = hashlib.md5(file_as_bytes(open(outf, 'rb'))).hexdigest()
            filename = os.path.basename(outf)

            #bucket = s3.Bucket(S3_BUCKET)
            #url_path = "https://s3.amazonaws.com/%s/results/%s" % (S3_BUCKET,  to_loc)
            s3_client.upload_file(outf, S3_BUCKET, "results/%s" % (to_loc))
            remote_manifest.append({'url' : real_url_path, "length": file_size, "filename" : filename, "md5" : md5sum})

        #print remote_manifest
        tmp_file_path = os.path.join(output_dir, "dbbag_history_info_{0}".format(identifier))
        #print tmp_file_path
        with open(tmp_file_path, 'w') as fp:
            json.dump(remote_manifest, fp)
        bdbag_name = "%s.%s" % (s3_dir_to_use, "outputs.bdbag")
        bdbag_path = os.path.join(output_dir, bdbag_name)
        if not os.path.isdir(bdbag_path):
            os.makedirs(bdbag_path, 0755)
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

        service_token = check_globus_token(options.goauth_token, "/home/galaxy/.globusgenomics/globus_creds", "/home/galaxy/.globusgenomics", "sulakhe", "identifier")
        minid_title = "GTEX %s annotations output results" % sample_name

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
                        'Title': minid_title,
                        'erc.who': "Dinanath Sulakhe",
                        'erc.what' : os.path.basename(json.dumps([s3_url])),
                        '_profile': "erc", 
                        "contentSize": os.path.getsize(archived_path)
                })
        }
        response_with_minid = ic.create_identifier(**kwargs)
        minid = response_with_minid['identifier']

        out_dict = { 'minid' : minid, 'BDbag_files' : remote_manifest, 'BDbag_path' : s3_url, 'title': minid_title}
        fh_out = open(options.output, "w")
        fh_out.write( json.dumps(out_dict))
        fh_out.close()

