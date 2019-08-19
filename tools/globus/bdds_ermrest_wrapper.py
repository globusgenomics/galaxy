import os
import functools
import urllib
import shutil
from argparse import ArgumentParser
from datetime import datetime, timedelta
import sys
import time
import requests
import urllib2
output = sys.stdout
from globusonline.transfer.api_client import Transfer, TransferAPIClient
from iobox.dams2bag import dams2bag
from iobox.bag2dams import bag2dams
import bagit
import uuid

TIMEOUT = 150000 
SLEEP_TIME = 30
GLOBUS_FETCH = "globusfetch.txt"
DEADLINE = 24 * 60

# Example download
# python ermrest-wrapper.py --type=download --token <TOKEN> --dst_ep kyle#analytics-vm --bag_path /home/kyle/tmp/globusbag --query_path "/A:=PPMI:PATIENT/ENROLL_STATUS=Enrolled"
#
# Example upload
#
# python ermrest-wrapper.py --type=upload --token <TOKEN> --bag_path /home/kyle/tmp/uploadbag --entity_path /entity/PPMI:BRAINMORPHSIGVEC --input_path PPMI/BRAINMORPHSIGVEC.csv 
#

def parse_cli():
    description = 'Fetch datasets from the BDDS Catalog'
    parser = ArgumentParser(description=description)
    parser.add_argument('--session', required=False)
    parser.add_argument('--token', required=False)
    parser.add_argument('--type', required=False)
    parser.add_argument('--cat_user', default='bdds')
    parser.add_argument('--cat_password', default='bddsdemo!')
    parser.add_argument('--cat_host', default='https://misd-vm-12.isi.edu')
    parser.add_argument('--cat_path', default='/ermrest/catalog/1')
    parser.add_argument('--bag_path', default='')
    parser.add_argument('--query_path', default='')
    #parser.add_argument('--output_name', default='') 
    #parser.add_argument('--output_format', default='') 
    parser.add_argument('--input_file', default='')
    parser.add_argument('--entity_path', default='')
    parser.add_argument('--input_path', default='')
    parser.add_argument('--dst_ep', default='')
    parser.add_argument('--dst_ep_path', default='')
    parser.add_argument('output_extra_files', nargs=1)
    parser.add_argument('output_primary', nargs=1)
    parser.add_argument('output_id', nargs=1)
    parser.add_argument('output_dir', nargs=1)
    return parser.parse_args()

def wait_for_task(client, task_id, timeout):
    status = "ACTIVE"
    while timeout and status == "ACTIVE":
        code, reason, data = client.task(task_id, fields="status")
        status = data["status"]
        time.sleep(SLEEP_TIME)
        timeout -= 1
    if status != "ACTIVE":
        print >>output,"Task %s complete!" % task_id
        return True
    else:
        print >>output,"Task still not complete after %d seconds" % timeout
        return False

def decode_ep(ep):
    ep = ep.replace("__pd__", "#")
    ep = ep.replace("__dq__", "");
    return ep

def create_transfer(client, fetch_file, src_ep, dst_ep):
    code, message, data = client.transfer_submission_id()
    submission_id = data["value"]
    #deadline = datetime.utcnow() + timedelta(minutes=DEADLINE)
    t = Transfer(submission_id, src_ep, dst_ep)#, deadline)
    for k,v in fetch_file.items():
        t.add_item(v['src_path'], v['dst_path'], False)
    code, reason, data = client.transfer(t)
    task_id = data["task_id"]
    return task_id

def split_ep(uri):
    res = uri.partition("/")
    return (res[0], res[2])

def decode_path(path):
    if path.startswith("X"):
        return "/~" + path[1:]

def start_transfer(args, fetch_map, src_ep, dst_ep):
     client = TransferAPIClient("XXX", goauth=args.token)
     _,_,_ = client.endpoint_autoactivate(src_ep, if_expires_in=600)
     _,_,_ = client.endpoint_autoactivate(dst_ep, if_expires_in=600)

     task_id = create_transfer(client, fetch_map, src_ep, dst_ep)

     # Spin on transfer
     timeout = int(TIMEOUT) 
     wait_for_task(client, task_id, timeout)

# Needed for the galaxy setup where users have restricted access
def change_file_permissions(path, user):
    for r,d,f in os.walk(path):
        os.chmod( r , 0777)

def parse_globus_fetch(args):
    bag_source = args.bag_path
    fetch_file_path = os.path.join(bag_source, GLOBUS_FETCH)
    data_dir_path = os.path.join(bag_source, "data")
    print "Checking for Globus fetch file %s" % fetch_file_path
    fetch_map = {}
    if not os.path.isfile(fetch_file_path):
        return
    with open(fetch_file_path, "r") as f:
        f.next()
        for l in f: 
             files = l.split('\t')
             src_ep, src_path = split_ep(files[0])
             fetch_map[files[0]] = {"src_ep" : src_ep, "src_path" : src_path, 
                                    "dst_ep" : args.dst_ep, "dst_path" : os.path.join(bag_source, files[2].strip())}
    change_file_permissions(data_dir_path, "USER?")
    start_transfer(args, fetch_map, src_ep, args.dst_ep)     

bag_metadata =  {
        "Source-Organization": "USC Information Sciences Institute, Informatics Systems Research Division",
        "Contact-Name": "Globus Genomics",
        "External-Description": "Auto generated bag",
        "Internal-Sender-Identifier": "USC-ISI-IRSD"
}

def create_catalog_config(args):
    return {"host": args.cat_host,
            "path": args.cat_path,
            "username": args.cat_user,
            "password": args.cat_password}

def create_queries_config(args):
    # Note: this is hard coded to the PPMI query types we expect
    # TODO work out the best way to pass in this info (not sure users can really 
    # be expected to do it themselves?
    return [
            {
                "query_path": "/entity/%s" % args.query_path,
                "schema_path": "/schema/PPMI/table/PATIENT",
                "output_name": "PPMI/PATIENT",
                "output_format": "csv"
            },
            {
                "query_path" : "/attribute/%s" % args.query_path + "/F:=DATA_FILES/URL:=F:uri,LENGTH:=F:bytes,FILENAME:=F:filepath", 
                "output_name" : "PPMI/BAM",
                "output_format" : "globusfetch"
            }
            ]

def create_entities_config(args):
    return [{
        "entity_path": args.entity_path,
        "input_path": args.input_path,
        "input_format": "csv"
        }]

def create_download_config(args):
    catalog = create_catalog_config(args)
    queries = create_queries_config(args)
    catalog["queries"] = queries
    return {"bag" : {
            "bag_path" : args.bag_path,
            "bag_metadata" : bag_metadata,
            },
            "catalog" : catalog}


def create_upload_config(args, bag_path):
    catalog = create_catalog_config(args)
    entities = create_entities_config(args)
    catalog["entities"] = entities
    return {"bag" : {
        "bag_path" : bag_path,
        "bag_metadata" : bag_metadata,
    },
    "catalog" : catalog}
import traceback
def main():
    args = parse_cli()

    if args.type.startswith("d"):
        print "Downloading %s" % args.query_path
        download_config = create_download_config(args)
        dams2bag.export_to_bag(download_config)
        parse_globus_fetch(args)
	# Move the downloaded bag into the extra files dir
	shutil.move(args.bag_path, args.output_extra_files[0])

    else:
        print "Uploading %s" % args.input_file
        bag_path = os.path.join("/scratch/galaxy/tmp", str(uuid.uuid4()))
        upload_file_path = os.path.join(bag_path, args.input_path)
        upload_dir = os.path.dirname(upload_file_path)

    	if os.path.exists(upload_dir):
            print "Directory exists %s" % upload_dir
            return 
	os.makedirs(upload_dir)
        print "Moving file %s to %s" % (args.input_file, upload_file_path)
        shutil.copy(args.input_file, upload_file_path)
        bag = bagit.make_bag(bag_path, bag_metadata)
        upload_config = create_upload_config(args, bag_path)
        print "Uploading config"
        print upload_config
        try:
            bag2dams.import_from_bag(upload_config)
        except Exception as e:
            traceback.print_exc

if __name__ == '__main__':
    main()

