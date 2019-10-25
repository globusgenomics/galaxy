#!/usr/bin/python
import os
import time
import sys
import re
import globus_sdk
from datetime import datetime
import requests.packages.urllib3
requests.packages.urllib3.disable_warnings()

import ast

from optparse import OptionParser
args=None
parser = OptionParser()

parser.add_option("--transfer-info", dest="transfer_info", help="Detailed Transfer Information")
parser.add_option("--transfer-direction", dest="transfer_direction", help="in or out")
parser.add_option("--extra-source-path", dest="extra_source_path", help="extra_source_path", default=None)

options, args = parser.parse_args(args)

tmp_file_dir = '/home/galaxy/.globusgenomics/tmp'
transfer_info_file_name = options.transfer_info
transfer_info_file_loc = os.path.join(tmp_file_dir, transfer_info_file_name) 
with open(transfer_info_file_loc, 'r') as f:
    lines = f.readlines()
    try:
        transfer_info_read = lines[0].strip()
    except:
        sys.exit("failed to read transfer_info")
os.remove(transfer_info_file_loc)

transfer_info = ast.literal_eval(transfer_info_read)
print transfer_info
transfer_direction = options.transfer_direction

transfer_client = globus_sdk.TransferClient(authorizer=globus_sdk.AccessTokenAuthorizer(transfer_info['goauth_token']))

destination_ep = transfer_info['to_endpoint']
source_ep = transfer_info['from_endpoint']

# Verify the endpoints
def get_ep(ep):
    # take a enpoint name or id, return the enpoint object, retuen None if not found
    try:
        # try as ep id
        ep_object = transfer_client.get_endpoint(ep)
    except:
        ep_object = None
    if ep_object == None:
        # try as ep name
        for i in transfer_client.endpoint_search(ep):
            if i["display_name"] == ep or i['canonical_name'] == ep:
                ep_object = transfer_client.get_endpoint(i['id'])
    if ep_object != None:
        # activate the endpoint if needed,
        # if auto-activate not enabled, fail job and ask user to activate
        r = transfer_client.endpoint_autoactivate(ep_object['id'], if_expires_in=3600)
        if r["code"] == "AutoActivateFailed":
            sys.exit("Endpoint %s needs to be manually activated.\nPlease log on to http://globus.org to activate your endpoint.\n" % (ep))
        if ep_object["is_globus_connect"]:
            if ep_object["gcp_connected"] is False:
                sys.exit("Endpoint %s is disconnected. Please connect your endpoint\n" % (ep))
        else:
            if ep_object['DATA'] != [] and ep_object['DATA'][0]['is_connected'] is False:
                sys.exit("Endpoint %s is disconnected. Please connect your endpoint\n" % (ep))

    return ep_object

s_ep = get_ep(source_ep)
if s_ep is None:
    sys.exit("Could not find endpoint: %s. Make sure you entered the correct endpoint." % source_ep)
d_ep = get_ep(destination_ep)
if d_ep is None:
    sys.exit("Cound not find destination endpoint %s. Are you sure you are allowed to transfer to this endpoint?" % destination_ep)


tdata = globus_sdk.TransferData(transfer_client, s_ep['id'], d_ep['id'])

def add_job_to_queue(source_path, destination_path):
    # Verify the source path
    try:
        # If this is a directory, transfer a directory
        transfer_client.operation_ls(s_ep['id'], path=source_path)
        transfer_type = "dir"
    except globus_sdk.exc.TransferAPIError as e:
        # this is a file or path does not exist
        if e[1] == "ExternalError.DirListingFailed.NotDirectory" or e[1] == "ClientError.NotFound":
            # transfer a file
            transfer_type = "file"
        else:
            # path does not exist
            sys.exit("The source path %s does not exist. Verify you have the correct path" % source_path)

    if transfer_type == "dir":
        tdata.add_item(source_path, destination_path, recursive=True)
    elif transfer_type == "file":
        tdata.add_item(source_path, destination_path)
    else:
        sys.exit("transfer failed to find source path")


for job in transfer_info['jobs']:

    source_path = job['from_path']
    destination_path = job['to_path']

    if transfer_direction == 'in':
        #destination_path = re.sub("/scratch","",destination_path,count=1)
        pass
    elif transfer_direction == 'out':
        #source_path = re.sub("/scratch","",source_path,count=1)
        try:
            transfer_client.operation_ls(d_ep['id'], path=destination_path)
            destination_path = os.path.join(destination_path, os.path.basename(source_path))
        except globus_sdk.exc.TransferAPIError as e:
            pass
        # transfer extra bam
        if options.extra_source_path != None:
            #extra_source_path = re.sub("/scratch","",options.extra_source_path,count=1)
            extra_source_path = options.extra_source_path
            extra_destination_path = destination_path + '.bai'
            add_job_to_queue(extra_source_path, extra_destination_path)
    else:
        sys.exit("transfer direction not specified")
    add_job_to_queue(source_path, destination_path)
    

transfer_result = transfer_client.submit_transfer(tdata)

# wait while transfer is active or exit with error message if there is an error
stat_lines = []
while not transfer_client.task_wait(transfer_result['task_id'], timeout=60):
    status = transfer_client.get_task(transfer_result['task_id'])
    #print status
    nice_status_short_description = None
    if status['nice_status_short_description'] is None:
        nice_status_short_description = "active - no errors"
    else:
        nice_status_short_description = status['nice_status_short_description']

    print("%s\t%s\t%s\t%s\t%s" % (datetime.now().strftime("%Y-%m-%d %H:%M"), transfer_result['task_id'], status['status'], status['bytes_transferred'], nice_status_short_description))

    # if there is a permissions denied error, kill the transfer job and report error
    if status['nice_status_short_description'] == "permission denied":
        transfer_client.cancel_task(transfer_result['task_id'])
        continue

status = transfer_client.get_task(transfer_result['task_id'])
print "\n".join(stat_lines)
if status['status'] != "SUCCEEDED":
    sys.exit("Transfer failed: View status log")
