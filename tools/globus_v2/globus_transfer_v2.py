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
import ConfigParser

from optparse import OptionParser
args=None
parser = OptionParser()


parser.add_option("--username", dest="username", help="user short")
parser.add_option("--globus-cred-file", dest="globus_cred_file", help="globus cred file path")
parser.add_option("--from-endpoint", dest="from_endpoint", help="from endpoint")
parser.add_option("--to-endpoint", dest="to_endpoint", help="to endpoint")
parser.add_option("-i", dest="items", help="from to paths in from::to::to in list", action="append")
parser.add_option("--goauth-token", dest="goauth_token", help="goauth_token")
parser.add_option("--transfer-direction", dest="transfer_direction", help="in or out")
parser.add_option("--extra-source-path", dest="extra_source_path", help="extra_source_path", default=None)

options, args = parser.parse_args(args)

#print options

jobs_list = []
if options.items != None:
    for i in options.items:
        job = {}
        i_tmp = i.split("::to::")
        job["from_path"] = i_tmp[0].strip()
        job["to_path"] = i_tmp[1].strip()
        jobs_list.append(job)

transfer_info = {'username': options.username.strip(), 
                 'globus_cred_file': options.globus_cred_file.strip(),
                 'goauth_token': options.goauth_token.strip(), 
                 'from_endpoint': options.from_endpoint.strip(),
                 'to_endpoint': options.to_endpoint.strip(),
                 'jobs': jobs_list}


if 'username' in transfer_info:
    print 'username: {0}'.format(transfer_info['username'])
if 'from_endpoint' in transfer_info:
    print 'from_endpoint: {0}'.format(transfer_info['from_endpoint'])
if 'to_endpoint' in transfer_info:
    print 'to_endpoint: {0}'.format(transfer_info['to_endpoint'])
if 'jobs' in transfer_info:
    print 'jobs: {0}'.format(transfer_info['jobs'])


# verify the token and refresh it if expired
access_token = transfer_info['goauth_token']

Config = ConfigParser.ConfigParser()
Config.read(transfer_info['globus_cred_file'])
client_id = Config.get('globus', 'client_id')
client_secret = Config.get('globus', 'client_secret')
auth_client = globus_sdk.ConfidentialAppAuthClient(client_id, client_secret)
if not auth_client.oauth2_validate_token(access_token)['active']:
    refresh_token_file_name = transfer_info['username'].strip() + "_refresh"
    refresh_token_file = os.path.join(os.path.join(os.path.dirname(transfer_info['globus_cred_file']), "tokens"), refresh_token_file_name)
    with open(refresh_token_file, 'r') as f:
        refresh_token = f.readlines()[0].strip()
    r = auth_client.oauth2_refresh_token(refresh_token)
    access_token = r.by_resource_server['transfer.api.globus.org']['access_token']


transfer_direction = options.transfer_direction

transfer_client = globus_sdk.TransferClient(authorizer=globus_sdk.AccessTokenAuthorizer(access_token))

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
            sys.exit("The source path {0} does not exist. Verify you have the correct path; {1}".format(source_path, vars(e)))

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
        # check dataset affiliated dir again

        dir_link = source_path[0:-4] + "_files"
        if os.path.islink(dir_link):
            dir_realpath = os.path.realpath(dir_link)
            source_path = dir_realpath
        elif os.path.isdir(dir_link):
            source_path = dir_link
        source_path = source_path.rstrip('/').rstrip('\\')
        
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
