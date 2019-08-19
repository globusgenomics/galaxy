#!/usr/bin/env python

# Copyright 2010 University of Chicago
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
Globus Online Transfer Tool

This is simply an example that produces a list of endpoints

python transfer.py USERNAME -k ~/.globus/userkey.pem \
           -c ~/.globus/usercert.pem \
           --source-ep=my#ep
           --source-path=/~/somefile.txt
           --destination-ep=other#ep
           --destination-path=/~/somefile.txt
           --deadline=5
"""
import sys, shutil, os
import time
from datetime import datetime, timedelta
from urllib2 import HTTPError
import random

import globusonline.transfer.api_client as transfer_api

# TransferAPIClient instance.
api = None
output = sys.stdout

DEBUG = False
DEFAULT_DEADLINE = 10
DEFAULT_TIMEOUT = DEFAULT_DEADLINE


def transfer(source_ep, source_filepath, destination_ep,
             destination_filepath, deadline, final_destination,
             path_type, extra_final_dir):
    """
    Do a bunch of API calls and display the results. Does a small transfer
    between tutorial endpoints, but otherwise does not modify user data.

    Uses module global API client instance.
    """
    # create extra data directory if transfering a directory to the galaxy endpoint
    #if "galaxy#" in destination_ep and path_type == 'directory' and not os.path.exists(extra_final_dir):
    #    os.makedirs(extra_final_dir)

    # See what is in the account before we make any submissions.
    print >>output,"=== Before tutorial ==="
    display_tasksummary(); print
    display_endpoints(); print
    # auto activate the endpoint, and display before/after.
    display_activation(source_ep)
    display_activation(destination_ep)
    print >>output,"=== Source Before transfer ==="
    display_ls(source_ep); print
    print >>output,"=== Destination Before transfer ==="
    display_ls(destination_ep); print
    # submit a transfer
    transfer_id = None
    while transfer_id is None:
        try:
            code, message, data = api.transfer_submission_id()
            transfer_id = data["value"]
        except HTTPError:
            time.sleep(random_number(60))
    if deadline == "":
        deadline = DEFAULT_DEADLINE
        timeout = DEFAULT_TIMEOUT
    else:
        timeout = deadline
    deadline = datetime.utcnow() + timedelta(minutes=int(deadline))
    transfer = transfer_api.Transfer(transfer_id, source_ep, destination_ep, deadline)
    # remove some dirs in the path before the transfer actaully happen. Mark X. 161102
    import re
    #destination_filepath = re.sub("/scratch","",destination_filepath,count=1)
    # if moving out of instance we need to remove /scratch from source_filepath
    # for now I will modify source path if source endpoint starts with galaxy#
    if "galaxy#" in source_ep:
        source_filepath = re.sub("/scratch","",source_filepath,count=1)
    else:
        destination_filepath = re.sub("/scratch","",destination_filepath,count=1)

    if path_type == 'directory':
        # get contents in directory
        # add each file or directory into transfer list
        #for item in list_dir(source_ep, path=source_filepath):
        #    src_file = "%s/%s" % (source_filepath, item[0])
        #    dest_file = "%s/%s" % (destination_filepath, item[0])
        #    if item[1] == 'file':
        #        transfer.add_item(src_file, dest_file)
        #    else:
        #        transfer.add_item(src_file, dest_file, recursive=True)
        transfer.add_item(source_filepath, destination_filepath, recursive=True)
    else:
        transfer.add_item(source_filepath, destination_filepath)

    task_id = None
    while task_id is None:
        try:
            code, reason, data = api.transfer(transfer)
            task_id = data["task_id"]
        except HTTPError:
            time.sleep(random_number(60))
    # see the new transfer show up
    print >>output,"=== After submit ==="
    display_tasksummary(); print
    display_task(task_id); print
    # wait for the task to complete, and see the summary and lists update
    # use user's deadline to specify timeout, or else timeout is fixed to
    # 120 seconds and GO transfer will turn green after 120 seconds,
    # even the transfer is not completed yet
    timeout = int(timeout) * 60  # change minutes to seconds
    if wait_for_task(task_id, timeout):
        print >>output,"=== After completion ==="
        display_tasksummary(); print
        display_task(task_id); print
        display_ls(destination_ep); print

    ## ASSUMPTION THAT ALL OUR ENDPOINT START WITH "galaxy#"
    #if "galaxy#" in destination_ep:
    #    # copy the file to the final destination path
    #    if path_type == 'file':
    #        shutil.copy(destination_filepath, final_destination)
    #    elif path_type == 'directory':
    #        shutil.copytree(destination_filepath, extra_final_dir)
    #        #delete_file = "%s/test.xyz" % extra_final_dir
    #        #os.remove(delete_file)
    #        directory_contents = get_filepaths(extra_final_dir)
    #        fh = open(final_destination, "w")
    #        fh.write("\n".join(directory_contents))
    #        fh.close()
    #
    #    # delete the temporary transferred files
    #    # submit a delete event
    #    delete_id = None
    #    while delete_id is None:
    #        try:
    #            code, message, data = api.transfer_submission_id()
    #            delete_id = data["value"]
    #        except HTTPError:
    #            time.sleep(random_number(60))
    #    if path_type == 'directory':
    #        delete = transfer_api.Delete(delete_id, destination_ep, deadline=deadline, recursive=True)
    #    else:
    #        delete = transfer_api.Delete(delete_id, destination_ep, deadline=deadline, recursive=True)
    #    delete_dir = "/".join(destination_filepath.split("/")[0:5])
    #    #delete.add_item(destination_filepath)
    #    delete.add_item(delete_dir)
    #    dtask_id = None
    #    while dtask_id is None:
    #        try:
    #            dcode, dreason, ddata = api.delete(delete)
    #            dtask_id =  ddata["task_id"]
    #        except HTTPError:
    #            time.sleep(random_number(60))
    #    if wait_for_task(dtask_id, timeout):
    #        print >>output,"Finished deleting temporary data"

def get_filepaths(directory):
    """
    This function will generate the file names in a directory 
    tree by walking the tree either top-down or bottom-up. For each 
    directory in the tree rooted at directory top (including top itself), 
    it yields a 3-tuple (dirpath, dirnames, filenames).
    """
    file_paths = []  # List which will store all of the full filepaths.

    # Walk the tree.
    for root, directories, files in os.walk(directory):
        for filename in files:
            # Join the two strings in order to form the full filepath.
            filepath = os.path.join(".", filename)
            file_paths.append(filepath)  # Add it to the list.

    return file_paths  # Self-explanatory.

def display_activation(endpoint_name):
    print >>output,"=== Endpoint pre-activation ==="
    display_endpoint(endpoint_name)
    print
    status = None
    while status is None:
        try:
            code, reason, result = api.endpoint_autoactivate(endpoint_name, if_expires_in=600)
            status = "success"
        except HTTPError:
            time.sleep(random_number(60))
    #if result.code.startswith("AutoActivationFailed"):
    #    print >>output,"Auto activation failed, ls and transfers will likely fail!"
    #print >>output,"result: %s (%s)" % (result.code, result.message)
    print >>output,"=== Endpoint post-activation ==="
    display_endpoint(endpoint_name)
    print


def display_tasksummary():
    if DEBUG:
        code, reason, data = api.tasksummary()
        print >>output,"Task Summary for %s:" % api.username
        for k, v in data.iteritems():
            if k == "DATA_TYPE":
                continue
            print >>output,"%3d %s" % (int(v), k.upper().ljust(9))


def display_tasks(max_age=None):
    """
    @param max_age: only show tasks requested at or after now - max_age.
    @type max_age: timedelta
    """
    if DEBUG:
        kwargs = {}
        if max_age:
            min_request_time = datetime.utcnow() - max_age
            # filter on request_time starting at min_request_time, with no
            # upper limit on request_time.
            kwargs["request_time"] = "%s," % min_request_time
        code, reason, tasks = api.task_list(**kwargs)
        print >>output,"Tasks for %s:" % api.username
        for task in tasks["DATA"]:
            print >>output,"Task %s:" % task["task_id"]
            _print_task(task)

def _print_task(data, indent_level=0):
    """
    Works for tasks and subtasks, since both have a task_id key
    and other key/values are printed by iterating through the items.
    """
    indent = " " * indent_level
    indent += " " * 2
    for k, v in data.iteritems():
        if k in ("DATA_TYPE", "LINKS"):
            continue
        print >>output,indent + "%s: %s" % (k, v)


def display_task(task_id, show_subtasks=True):
    if DEBUG:
        code, reason, data = api.task(task_id)
        print >>output,"Task %s:" % task_id
        _print_task(data, 0)
        if show_subtasks:
            code, reason, data = api.subtask_list(task_id)
            subtasks = data["DATA"]
            for t in subtasks:
                print >>output,"  subtask %s:" % t["task_id"]
                _print_task(t, 4)


def _count_queue_transfers():
    q_count = None
    while q_count is None:
        try:
            code, reason, data = api.task_list(filter="status:ACTIVE,INACTIVE")
            q_count = len(data['DATA'])
            return q_count
        except:
            time.sleep(random_number(60))
        

def wait_for_task(task_id, timeout):
    status = "ACTIVE"
    while timeout and status == "ACTIVE":
        mystat = None
        while mystat is None:
            try:
                code, reason, data = api.task(task_id, fields="status")
                mystat = "yes"
            except:
                time.sleep(random_number(60))
        status = data["status"]
        time.sleep(10)
        timeout = 1
    if status != "ACTIVE":
        print >>output,"Task %s complete!" % task_id
        return True
    else:
        print >>output,"Task still not complete after %d seconds" % timeout
        return False


def display_endpoint(endpoint_name):
    if DEBUG:
        code, reason, data = api.endpoint(endpoint_name)
        _print_endpoint(data)


def unicode_(data):
    """
    Coerce any type to unicode, assuming utf-8 encoding for strings.
    """
    if isinstance(data, unicode):
        return data
    if isinstance(data, str):
        return unicode(data, "utf-8")
    else:
        return unicode(data)

def list_dir(endpoint_name, path=""):
    ls_list = []
    print endpoint_name
    print path
    status = None
    while status is None:
        try:
            code, reason, data = api.endpoint_ls(endpoint_name, path)
            status = "SUCCESS"
        except HTTPError:
            time.sleep(random_number(60))
    # Server returns canonical path; "" maps to the users default path,
    # which is typically their home directory "/~/".
    path = data["path"]
    print path
    for f in data["DATA"]:
        print f
        ls_list.append([f['name'], f['type']])
    return ls_list

def display_ls(endpoint_name, path=""):
    if DEBUG:
        try:
            code, reason, data = api.endpoint_ls(endpoint_name, path)
        except Exception:
            print >>output, "Unable to list %s on %s:" % (path, endpoint_name)
            return
        # Server returns canonical path; "" maps to the users default path,
        # which is typically their home directory "/~/".
        path = data["path"]
        print >>output,"Contents of %s on %s:" % (path, endpoint_name)
        headers = "name, type, permissions, size, user, group, last_modified"
        headers_list = headers.split(", ")
        print >>output, headers
        for f in data["DATA"]:
            print >>output,", ".join([unicode_(f[k]) for k in headers_list])


def _print_endpoint(ep):
    name = ep["canonical_name"]
    print >>output,name
    if ep["activated"]:
        print >>output,"  activated (expires: %s)" % ep["expire_time"]
    else:
        print >>output,"  not activated"
    if ep["public"]:
        print >>output,"  public"
    else:
        print >>output,"  not public"
    if ep["myproxy_server"]:
        print >>output,"  default myproxy server: %s" % ep["myproxy_server"]
    else:
        print >>output,"  no default myproxy server"
    servers = ep.get("DATA", "")
    print >>output,"  servers:"
    for s in servers:
        print >>output,"    " + s["uri"],
        if s["subject"]:
            print >>output," (%s)" % s["subject"]
        else:
            print >>output, ""


def display_endpoints():
    if DEBUG:
        code, reason, endpoints = api.endpoint_list(limit=100)
        print >>output, ("Found %d endpoints for user %s:" %
                         (endpoints["length"], api.username))
        for ep in endpoints["DATA"]:
            _print_endpoint(ep)

def random_number(maxi):
    return random.randint(0,maxi)

def process_args(args=None, parser=None):
    from optparse import OptionParser
    if not parser:
        usage = "usage: %prog username -k KEY_FILE -c CERT_FILE"

        parser = OptionParser(usage=usage)
    parser.add_option("-c", "--cert", dest="cert_file",
                      help="client cert file", metavar="CERT_FILE")
    parser.add_option("-k", "--key", dest="key_file",
                      help="client key file", metavar="KEY_FILE")
    parser.add_option("-b", "--base-url", dest="base_url",
                      help="alternate base URL", metavar="URL")
    parser.add_option("-o", "--output", dest="output",
                      help="write log output to PATH", metavar="PATH")
    parser.add_option("--source-ep", dest="source_ep",
                      help="Endpoint to transfer from")
    parser.add_option("--source-path", dest="source_path",
                      help="Source endpoint filepath to transfer")
    parser.add_option("--extra-source-path", dest="extra_source_path",
                      help="Source endpoint filepath to transfer for BAM Index file")
    parser.add_option("--destination-ep", dest="destination_ep",
                      help="Endpoint to transfer to")
    parser.add_option("--destination-path", dest="destination_path",
                      help="Destination endpoint filepath to transfer")
    parser.add_option("-d", "--deadline", dest="deadline",
                      help="Deadline for transfer in minutes.")
    parser.add_option("-g", "--galaxy-dataset-id", dest="galaxy_dataset_id",
                      help="Galaxy Dataset Id For This Transfer")
    parser.add_option('-a', '--goauth-token', dest="goauth_token",
                      help="Use the Globus Access Token as the authentication" \
                      " method for this transfer")
    parser.add_option("--final", dest="final_dest",
                      help="final resting spot", metavar="FINAL")
    parser.add_option("--final_extra", dest="final_extra_path",
                      help="final extra directory path when transferring a directory")
    parser.add_option("--type", dest="path_type",
                      help="directory or file")

    parser.set_defaults(base_url=transfer_api.DEFAULT_BASE_URL)
    parser.set_defaults(output=sys.stdout)
    options, args = parser.parse_args(args)
    if len(args) != 1:
        parser.error("username arguments is required")
    if options.output != sys.stdout:
        options.output = open(options.output, "w")
    return options, args


if __name__ == '__main__':
    options, args = process_args()
    apikwargs = {'base_url':  options.base_url,}
    if options.goauth_token is None:
        apikwargs['cert_file'] = options.cert_file,
        apikwargs['key_file'] = options.key_file
    else:
        apikwargs['goauth'] = options.goauth_token.strip()
    api = transfer_api.TransferAPIClient(args[0], **apikwargs)

    ## wait for submission until the user only has less than 3 active transfers going on
    while _count_queue_transfers() >= 3:
        time.sleep(random_number(60))

    output = options.output
    display_endpoints()
    transfer(options.source_ep, options.source_path,
             options.destination_ep, options.destination_path,
             options.deadline, options.final_dest, options.path_type, options.final_extra_path)

    if options.extra_source_path != None:
        bai_dest_path = "%s.bai" % options.destination_path
        transfer(options.source_ep, options.extra_source_path,
             options.destination_ep, bai_dest_path,
             options.deadline, options.final_dest, options.path_type, options.final_extra_path)

