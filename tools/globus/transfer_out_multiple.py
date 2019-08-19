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
           --source-ep=my#ep
           --source-path=/~/somefile.txt
           --destination-ep=other#ep
           --destination-path=/~/somefile.txt
           --deadline=5
"""
import sys, shutil
import time
from datetime import datetime, timedelta
import random, os

import globusonline.transfer.api_client as transfer_api

# TransferAPIClient instance.
api = None
output = sys.stdout

DEBUG = False
DEFAULT_DEADLINE = 10
DEFAULT_TIMEOUT = DEFAULT_DEADLINE


def transfer(source_ep, destination_ep, pairs, deadline, delete_after_transfer=False):
#def transfer(source_ep, source_filepath, destination_ep,
#             destination_filepath, deadline, final_destination):
    """
    Do a bunch of API calls and display the results. Does a small transfer
    between tutorial endpoints, but otherwise does not modify user data.

    Uses module global API client instance.
    """
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
    code, message, data = api.transfer_submission_id()
    transfer_id = data["value"]
    if deadline == "":
        deadline = DEFAULT_DEADLINE
        timeout = DEFAULT_TIMEOUT
    else:
        timeout = deadline
    deadline = datetime.utcnow() + timedelta(minutes=int(deadline))
    transfer = transfer_api.Transfer(transfer_id, source_ep, destination_ep, deadline)
    import re
    for pair in pairs:
        transfer.add_item(re.sub("/scratch", "",pair[0], count=1), pair[1])

    #transfer.add_item(source_filepath, destination_filepath)
    # make sure there are no more than 3 active transfers going

    while can_submit() == 0:
        print "running more than needed"
        time.sleep(random_number(300))
    code, reason, data = api.transfer(transfer)
    task_id = data["task_id"]
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
        if delete_after_transfer:
            for pair in pairs:
                os.remove(pair[0])

def can_submit():
    transfer_count = _count_globus_queue_transfers()

    if transfer_count >= 3:
        return 0
    else:
        return 1

def _count_globus_queue_transfers():
    q_count = None
    while q_count is None:
        try:
            code, reason, data = api.task_list(filter="status:ACTIVE,INACTIVE")
            q_count = len(data['DATA'])
            return q_count
        except:
            time.sleep(random_number(60))
    return q_count

def random_number(maxi):
    return random.randint(0,maxi)

def display_activation(endpoint_name):
    print >>output,"=== Endpoint pre-activation ==="
    display_endpoint(endpoint_name)
    print
    code, reason, result = api.endpoint_autoactivate(endpoint_name, if_expires_in=600)
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
    if status != "ACTIVE" and status != "FAILED":
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


def process_args(args=None, parser=None):
    from optparse import OptionParser
    if not parser:
        usage = "usage: %prog username -k KEY_FILE -c CERT_FILE"

        parser = OptionParser(usage=usage)
    parser.add_option("-b", "--base-url", dest="base_url",
                      help="alternate base URL", metavar="URL")
    parser.add_option("-o", "--output", dest="output",
                      help="write log output to PATH", metavar="PATH")
    parser.add_option("--local-ep", dest="local_ep",
                      help="Galaxy instance endpoint name")
    parser.add_option("--source-ep", dest="source_ep",
                      help="Endpoint to transfer from")
    parser.add_option("--destination-ep", dest="destination_ep",
                      help="Endpoint to transfer to")
    parser.add_option("-d", "--deadline", dest="deadline",
                      help="Deadline for transfer in minutes.")
    parser.add_option("-g", "--galaxy-dataset-id", dest="galaxy_dataset_id",
                      help="Galaxy Dataset Id For This Transfer")
    parser.add_option('-a', '--goauth-token', dest="goauth_token",
                      help="Use the Globus Access Token as the authentication" \
                      " method for this transfer")
    parser.add_option("--dataset", dest="datasets", action="append", type="string", nargs=3,
                      help="source_path destination_tmp_path final_path")
    parser.add_option("--delete-output", dest="delete_outputs", action="store_true",
                      help="delete files after succesful transfer")
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
    apikwargs['goauth'] = options.goauth_token.strip()
    api = transfer_api.TransferAPIClient(args[0], **apikwargs)
    output = options.output
    delete_option = False
    if options.delete_outputs:
        delete_option = True
    display_endpoints()
    # some users just want to create an empty file. This gives them the option to create an empty file
    # by "transferring" a file with the string "empty" in it.

    pairs = []
    for (src_path, final_path, log) in options.datasets:
        pairs.append([src_path, final_path, log])

    transfer(options.source_ep, options.destination_ep, pairs, options.deadline, delete_after_transfer=delete_option)
