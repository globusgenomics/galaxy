import os
import functools
import urllib
from argparse import ArgumentParser
from globusonline.transfer.api_client import Transfer, TransferAPIClient
from datetime import datetime, timedelta
import sys
import time
import requests
output = sys.stdout

TIMEOUT = 10 # mins

def parse_cli():
    description = 'Fetch datasets from the Globus Online Catalog'
    parser = ArgumentParser(description=description)
    parser.add_argument('--session', required=True)
    parser.add_argument('--token', required=True)
    #parser.add_argument('--host', required=True)
    parser.add_argument('--epname', default='')
    parser.add_argument('--eppath', default='')
    parser.add_argument('--dest_epname', default='')
    parser.add_argument('--dest_eppath', default='')
    parser.add_argument('output_primary', nargs=1)
    parser.add_argument('output_id', nargs=1)
    parser.add_argument('output_dir', nargs=1)
    return parser.parse_args()



def wait_for_task(client, task_id, timeout):
    status = "ACTIVE"
    while timeout and status == "ACTIVE":
        code, reason, data = client.task(task_id, fields="status")
        status = data["status"]
        time.sleep(10)
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

def start_transfer(client, src_ep, src_path, dest_ep, dest_path):
    code, message, data = client.transfer_submission_id()
    submission_id = data["value"]
    deadline = datetime.utcnow() + timedelta(minutes=10)
    t = Transfer(submission_id, src_ep, dest_ep, deadline)
    t.add_item(src_path, dest_path, True)
    code, reason, data = client.transfer(t)
    task_id = data["task_id"]
    return task_id

def main():
    args = parse_cli()
    src_ep = decode_ep(args.epname)
    src_path = args.eppath
    output_primary = args.output_primary[0]
    output_id = args.output_id[0]
    output_dir = args.output_dir[0]
    token = args.token
    client = TransferAPIClient("XXX", goauth=token)
    dst_ep = decode_ep(args.dest_epname)
    dst_path = args.dest_eppath

    print "%s::%s -->> %s::%s" % (src_ep, src_path, dst_ep, dst_path)
    # Activate EPs
    _,_,_ = client.endpoint_autoactivate(src_ep, if_expires_in=600)
    _,_,_ = client.endpoint_autoactivate(dst_ep, if_expires_in=600)

    # Get dataset dirs to copy
    task_id = start_transfer(client, src_ep, src_path, dst_ep, dst_path)

    # Spin on transfer
    timeout = int(TIMEOUT) * 60  # 60 seconds
    wait_for_task(client, task_id, timeout)

    with open(output_primary, 'w') as ofile:
        ofile.write(dst_path)
        
    #cred = fetch_credentials(args.host, args.session, args.token)
    # we are going to ignore catalog for the demo purposes.
    #log_msg = ('Globus Online Catalog Tool Information\n\n' 
    #           'User: %s\n'
    #           'Catalog: %s\n'
    #           'Dataset: %s\n'
    #           'Selected files: %s\n') % \
    #           (cred['gouser'], 'BIRN Catalog', dataset_ref, files)
    #with open(output_log, 'w') as outlog:
    #    outlog.write(log_msg)
    #    outlog.write(('-' * 80) + '\n')
    #    make_transfer(cred['gouser'], cred['gotoken'], 
    #                  catalog, dataset_ref, outlog.write, output_id,
    #                  output_dir, files)


if __name__ == '__main__':
    main()
