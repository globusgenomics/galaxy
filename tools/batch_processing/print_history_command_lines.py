#!/usr/bin/env python

"""
A wrapper script for deleting a history

python ./print_history_command_lines.py --userkey 7fb0b744057f607cf89a95b66e877601 --url 'https://dev.globusgenomics.org' --history_id 939393933948

"""
import glob, re, json, time, sys, optparse, os, random
from bioblend import galaxy
import requests
requests.packages.urllib3.disable_warnings()

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]

def search_for_line(file_path, pbs_list):
    for pbs in pbs_list:
        for line in open(pbs):
            if file_path in line:
                return pbs

def get_all_files_in_command(command):
    files = []
    for i in command.split():
        if "/scratch/galaxy/files" in i or "/scratch/go" in i:
            i = i.replace('"', '')
            i = i.replace(';', '')
            if "=" in i:
                files.append(i.split("=")[1])
            else:
                files.append(i)
    return files


parser = optparse.OptionParser()
parser.add_option( '--userkey', dest='user_api_key', help="The user api key" )
parser.add_option( '--url', dest='galaxyurl', help='Base URL' )
parser.add_option( '--output', dest='output_file', help='output file' )
parser.add_option( '--history_id', dest='history_id', help='history_id')
(options, args) = parser.parse_args()


url = options.galaxyurl
if "http:" in url:
   url = url.replace("http", "https")

user_api_key = options.user_api_key
history_id = options.history_id

# get an api handle and the user's id
pbs_path = "/opt/galaxy/database/pbs"
gi = galaxy.GalaxyInstance(url=url, key=user_api_key)

complete_options_output_list = list()

message = ""
seen = []
pbs_files = glob.glob("%s/*.sh" % (pbs_path))
#pbs_files.sort(key=natural_keys)
for history_item in gi.histories.show_history(history_id, contents=True, details="all"):
    dataset_metadata = gi.histories.show_dataset(history_id, history_item['id'])
    if dataset_metadata['state'] != 'error' and dataset_metadata['state'] == 'ok' and dataset_metadata['deleted'] == False :
        pbs_file = search_for_line(dataset_metadata['file_name'], pbs_files)
        if pbs_file is None:
            continue
        with open(pbs_file, 'r') as f:
            lines = f.readlines()
            cmd = lines[-1].split("cd /opt/galaxy")[0]
            stime = os.path.getctime(pbs_file)
            if cmd not in seen:
                #print "%s\n%s\n" % (dataset_metadata['name'], cmd)
                complete_options_output_list.append([dataset_metadata['file_name'], dataset_metadata['name'], pbs_file, cmd, stime])
                seen.append(cmd)

for tup in sorted(complete_options_output_list, key=lambda x: x[-1]):
    print "%s\n%s\n" % (tup[1], tup[3])
