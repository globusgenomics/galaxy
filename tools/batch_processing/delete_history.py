#!/usr/bin/env python

"""
A wrapper script for deleting a history

python ./delete_history.py --input /scratch/dev/galaxy/files/002/dataset_2662.dat --email 'arodri7@globusonline.org' --userkey 7fb0b744057f607cf89a95b66e877601 --url 'https://dev.globusgenomics.org'

"""
import json, time, sys, optparse, os, random
from bioblend import galaxy
import requests
requests.packages.urllib3.disable_warnings()

parser = optparse.OptionParser()
parser.add_option( '--input', dest='dataset_path', help='The table file containing the parameters for the workflow to submit' )
parser.add_option( '--userkey', dest='user_api_key', help="The user api key" )
parser.add_option( '--url', dest='galaxyurl', help='Base URL' )
parser.add_option( '--email', dest='email', help='users email' )
parser.add_option( '--output', dest='output_file', help='output file' )
parser.add_option( '--history_id', dest='history_id', help='history_id')
(options, args) = parser.parse_args()


url = options.galaxyurl
if "http:" in url:
   url = url.replace("http", "https")

admin_api_key = '264c4c2466858acbaf323fc7b0ade690'
user_api_key = options.user_api_key
if len(user_api_key) < 5:
    sys.exit(errno.EACCES + ': Please create a user api key')

email = options.email
user_id = None
dataset_path = options.dataset_path
history_id = options.history_id

# get an api handle and the user's id
gi = galaxy.GalaxyInstance(url=url, key=user_api_key)
for user_dict in gi.users.get_users():
    if email == user_dict['email']:
        user_id = user_dict['id']

print "User:\t%s\nID:\t%s\n" % (email, user_id)

flag = 0

#print history['id']
#print history
#print gi.histories.show_history(history['id'], contents=True)
# NEED TO REWORK THIS SECTION!!
message = ""
for history_item in gi.histories.show_history(history_id, contents=True):
    dataset_metadata = gi.histories.show_dataset(history_id, history_item['id'])
    print dataset_metadata
    if dataset_metadata['state'] == 'error' and dataset_metadata['deleted'] == False:
        flag = 2
        message = dataset_metadata
        break
    elif dataset_metadata['deleted'] == False:
        print "checking:"
        print dataset_metadata
        while dataset_metadata['state'] != 'ok' and dataset_metadata['state'] != 'error' and dataset_metadata['file_name'] != options.output_file:
            flag = 1
            print dataset_metadata
            time.sleep(60)
            dataset_metadata = gi.histories.show_dataset(history_id, history_item['id'])
        flag = 0

print "FLAG: %s" % flag
## delete the history, THIS WORKS
if flag == 0:
    print "deleting"
    print gi.histories.delete_history(history_id, purge=True)
elif flag == 2:
    print "History could not be deleted. There is an error in the workflow:"
    print message
