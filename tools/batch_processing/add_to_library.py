#!/usr/bin/env python

"""
A wrapper script for moving a dataset to a shared library 
If the library does not exist, it will create on and give access only to the correct user

python ./add_to_library.py --input /scratch/dev/galaxy/files/002/dataset_2662.dat --email 'arodri7@globusonline.org' --userkey 7fb0b744057f607cf89a95b66e877601 

"""
import json, time, sys, optparse, random
from bioblend import galaxy
import requests
requests.packages.urllib3.disable_warnings()

#wait_time = random.randrange(0,600)
#time.sleep(wait_time)

parser = optparse.OptionParser()
parser.add_option( '--input', dest='dataset_paths', action="append", help='The table file containing the parameters for the workflow to submit' )
parser.add_option( '--userkey', dest='user_api_key', help="The user api key" )
parser.add_option( '--url', dest='galaxyurl', help='Base URL' )
parser.add_option( '--email', dest='email', help='users email' )
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
dataset_paths = options.dataset_paths
history_id = options.history_id

#new_datasetname = 'MHUIC_snpeff'

gi_admin = galaxy.GalaxyInstance(url=url, key=admin_api_key)
gi =  galaxy.GalaxyInstance(url=url, key=user_api_key)
for user_dict in gi.users.get_users():
    if email == user_dict['email']:
        user_id = user_dict['id']

print "User:\t%s\nID:\t%s\n" % (email, user_id)

# get the user's histories and then look for the dataset with the correct path
histories = gi.histories.get_histories()
#print histories
history = gi.histories.get_histories(history_id=history_id)[0]
dataset = None
hdas = []
break_flag = 0
for dataset_dict in gi.histories.show_history(history_id, contents=True, details='all'):
    #print dataset_dict
    #print gi.histories.show_dataset(history_id, dataset_dict['id'])
    if 'file_name' in gi.histories.show_dataset(history_id, dataset_dict['id']):
        for dataset_path in dataset_paths:
            if gi.histories.show_dataset(history_id, dataset_dict['id'])['file_name'] == dataset_path:
                dataset = dataset_dict
                hdas.append(dataset_dict['id'])
                print gi.histories.show_dataset(history_id, dataset_dict['id'])['file_name']

# generate the new library name from the history name.
# Split the history name by the ~ characters and use the 1 and 2nd value.
new_libname = "%s-%s" % (history['name'].split('~')[0], history['name'].split('~')[1])
lib_users = [user_id]

# Check if user has an existing library with the name. If not, create one and add dataset
libraries = gi.libraries.get_libraries(name=new_libname)
if len(libraries) < 1:
    library = gi_admin.libraries.create_library(new_libname, description=None, synopsis=None)
    gi_admin.libraries.set_library_permissions(library['id'], access_in=lib_users, modify_in=lib_users, add_in=lib_users, manage_in=lib_users)
else:
    library = libraries[0]

## rename dataset in history. THIS DOES NOT WORK YET ON THIS VERSION OF GALAXY
#payload = {'name': new_datasetname}
#output = gi.histories.update_dataset(history['id'], dataset['id'], name=new_datasetname, test='something')
#print output

## add the dataset to the library
#hdas = [dataset['id']]
gi.libraries.upload_from_history(library['id'], hdas)

## delete the history, THIS WORKS
#print gi.histories.delete_history(history['id'], purge=True)

## current version of Galaxy cannot get user's API key or create one
##print gi.users.create_user_apikey(user_id)


