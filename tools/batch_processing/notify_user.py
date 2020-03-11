#!/usr/bin/env python

"""
A wrapper script for deleting a history

python ./delete_history.py --input /scratch/dev/galaxy/files/002/dataset_2662.dat --email 'arodri7@globusonline.org' --userkey 7fb0b744057f607cf89a95b66e877601 --url 'https://dev.globusgenomics.org'

"""
import json, time, sys, optparse, os, random, smtplib
from bioblend import galaxy
from email.mime.text import MIMEText
import requests
requests.packages.urllib3.disable_warnings()

parser = optparse.OptionParser()
parser.add_option( '--input', dest='dataset_path', help='The table file containing the parameters for the workflow to submit' )
parser.add_option( '--userkey', dest='user_api_key', help="The user api key" )
parser.add_option( '--globus-email', dest='globus_email', help="The user email" )
parser.add_option( '--url', dest='galaxyurl', help='Base URL' )
parser.add_option( '--email', dest='email', help='users email' )
parser.add_option( '--output', dest='output_file', help='output file' )
parser.add_option( '--history_id', dest='history_id', help='history_id')
(options, args) = parser.parse_args()

globus_email = options.globus_email.replace("__at__", "@")

url = options.galaxyurl
if "http:" in url:
   url = url.replace("http", "https")

if "dev1" in url:
   url = url.replace("dev1", "dev")

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

print "User:\t%s\nID:\t%s\n" % (globus_email, user_id)

flag = 0

#print history['id']
#print history
#print gi.histories.show_history(history['id'], contents=True)
# NEED TO REWORK THIS SECTION!!
message = ""
history = gi.histories.show_history(history_id)
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

## send email
receivers = [globus_email]
subject = "Workflow Status: {0}".format(history['name'])
body_message = ""
if flag == 0:
    body_message += "Your jobs have completed without errors.\nHISTORY NAME: %s\nHISTORY ID: %s\n\n" % (history['name'], history_id)
    print "View your history at:  %s/history/view/%s" % (url, history_id)
elif flag == 2:
    body_message += "There were errors in some of your jobs in your history.\nHISTORY NAME: %s\nHISTORY ID: %s\n\n%s\n" % (history['name'], history_id, message)
    body_message += "View your history at:  %s/history/view/%s" % (url, history_id)
    print "History could not be deleted. I should send email to:  %s. There is an error in the workflow: %s" % (globus_email, message)

body_message += "\nVisit you GlobusGenomics instance at %s/history/view/%s to review your work.\n\nSincerely,\n\nAdmin user" % (url, history_id)


def send_email(to_addresses, subject, email_content):
    import boto3
    import ConfigParser
    Config = ConfigParser.ConfigParser()
    Config.read("/home/galaxy/.globusgenomics/aws_creds")
    ses_aws_access_key_id = Config.get('ses_iam', 'aws_access_key_id')
    ses_aws_secret_access_key = Config.get('ses_iam', 'aws_secret_access_key')

    ses_client = boto3.client("ses", 
                        region_name="us-east-1", 
                        aws_access_key_id=ses_aws_access_key_id,
                        aws_secret_access_key=ses_aws_secret_access_key)
    response = ses_client.send_email(
            Source="server@globusgenomics.org",
            Destination={"ToAddresses": to_addresses},
            Message={
                "Subject": {
                    "Data": subject,
                    "Charset": "UTF-8"
                },
                "Body": {
                    "Text": {
                        "Data": email_content,
                        "Charset": "UTF-8"
                    }
                }
            }
        )

send_email(receivers, subject, body_message)
