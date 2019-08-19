#!/usr/bin/python
#!/usr/bin/python
# -*- coding: utf-8 -*-

import time, optparse, os, shutil, subprocess, sys, glob, tempfile, psycopg2, socket, re
#from galaxy.web import url_for
from bioblend.galaxy import GalaxyInstance
import json
import requests
requests.packages.urllib3.disable_warnings()
import urllib
from shutil import copyfile, rmtree
import requests.packages.urllib3
requests.packages.urllib3.disable_warnings()
sys.path.append("/mnt/galaxyTools/tools/pymodules/python2.7/lib/python/globus_sdk-1.5.0-py2.7.egg")
import globus_sdk
from datetime import datetime
import hashlib
import uuid

def stop_err( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit()

def is_complete (historyID, gi):
    status = gi.histories.get_status(historyID)
    if status['percent_complete'] == "100":
        return True
    elif status['state'] == 'ok':
        return True
    else:
        return False

parser = optparse.OptionParser()
parser.add_option( '-k', '--api-key', dest='api_key', help='Mandatory. Need to get your galaxy API key by visiting the Galaxy webpage http://dev.globusgenomics.edu' )
parser.add_option( '-u', '--url', dest='url_name', help='Galaxy URL' )
parser.add_option( '--input', dest='table_file', help='The table file containing the parameters for the workflow to submit' )
parser.add_option( '-t', '--token-auth', dest="goauth_token", help='Globus auth token' )
(options, args) = parser.parse_args()

url = options.url_name
if "http:" in options.url_name:
   url = options.url_name.replace("http", "https")

# get the database, user, host, password, url
#print os.getcwd()

try:
    key = options.api_key
    if len(key) > 0:
        # get an API handle
        gi = GalaxyInstance(url=url, key=key)

        # monitor the jobs
        # get the list of histories to monitor
        monitor_meta = {}
        fh = open(options.table_file, "r")
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("SUBMITTED"):
                status,sampleName,wfName,wfID,historyName,historyID = line.split("\t")
                monitor_meta[historyID] = { 'sampleName':sampleName, 'wfName':wfName, 'wfID':wfID, 'historyName':historyName }
               
        completed_meta = {}
        #summary_output = "%s/summary_output.txt" % options.output_dir
        #fh_out = open(summary_output, "w")
        snps_minids = []
        indels_minids = []
        missing_minids = []
        while len(monitor_meta) != len(completed_meta):
            for hKey in monitor_meta:
                print "Analyzing HKEY %s: %s" % (hKey, monitor_meta[hKey]['sampleName']) 
                if hKey in completed_meta.keys():
                    continue
                elif is_complete(hKey, gi) == True:
                    completed_meta[hKey] = {'status' : 'complete', 'sampleName': monitor_meta[hKey]['sampleName']}

            if len(monitor_meta) != len(completed_meta):
                time.sleep(600)

except Exception as e:
    print "There is some kind of issue here: %s" % (e)
    sys.exit()



