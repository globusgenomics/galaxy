#!/usr/bin/python
#!/usr/bin/python
# -*- coding: utf-8 -*-

import time, optparse, os, shutil, subprocess, sys, glob, tempfile, psycopg2, socket, re
from bioblend.galaxy import GalaxyInstance
#from galaxy.web import url_for
import requests
requests.packages.urllib3.disable_warnings()

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
parser.add_option(  '--history', dest="history_id", help="the history api id" ) 
parser.add_option( '--out', dest='log', help='Log file' )
(options, args) = parser.parse_args()

url = options.url_name
if "http:" in options.url_name:
   url = options.url_name.replace("http", "https")

# get the database, user, host, password, url
#print os.getcwd()

key = options.api_key
if len(key) > 0:
    # get an API handle
    gi = GalaxyInstance(url=url, key=key)
    # monitor the jobs
    # get the list of histories to monitor
    completed_meta = {}
    fh_out = open(options.log, "w")
    # get metadata about the completed jobs
    for dataset in gi.histories.show_history(options.history_id, contents=True):
        #print dataset
        meta = gi.datasets.show_dataset(dataset['id'])
        job_meta = gi.jobs.show_job(meta['creating_job'], full_details=True)
        #print meta
        fh_out.write(str(job_meta['command_line']) + "\n")
        fh_out.write(str(job_meta['stdout']) + "\n")
        fh_out.write("================================END OF JOB=========================\n\n")
fh_out.close()
