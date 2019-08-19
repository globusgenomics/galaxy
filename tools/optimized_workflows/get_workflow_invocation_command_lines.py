#!/usr/bin/env python

import os, sys, urllib2, csv, string, time, optparse, re, socket, errno
from bioblend.galaxy import GalaxyInstance
from itertools import *
sys.path.insert( 0, os.path.dirname( __file__ ) )
import requests
requests.packages.urllib3.disable_warnings()

# Get all command lines from a succesful workflow invocation
# Inputs: workflow id, workflow invocation id

def request (url, key):
    kwargs = dict()
    kwargs['params'] = {'key': key}
    kwargs.setdefault('verify', True)
    r = requests.get(url, **kwargs)
    if r.status_code == 200:
        return r.json()
    # @see self.body for HTTP response body
    raise ConnectionError("Unexpected response from galaxy: %s" %
                          r.status_code, body=r.text) 

def get_workflow (base_url, wfid, key):
    url = "%s/api/workflows/%s" % (base_url, wfid)
    wf = request(url, key)
    return wf

def get_workflow_invocations (base_url, wfid, key):
    url = "%s/api/workflows/%s/invocations" % (base_url, wfid)
    invocations = request(url, key)
    return invocations

def get_workflow_invocation (base_url, wfid, key, invocation_id):
    url = "%s/api/workflows/%s/invocations/%s" % (base_url, wfid, invocation_id)
    invocation = request(url, key)
    return invocation

def get_job (base_url, job_id, key):
    url = "%s/api/jobs/%s" % (base_url, job_id)
    job = request(url, key)
    return job

def get_dataset (base_url, hda, key):
    url = "%s/api/datasets/%s" % (base_url, hda)
    dataset = request(url, key)
    return dataset

def get_cl (job):
    return job['command_line']

def get_job_inputs (job):
    return job['inputs']

def get_job_outputs (job):
    return job['outputs']

def get_filename (dataset):
    return dataset['file_name']


parser = optparse.OptionParser()
parser.add_option( '-k', '--api-key', dest='api_key', help='Mandatory. Need to get your galaxy API key by visiting the Galaxy webpage http://dev.globusgenomics.edu' )
parser.add_option( '-u', '--url', dest='url_name', help='Galaxy URL' )
parser.add_option( '-w', '--workflow', dest='workflow_id', help='workflow_id' )
parser.add_option( '-i', '--invocation', dest='invocation_id', help='invocation_id' )
(options, args) = parser.parse_args()

try:
    options.api_key
    options.workflow_id
    options.url_name
    options.invocation_id
except NameError: 
    print 'usage: %s --workflow <workflow_id> --api-key <api_key> --url <galaxy_url> --invocation <invocation_id>' % os.path.basename( sys.argv[0] )
    sys.exit(1)

if options.api_key == "None":
    print 'ERROR: User does not have a valid API Key. Please create one by clicking on the User Menu directory'
    sys.exit(1)

if "http:" in options.url_name:
   options.url_name = options.url_name.replace("http", "https")

# get workflow
workflow = get_workflow(options.url_name, options.workflow_id, options.api_key)

# get invocation
invocation = get_workflow_invocation(options.url_name, options.workflow_id, options.api_key, options.invocation_id)

# get workflow input datasets
wf_inputs = []
for wf_input in invocation['inputs']:
    #print invocation['inputs'][wf_input]
    wf_inputs.append(invocation['inputs'][wf_input]['id'])

print wf_inputs

# loop through invocation jobs
for step in invocation['steps']:
    if step['job_id']: 
        job = get_job(options.url_name, step['job_id'], options.api_key)

        # get job command lines
        cl = get_cl(job)
        print cl

        # get job inputs
        inputs = get_job_inputs(job)
        print inputs
        for obj in inputs:
            dataset = get_dataset (options.url_name, inputs[obj]['id'], options.api_key)
            filename = get_filename (dataset)
            print filename

        # get job outputs
        outputs = get_job_outputs(job)
        print outputs
        for obj in outputs:
            dataset = get_dataset (options.url_name, outputs[obj]['id'], options.api_key)
            filename = get_filename (dataset)
            print filename


        print "\n"

# associate the inputs with the command lines

# return a json object of workflow command lines
