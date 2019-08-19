#!/usr/bin/env python

from bioblend.galaxy import GalaxyInstance
from itertools import *
import errno, csv, urllib2, json, time, sys, optparse, os, tempfile, subprocess, shutil, psycopg2, socket, re
from binascii import unhexlify
from string import Template
import requests
requests.packages.urllib3.disable_warnings()

def __main__():
    print "Start time:"
    print time.strftime('%d/%m/%Y %H:%M:%S\n', time.localtime(time.time()))

    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '-w', '--workflow-name', dest='workflow_name', action='store', type="string", help='Workflow name' )
    parser.add_option( '--history-id', dest='history_id', action='store', type="string", help='History id' )
    parser.add_option( '--history-item', dest='history_item', action='store', type="string", help='Item in History name' )
    parser.add_option( '-i', '--input-dir', dest='input_dir', action='store', type="string", help='Input FASTQ directory' )
    parser.add_option( '-o', '--output', dest='output_batch', action='store', type="string", help='output_batch' )
    parser.add_option( '-k', '--key', dest='api_key', help='Mandatory. Need to get your galaxy API key by visiting the Galaxy webpage http://dev.globusgenomics.edu' )
    parser.add_option( '--url', dest='url_name', help='Galaxy URL' )
    parser.add_option( '--username', dest='username', help='The user name from the postgresql table associated with the user' )
    (options, args) = parser.parse_args()

    # get api_url
    api_url = options.url_name
    if "http:" in api_url:
        api_url = api_url.replace("http", "https")

    # get user's API key
    api_key = options.api_key

    # connect to Galaxy vi API
    gi = GalaxyInstance(url=api_url, key=api_key)
    project_name = "QC_batch_per_sample"
    wf = gi.workflows.get_workflows(name=options.workflow_name)
    workflow_id = wf[0]['id']
    
    # get the param_id
    wf_input_id = None
    tool_id = None
    wf_steps = gi.workflows.show_workflow(workflow_id)
    for step_id in wf_steps['steps']:
        for param_name in wf_steps['steps'][step_id]['input_steps']:
            if 'input_type_selector' in param_name:
                wf_input_id = step_id
                tool_id = wf_steps['steps'][step_id]['tool_id']

    f1=open(options.output_batch, 'w+')
    f1.write("### METADATA\n")
    f1.write("#######################################\n")
    f1.write("Workflow Name\t%s\n" % options.workflow_name)
    f1.write("Workflow id\t%s\n" % workflow_id)
    f1.write("Project Name\t%s\n" % project_name)
    f1.write("#######################################\n\n")
    f1.write("###TABLE DATA\n")
    f1.write("#######################################\n")
    f1.write("SampleName\t##SourceType::SourceName::InputDir\t##Param::%s::%s::input_type_selector::input_sample_name\n" % (wf_input_id, tool_id))

    # walk through directory. Sample name is directory name, Parameter is the path to directory
    history = gi.histories.get_histories(history_id=options.history_id)[0]
    #print history
    history_name = history["name"]
    source = "history::%s::%s" % (history_name, options.history_item)
    for item in os.listdir(options.input_dir):
        if item == "tmp-illumina-dir":
            continue
        # is it a directory:
        if os.path.isdir(options.input_dir + "/" + item):
            f1.write("%s\t%s\t%s\n" % (item, source, item) )

    f1.close()

if __name__=="__main__": __main__()

