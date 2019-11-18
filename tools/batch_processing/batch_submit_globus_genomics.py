#!/usr/bin/env python

import os, sys, csv, time, optparse, re, random
import string, socket, errno
from bioblend.galaxy import GalaxyInstance
from itertools import *
sys.path.insert( 0, os.path.dirname( __file__ ) )
#import globusonline.transfer.api_client as transfer_api
import globus_sdk
import requests
requests.packages.urllib3.disable_warnings()

#from common import get, submit, display

def random_number(maxi):
       return random.randint(0,maxi)

def get_user(gi):
    return gi.users.get_current_user()['username']

def _count_globus_queue_transfers(gi, token):
    #username = get_user(gi)
    #base_url=transfer_api.DEFAULT_BASE_URL
    #apikwargs = {'base_url':  base_url, 'goauth' : token}
    api = globus_sdk.TransferClient(authorizer=globus_sdk.AccessTokenAuthorizer(token))
    #api = transfer_api.TransferAPIClient(username, **apikwargs)
    q_count = None
    while q_count is None:
        try:
            tasks = api.task_list(num_results=50)
            q_count = 0
            for task in tasks:
                if task['status'] == "ACTIVE" or task['status'] == "INACTIVE":
                    q_count += 1
            return q_count
        except:
            time.sleep(random_number(60))
    return q_count

def get_user_running_history_count (gi):
    count = 0
    for history in gi.histories.get_histories():
        if gi.histories.get_status(history['id'])['state'] == 'running' or gi.histories.get_status(history['id'])['state'] == 'queued':
            count += 1
    return count

def get_user_threshold_parallel_jobs(t_file):
    fh = open(t_file, "r")
    line = fh.readline().rstrip("\n")
    fh.close()
    return line

def can_submit(t_file, gi, token):
    transfer_count = 0
    if len(token) > 0 and token is not None and token != "None":
        transfer_count = _count_globus_queue_transfers(gi, token)
    running_jobs = get_user_running_history_count(gi)
    threshold_jobs = get_user_threshold_parallel_jobs(t_file)

    if int(running_jobs) > int(threshold_jobs) or transfer_count >= 3:
        return 0
    else:
        return 1

def create_history_name (workflow_name, sampleName, project_name):
    current_time = time.strftime("%a_%b_%d_%Y_%-I:%M:%S_%p",time.localtime(time.time()))
    history_name = "%s~%s~%s~%s" % (project_name, sampleName, workflow_name, current_time)
    return history_name

def get_src_id (mySrcName, srcItems):
    src_id = None
    for item in srcItems:
        if item['name'] == mySrcName and item['deleted'] is False:
            src_id = item['id']
            return item['id']
    return None

def get_workflow (gi, workflow_name):
    # get id for your desired workflow

    workflows = gi.workflows.get_workflows()

    for __workflow__ in workflows:
        if __workflow__['name'] == workflow_name:
            workflow_id = __workflow__['id']
            break

    workflow = gi.workflows.show_workflow(workflow_id)
    if not workflow:
        if return_formatted == True:
            print "Workflow %s not found, terminating." % workflow
            sys.exit(1)
    else:
        return workflow

def get_file_id (myFileNameList, src, src_id, gi):
    myFileName = myFileNameList[-1]
    file_id = None
    if src == 'ld':
        contents = gi.libraries.show_library(src_id, contents=True)
        for f in contents:
            if f['type'] == 'file' and myFileName in f['name']:
                file_id = f['id']
                return f['id']

    else:
        contents = gi.histories.show_history(src_id, contents=True)
        for f in contents:
            if f['type'] == 'file' and myFileName in f['name'] and gi.histories.show_dataset(src_id, f['id'])['deleted'] is False:
                file_id = f['id']
                return f['id']

    return None

def isBadLine(line):
    return line.startswith('#')

# usage: sys.argv[0] -k <api_key> -u <galaxy_url> -i <file_with_sample_info_and_parameters>
# example:  python submit_batch_workflow.py -k 32a9f8d74b8a3920ea2e89ed1ebeb27b -u http://localhost:8080 -i ~/mytest_sample_test.txt
################################################################################

#Parse Command Line
#api_url = "http://cvrg.galaxycloud.org/api/"

parser = optparse.OptionParser()
parser.add_option( '-k', '--api-key', dest='api_key', help='Mandatory. Need to get your galaxy API key by visiting the Galaxy webpage http://dev.globusgenomics.edu' )
parser.add_option( '-u', '--url', dest='url_name', help='Galaxy URL' )
parser.add_option( '-i', '--input', dest='inputBatch', help='Table containing the samples you wish to run through workflow.' )
parser.add_option( '-t', '--token', dest='goauth', help='GO authentication token for GO tools' )
(options, args) = parser.parse_args()

try:
    options.api_key
    options.inputBatch
    options.url_name
except NameError:
    print 'usage: %s --token <GO_auth_token> --api-key <api_key> --url <galaxy_url> --input <file_with_read_information>' % os.path.basename( sys.argv[0] )
    sys.exit(1)

if options.api_key == "None":
    print 'ERROR: User does not have a valid API Key. Please create one by clicking on the User Menu directory'
    sys.exit(1)

if "http:" in options.url_name:
   options.url_name = options.url_name.replace("http", "https")

batchRun = {'samples' : [] }

# Make api connection with Galaxy api server
gi = GalaxyInstance(url=options.url_name, key=options.api_key)

# Parse the input file
table_data_lines = []
table_data_flag = 0
header_flag = 0
header_column_count = 1
with open(options.inputBatch, 'rU') as input_file:
    lines = (line.rstrip("\n") for line in input_file)
    lines = (line for line in lines if line)
    for line in dropwhile(isBadLine, lines):
        if line.startswith('Workflow Name'):
            key, value = line.split("\t")
            batchRun['workflow_name'] = value
            print line
        elif line.startswith('Workflow id'):
            key, value = line.split("\t")
            batchRun['workflow_id'] = value
            print line
        elif line.startswith('Project Name'):
            key, value = line.split("\t")
            batchRun['project_name'] = value
            print line
        elif line.startswith('###TABLE DATA'):
            table_data_flag = 1
        elif table_data_flag == 1 and not line.startswith('####') and header_flag == 0:
            table_data_lines.append(line)
            print line
            header = line
            for col in header.split("\t"):
                if "##" in col:
                    header_column_count += 1
            header_flag = 1
        elif table_data_flag == 1 and not line.startswith('####') and header_flag == 1:
            line_column_count = len(line.split("\t"))
            if header_column_count != line_column_count:
                print 'ERROR: You are missing a parameter in your table for sample:\n%s\n%s' % (table_data_lines[0], line)
                sys.exit(1)
            table_data_lines.append(line)
            print line

print "STATUS\tSampleName\tWorkflowName\tWorkflowID\tHistoryName\tHistoryID"

t_file = "/opt/galaxy/tools/batch_processing/parallel_count_threshold.txt"
reader = csv.DictReader(table_data_lines, delimiter='\t')
for line in reader:
    ds_map = {}
    parameters = {}
    sampleName = None

    for headerSpot, value in line.iteritems():
        if "##" in headerSpot:
            annotation, key = headerSpot.split("##")
        else:
            key = headerSpot
        if key.startswith('SourceType::'):
            subkeys = key.split("::")
            key_values = value.split("::")
            if len(key_values) < 2:
                print 'ERROR: Your SourceType items should be separated by double colon \"::\" characters not single colon \":\" character.'
                sys.exit(1)

            #ds_map.update({subkeys[2]: { 'srcType': key_values[0], 'srcName': key_values[1], 'fileName': key_values[2] }})
            ds_map.update({subkeys[2]: { 'srcType': key_values[0], 'srcName': key_values[1], 'fileName': key_values[2:] }})
            if ds_map[subkeys[2]]['srcType'] == 'library':
                ds_map[subkeys[2]]['srcType'] = 'ld'
            elif ds_map[subkeys[2]]['srcType'] == 'history':
                ds_map[subkeys[2]]['srcType'] = 'hda'

            # get id from library
            if ds_map[subkeys[2]]['srcType'] == "ld":
                libs = gi.libraries.get_libraries()
                src_id = get_src_id(ds_map[subkeys[2]]['srcName'], libs)
                if src_id is None:
                    print >> sys.stderr, 'ERROR: Your library name does not exist:\n%s\n%s\n%s Does not exist' % (table_data_lines[0], line, ds_map[subkeys[2]]['srcName'])
                    sys.exit(1)
                file_id = get_file_id (ds_map[subkeys[2]]['fileName'], ds_map[subkeys[2]]['srcType'], src_id, gi)
                if file_id is None:
                    print >> sys.stderr, 'ERROR: Your library file does not exist:\n%s\n%s\n%s Does not exist' % (table_data_lines[0], line, ds_map[subkeys[2]]['fileName'])
                    sys.exit(1)
                ds_map[subkeys[2]].update({'src_id':src_id, 'file_id':file_id})

            # get id from history
            elif ds_map[subkeys[2]]['srcType'] == "hda":
                histories = gi.histories.get_histories()
                src_id = get_src_id(ds_map[subkeys[2]]['srcName'], histories)
                if src_id is None:
                    print >> sys.stderr, 'ERROR: Your History does not exist:\n%s\n%s\n%s Does not exist' % (table_data_lines[0], line, ds_map[subkeys[2]]['srcName'])
                    sys.exit(1)
                #print src_id
                file_id = get_file_id (ds_map[subkeys[2]]['fileName'], ds_map[subkeys[2]]['srcType'], src_id, gi)
                if file_id is None:
                    print >> sys.stderr, 'ERROR: Your History file does not exist:\n%s\n%s\n%s Does not exist' % (table_data_lines[0], line, ds_map[subkeys[2]]['fileName'])
                    sys.exit(1)
                #print file_id
                ds_map[subkeys[2]].update({'src_id':src_id, 'file_id':file_id})

        elif key.startswith('Param::'):
            param_values = key.split("::")
            key_type = param_values.pop(0)
            key_toolid = param_values.pop(0)
            key_toolName = param_values.pop(0)
            #key_type, key_toolid, key_toolName, key_toolParam = key.split("::")
            key_parameter = "|".join(param_values)

            if str(key_toolid) in parameters.keys():
                parameters[str(key_toolid)][key_parameter] =  value
            else:
                parameters[str(key_toolid)] = {key_parameter: value}

            if "globus" in key_toolName:
                parameters[str(key_toolid)]['goauth'] = options.goauth
                parameters[str(key_toolid)]['username'] = get_user(gi)

        elif key.startswith('SampleName'):
            sampleName = value

    # create history, move file, start workflow:

    # associate the file name with the file_id

    workflow = get_workflow(gi, batchRun['workflow_name'])
    #workflow = gi.workflows.show_workflow(batchRun['workflow_id'])

    # create history
    history_name = create_history_name(batchRun['workflow_name'], sampleName, batchRun['project_name'])
    history = gi.histories.create_history(name=history_name)
    data = {}
    data[ 'name' ] = history_name

    for step_id, ds_in in workflow['steps'].iteritems():
        if ds_in['tool_id'] == 'delete_history' or ds_in['tool_id'] == 'add_to_library':
            parameters[str(step_id)] = {'historyid': history['id'] }

    # submit workflow with files in library into new history
    wf_data = {}
    wf_data['workflow_id'] = workflow['id']
    wf_data['history'] = history_name
    wf_data['ds_map'] = {}
    wf_data['parameters'] = parameters
    for step_id, ds_in in workflow['inputs'].iteritems():
        #print ds_in
        #print step_id
        ds_in['label']
        wf_data['ds_map'][step_id] = {'src':ds_map[ds_in['label']]['srcType'], 'id':ds_map[ds_in['label']]['file_id']}


        #print "      Step Input ID: %s" % step_id
        #print "         Input label: %s" % ds_in['label']

    try:
        # check if user has reached its parallel job limit
        while can_submit(t_file, gi, options.goauth) == 0:
            time.sleep(60)

        gi_connect = GalaxyInstance(url=options.url_name, key=options.api_key)
        #res = gi_connect.workflows.run_workflow(wf_data['workflow_id'], wf_data['ds_map'], params=parameters, history_id=history['id'], import_inputs_to_history=True)
        res = gi_connect.workflows.invoke_workflow(wf_data['workflow_id'], inputs=wf_data['ds_map'], params=parameters, history_id=history['id'], import_inputs_to_history=False)

        print "SUBMITTED\t%s\t%s\t%s\t%s\t%s" % (sampleName, workflow['name'], workflow['id'], history_name, res['history_id'])
    except Exception, e:
        print "There was an error: %s" % str(e)
        sys.exit(1)

#    if res:
#        # you can loop through items in res to view the output files
#        #print res
#        print "   History Name: %s" % history_name
#        #print "   ID: %s" % res['history']
#    print "Submitted Workflow history: %s" %  history_name
    time.sleep(10)


