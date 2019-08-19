import sys, os, shutil, urllib, zipfile
from minid_client import minid_client_api
from bioblend.galaxy import GalaxyInstance
from itertools import *
import errno, csv, urllib2, json, time, optparse, tempfile, subprocess, psycopg2, socket, re
from binascii import unhexlify
from string import Template
import requests
requests.packages.urllib3.disable_warnings()

def cleanup(name, dirname):
    os.remove(name)
    shutil.rmtree(dirname)

def get_minid_data(minid, outmeta):
    entities = minid_client_api.get_entities("http://minid.bd2k.org/minid", minid, False)
    locations = entities[minid]['locations']

    if len(locations) < 1:
        print "No locations found for minid %s"  % sys.argv[1]
           
    name = locations[0]["uri"].rsplit('/', 1)[1]
    testfile = urllib.URLopener()
    testfile.retrieve(locations[0]["uri"], name)

    zip = zipfile.ZipFile(name)
    zip.extractall()

    # search for the file in the data directory
    dirname = name.split(".")[0:-1][0]
    metafile = "%s/data/metadata.tsv" % dirname

    if os.path.exists(metafile):
        shutil.copyfile(metafile, outmeta)
        cleanup(name, dirname)
    else:
        cleanup(name, dirname)
        sys.exit(1)


def __main__():
    print "Start time:"
    print time.strftime('%d/%m/%Y %H:%M:%S\n', time.localtime(time.time()))

    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '-m', '--minid', dest='minid', action='store', type="string", help='Minid containing tissue data' )
    parser.add_option( '--output-metadata', dest='outmeta', action='store', type="string", help='output metadata' )
    parser.add_option( '-w', '--workflow-name', dest='workflow_name', action='store', type="string", help='Workflow name' )
    parser.add_option( '-I', '--input-dataset', dest='inputDatasets', action='append', type="string", nargs=5, help='"datasource" "sourceName" "file_name" "workflow input label" "input_directory"' )
    parser.add_option( '-P', '--input-parameters', dest='inputParameters', action='append', type="string", nargs=4, help='"param_id" "param_toolname" "param_name" "param_value"' )
    parser.add_option( '-o', '--output', dest='output_batch', action='store', type="string", help='output_batch' )
    parser.add_option( '-k', '--key', dest='api_key', help='Mandatory. Need to get your galaxy API key by visiting the Galaxy webpage http://dev.globusgenomics.edu' )
    parser.add_option( '--url', dest='url_name', help='Galaxy URL' )
    parser.add_option( '--username', dest='username', help='The user name from the postgresql table associated with the user' )
    parser.add_option( '-g', '--group', dest="group", help="Name of group in the bag to analyze" )
    (options, args) = parser.parse_args()

    # get the minid data
    get_minid_data(options.minid, options.outmeta)
 
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
    header = "SampleName"

    ## Add input datasets from the workflow to the header
    ordered_inputs = []
    wf_steps = gi.workflows.show_workflow(workflow_id)
    inputLabels = {}
    for input_id in wf_steps['inputs']:
        header += "\t##SourceType::SourceName::%s" % wf_steps['inputs'][input_id]['label']
        inputLabels[wf_steps['inputs'][input_id]['label']] = ""
        ordered_inputs.append(wf_steps['inputs'][input_id]['label'])

    ## Add parameter inputs to the header
    ## Figure out which are the Runtime parameters from the json object of the workflow
    tool_dict = {}
    for step_id in wf_steps['steps']:
        if wf_steps['steps'][step_id]['tool_id'] in tool_dict:
            tool_dict[wf_steps['steps'][step_id]['tool_id']].append(int(step_id))
            tool_dict[wf_steps['steps'][step_id]['tool_id']] = sorted(tool_dict[wf_steps['steps'][step_id]['tool_id']])
        else:
            tool_dict[wf_steps['steps'][step_id]['tool_id']] = [int(step_id)]

    wf_json = gi.workflows.export_workflow_json(workflow_id)
    runtime_dict = {"__class__": "RuntimeValue"}

    ## map the ids from the wf_json dict to wf_steps dict
    step_maps = {}
    sorted_steps = {int(k) for k,v in wf_json['steps'].items()}
    for step_id in sorted_steps:
        tool_id =  wf_json['steps'][str(step_id)]['tool_id']
        step_maps[step_id] = tool_dict[tool_id].pop(0)

    params_header = {}
    for step_id in wf_json['steps']:
        tool_id = wf_json['steps'][step_id]['tool_id']
        if "RuntimeValue" in wf_json['steps'][step_id]['tool_state']:
            sub_dict = json.loads(wf_json['steps'][step_id]['tool_state'])
            for param in sub_dict:
                if type(sub_dict[param]) is unicode and "RuntimeValue" in sub_dict[param]:
                    sub_sub_dict = json.loads(sub_dict[param])
                    if sub_sub_dict == runtime_dict:
                        params_header[str(step_maps[int(step_id)])] = "##Param::%s::%s::%s" % (step_maps[int(step_id)], tool_id, param) 
                    else:
                        for sub_param in sub_sub_dict:
                            if type(sub_sub_dict[sub_param]) is dict:
                                if sub_sub_dict[sub_param] == runtime_dict:
                                    params_header[str(step_maps[int(step_id)])] = "##Param::%s::%s::%s::%s" % (step_maps[int(step_id)], tool_id, param, sub_param)

    ordered_samples = params_header.values()
    header += "\t" + "\t".join(params_header.values())
    f1=open(options.output_batch, 'w+')
    f1.write("### METADATA\n")
    f1.write("#######################################\n")
    f1.write("Workflow Name\t%s\n" % options.workflow_name)
    f1.write("Workflow id\t%s\n" % workflow_id)
    f1.write("Project Name\t%s\n" % project_name)
    f1.write("#######################################\n\n")
    f1.write("###TABLE DATA\n")
    f1.write("#######################################\n")
    f1.write("%s\n" % header)

    # get any static value parameters
    inputPs = {}
    if options.inputParameters:
        for (param_id, param_toolname, param_name, param_value) in options.inputParameters:
            id1 = "%s::%s::%s" % (param_id, param_toolname, param_name)
            inputPs[id1] = param_value

    # walk through directory. Sample name is directory name, Parameter is the path to directory
    sample_names = {}
    if options.minid:
        for param in sorted(params_header.values()):
            values = param.split("::")
            tool_id = values[1]
            tool_name = values[2]
            param_name = values[-1]
            id1 = "%s::%s::%s" % (tool_id, tool_name, param_name)
            #fh = open(options.outmeta, "r")
            #header = fh.readline()
            fh = ["GM06990", "GM10248", "GM12864", "GM12865", "GM12878", "GM12891", "GM12892", "GM13976", "GM18507","GM19238", "GM19239", "GM19240" ]
            for line in fh:
                item = line.split("\t")[0]
                if id1 in inputPs:
                    if item in sample_names:
                        sample_names[item].update( {param : inputPs[id1]} )
                    else:
                        sample_names[item] = { param : inputPs[id1] }
                else:
                    if item in sample_names:
                        sample_names[item].update( {param : item} )
                    else:
                        sample_names[item] = { param : item }

    if options.inputDatasets:
        for (sourceType, sourceName, fileName, label, inputDirectory) in options.inputDatasets:
            if label in inputLabels:
                if sourceType == "history":
                    history = gi.histories.get_histories(history_id=sourceName)[0]
                    history_name = history["name"]
                    inputLabels[label] = "history::%s::%s" % (history_name, fileName)

                    ## If there are "Param" type parameters in the workflow then loop through items in the inputDirectory
                    ## and use the sample names as parameter values
                    if len(params_header.values()) > 0:
                        for param in sorted(params_header.values()):
                            values = param.split("::")
                            tool_id = values[1]
                            tool_name = values[2]
                            param_name = values[-1]
                            id1 = "%s::%s::%s" % (tool_id, tool_name, param_name)
                            for item in os.listdir(inputDirectory):
                                if id1 in inputPs:
                                    if item in sample_names:
                                        sample_names[item].update( {param : inputPs[id1]} )
                                    else:
                                        sample_names[item] = { param : inputPs[id1] }
                                else:
                                    full_name = None
                                    if options.group != "None":
                                        full_name = "%s/%s" % (options.group, item)
  
                                    if os.path.isdir(inputDirectory + "/" + item):
                                        if item in sample_names:
                                            sample_names[item].update( {param : full_name} )
                                        else:
                                            sample_names[item] = { param : full_name }

                elif sourceType == "library":
                    inputLabels[label] = "library::%s::%s" % (sourceName, fileName)
            else:
                sys.exit("Label %s does not exist in workflow for Input Datasets. Make sure you have the correct label." % label)


    ## fill out the table with a row for each sample
    if len(sample_names) == 0:
        sample_names["Run1"] = {}

    for sample in sample_names:
        row = [sample]
        for label in ordered_inputs:
            row.append(inputLabels[label])
        for column in ordered_samples:
            row.append(sample_names[sample][column])
        f1.write("%s\n" % "\t".join(row) )

    f1.close()

if __name__=="__main__": __main__()

