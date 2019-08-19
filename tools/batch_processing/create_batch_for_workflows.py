#!/usr/bin/env python

from bioblend.galaxy import GalaxyInstance
from itertools import *
import errno, csv, urllib2, json, time, sys, optparse, os, tempfile, subprocess, shutil, psycopg2, socket, re
from binascii import unhexlify
from string import Template
import requests
requests.packages.urllib3.disable_warnings()

def _get_samples(inputDirectory,datatype):
    # assume paired end fastq files as input
    if datatype != "bam":
        forward = []
        reverse = []
        count = 0
        for infile in sorted(os.listdir(inputDirectory)):
            if infile.endswith(".fastq.gz"):
                if "_R1_" in infile or "_1_" in infile or "_1.fastq.gz" in infile:
                    forward.append(infile)
                elif "_R2_" in infile or "_2_" in infile or "_2.fastq.gz" in infile:
                    reverse.append(infile)
        sample_pairs = list(zip(forward, reverse))
        sample_names = []
        for pair in sample_pairs:
            sample_names.append(_get_sampleName(pair[0], pair[1]))
        return sample_names
    else:
        bams = []
        for infile in sorted(os.listdir(inputDirectory)):
            print infile
            if infile.lower().endswith("bam"):
                print infile
                bams.append(os.path.basename(infile).split(".")[0])
        print bams
        return bams

def _get_sampleName(forward, reverse):
    common = []
    count = 0
    f_vals = list(forward)
    r_vals = list(reverse)
    if "_R1_" in forward:
        values = forward.split("_R1")
        return values[0]
    else:
        for i in f_vals:
            if f_vals[count] == r_vals[count]:
                common.append(f_vals[count])
            else:
                return "".join(common)
            count += 1

def __main__():
    print "Start time:"
    print time.strftime('%d/%m/%Y %H:%M:%S\n', time.localtime(time.time()))

    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '-w', '--workflow-name', dest='workflow_name', action='store', type="string", help='Workflow name' )
    parser.add_option( '-I', '--input-dataset', dest='inputDatasets', action='append', type="string", nargs=7, help='"datasource" "sourceName" "file_name" "workflow input label" "input_directory" "multigroup_samples" "datatype"' )
    parser.add_option( '-P', '--input-parameters', dest='inputParameters', action='append', type="string", nargs=4, help='"param_id" "param_toolname" "param_name" "param_value"' )
    parser.add_option( '-o', '--output', dest='output_batch', action='store', type="string", help='output_batch' )
    parser.add_option( '-k', '--key', dest='api_key', help='Mandatory. Need to get your galaxy API key by visiting the Galaxy webpage http://dev.globusgenomics.edu' )
    parser.add_option( '--url', dest='url_name', help='Galaxy URL' )
    parser.add_option( '--username', dest='username', help='The user name from the postgresql table associated with the user' )
    parser.add_option( '-g', '--group', dest="group", help="Name of group in the bag to analyze" )
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
    step_inputs = {}
    for step_id in wf_steps['steps']:
        if wf_steps['steps'][step_id]['tool_id'] in tool_dict:
            tool_dict[wf_steps['steps'][step_id]['tool_id']].append(int(step_id))
            tool_dict[wf_steps['steps'][step_id]['tool_id']] = sorted(tool_dict[wf_steps['steps'][step_id]['tool_id']])
        else:
            tool_dict[wf_steps['steps'][step_id]['tool_id']] = [int(step_id)]
        # get the input connections
        for key in wf_steps['steps'][step_id]['input_steps']:
            if step_id in step_inputs:
                step_inputs[step_id].append(key)
            else:
                step_inputs[step_id] = [key]

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
                if type(sub_dict[param]) is unicode and "RuntimeValue" in sub_dict[param] and param not in step_inputs[step_id] and "file" not in param:
                    sub_sub_dict = json.loads(sub_dict[param])
                    if sub_sub_dict == runtime_dict:
                        if str(step_maps[int(step_id)]) in params_header:
                            params_header[str(step_maps[int(step_id)])].append("##Param::%s::%s::%s" % (step_maps[int(step_id)], tool_id, param))
                        else:
                            params_header[str(step_maps[int(step_id)])] = [ "##Param::%s::%s::%s" % (step_maps[int(step_id)], tool_id, param) ]
                    else:
                        for sub_param in sub_sub_dict:
                            if type(sub_sub_dict[sub_param]) is dict:
                                if sub_sub_dict[sub_param] == runtime_dict and sub_param not in step_inputs[step_id] and "file" not in sub_param:
                                    if str(step_maps[int(step_id)]) in params_header:
                                        params_header[str(step_maps[int(step_id)])].append("##Param::%s::%s::%s::%s" % (step_maps[int(step_id)], tool_id, param, sub_param))
                                    else:
                                        params_header[str(step_maps[int(step_id)])] = [ "##Param::%s::%s::%s::%s" % (step_maps[int(step_id)], tool_id, param, sub_param) ]
                                else:
                                    sub_sub_sub_dict = sub_sub_dict[sub_param]
                                    for sub_sub_param in sub_sub_sub_dict:
                                        if sub_sub_sub_dict[sub_sub_param] == runtime_dict and sub_sub_param not in step_inputs[step_id] and "file" not in sub_sub_param:
                                            if str(step_maps[int(step_id)]) in params_header:
                                                params_header[str(step_maps[int(step_id)])].append("##Param::%s::%s::%s::%s::%s" % (step_maps[int(step_id)], tool_id, param, sub_param, sub_sub_param))
                                            else:
                                                params_header[str(step_maps[int(step_id)])] = [ "##Param::%s::%s::%s::%s::%s" % (step_maps[int(step_id)], tool_id, param, sub_param, sub_sub_param) ]
                                        else:
                                            sub_sub_sub_sub_dict = sub_sub_sub_dict[sub_sub_param]
                                            if type(sub_sub_sub_dict[sub_sub_param]) is dict:
                                                for sub_sub_sub_param in sub_sub_sub_sub_dict:
                                                    if sub_sub_sub_sub_dict[sub_sub_sub_param] == runtime_dict and sub_sub_sub_param not in step_inputs[step_id] and "file" not in sub_sub_sub_param:
                                                        if str(step_maps[int(step_id)]) in params_header:
                                                            params_header[str(step_maps[int(step_id)])].append("##Param::%s::%s::%s::%s::%s::%s" % (step_maps[int(step_id)], tool_id, param, sub_param, sub_sub_param, sub_sub_sub_param))
                                                        else:
                                                            params_header[str(step_maps[int(step_id)])] = [ "##Param::%s::%s::%s::%s::%s::%s" % (step_maps[int(step_id)], tool_id, param, sub_param, sub_sub_param, sub_sub_sub_param) ]
    ordered_samples = []
    for tool_params in params_header.values():
        ordered_samples.extend(tool_params)
        header += "\t" + "\t".join(tool_params)

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
    if options.inputDatasets:
        for (sourceType, sourceName, fileName, label, inputDirectory, multisample_dir, datatype) in options.inputDatasets:
            if label in inputLabels:
                if sourceType == "history":
                    history = gi.histories.get_histories(history_id=sourceName)[0]
                    history_name = history["name"]
                    inputLabels[label] = "history::%s::%s" % (history_name, fileName)

                    ## If there are "Param" type parameters in the workflow then loop through items in the inputDirectory
                    ## and use the sample names as parameter values
                    if len(params_header.values()) > 0:
                        for tool_params in sorted(params_header.values()):
                            for param in tool_params:
                                values = param.split("::")
                                tool_id = values[1]
                                tool_name = values[2]
                                #param_name = values[-1]
                                param_name = "::".join(values[3:])
                                id1 = "%s::%s::%s" % (tool_id, tool_name, param_name)
                                #print id1
                                if multisample_dir == "yes":
                                    iterator = multi_sampleNames = _get_samples(inputDirectory,datatype)
                                else:
                                    iterator = os.listdir(inputDirectory)
                                for item in iterator:
                                    if id1 in inputPs:
                                        if "GP::" in inputPs[id1]:
                                            if item in sample_names:
                                                sample_names[item].update( {param : item} )
                                            else:
                                                sample_names[item] = { param : item }
                                        else:
                                            if item in sample_names:
                                                sample_names[item].update( {param : inputPs[id1]} )
                                            else:
                                                sample_names[item] = { param : inputPs[id1] }
                                    else:
                                        full_name = None
                                        if options.group != "None":
                                            full_name = "%s/%s" % (options.group, item)
                                        else:
                                            full_name = "%s" % (item)

                                        if multisample_dir == "yes":
                                            if item in sample_names:
                                                sample_names[item].update( {param : full_name} )
                                            else:
                                                sample_names[item] = { param : full_name }
                                        else:
                                            if os.path.isdir(inputDirectory + "/" + item):
                                                if item in sample_names:
                                                    sample_names[item].update( {param : full_name} )
                                                else:
                                                    sample_names[item] = { param : full_name }

                elif sourceType == "library":
                    inputLabels[label] = "library::%s::%s" % (sourceName, fileName)
            else:
                sys.exit("Label %s does not exist in workflow for Input Datasets. Make sure you have the correct label." % label)


    print sample_names
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
