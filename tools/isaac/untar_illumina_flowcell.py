#!/usr/bin/env python

"""
A wrapper script for untarring illumina flowcell tar data 
and creating a flowcell library and adding the RunInfo.xml file to the new library
"""

from bioblend.galaxy import GalaxyInstance
from itertools import *
import errno, csv, urllib2, json, time, sys, optparse, os, tempfile, subprocess, shutil, psycopg2, socket, re
from binascii import unhexlify
from string import Template
import requests
requests.packages.urllib3.disable_warnings()
#sys.path.insert( 0, "/nfs/software/galaxy/scripts/api" )
#from common import submit, display

def __main__():
    print "Start time:"
    print time.strftime('%d/%m/%Y %H:%M:%S\n', time.localtime(time.time()))

    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '-f', '--flowcell_name', dest='flowcell_name', action='store', type="string", help='Flowcell sample name' )
    parser.add_option( '-t', '--tar_file', dest='tar_file', action='store', type="string", help='Illumina flowcell tar file' )
    parser.add_option( '', '--output_dir', dest='output_dir', action='store', type="string", default=None, help='output_directory' )
    parser.add_option( '', '--output_config', dest='output_config', action='store', type="string", default=None, help='output config file' )
    parser.add_option( '--user', dest='userid', help='The user id from the postgresql table associated with the user' )
    parser.add_option( '--username', dest='username', help='The user name from the postgresql table associated with the user' )
    parser.add_option( '--rootdir', dest='rootdir', help='The root directory for galaxy installation' )
    parser.add_option( '-k', '--key', dest='api_key', help='Mandatory. Need to get your galaxy API key by visiting the Galaxy webpage http://dev.globusgenomics.edu' )
    parser.add_option( '--url', dest='url_name', help='Galaxy URL' )
    (options, args) = parser.parse_args()

    # make tmpdir and output dir
    if not os.path.isdir(options.output_dir):
        os.mkdir(options.output_dir)

    cmd = "tar xvf %s -C %s" % ( options.tar_file, options.output_dir )

    #set up stdout and stderr output options
    stdout_name = tempfile.NamedTemporaryFile( prefix = "stdout" ).name
    stderr_name = tempfile.NamedTemporaryFile( prefix = "stderr" ).name
   
    print cmd
    buffsize = 1048576
 
    proc = subprocess.Popen( args=cmd, stdout=open( stdout_name, 'wb' ), stderr=open( stderr_name, 'wb' ), shell=True )

    return_code = proc.wait()
    if return_code:
        tmp_stderr = open( stderr_name, 'rb' )
        stderr = ''
        try:
            while True:
                stderr += tmp_stderr.read( buffsize )
                if not stderr or len( stderr ) % buffsize != 0:
                    break
        except OverflowError:
            pass
        tmp_stderr.close()
        if return_code != 0:
            raise Exception, stderr

    # create the configuration output file with location for the RunInfo.xml, Flowcell Name and any other information that could be useful downstream
    xml_location = None
    for root, dirs, files in os.walk(options.output_dir):
        if 'RunInfo.xml' in files:
            #xml_location = os.path.join(root, 'RunInfo.xml')
            xml_location = root
            break
    print "RunInfo.xml is at: %s" % xml_location
    #f1=open(options.output_config, 'w+')
    #f1.write('FLOWCELL_NAME\t%s\n' % options.flowcell_name)
    #f1.write('RUNINFO\t%s\n' % xml_location)
    #f1.close()

    # get api_url
    api_url = options.url_name
    if "http:" in api_url:
        api_url = api_url.replace("http", "https")

    # get user's API key
    api_key = options.api_key

    # connect to Galaxy vi API
    gi = GalaxyInstance(url=api_url, key=api_key)

    # create shared library if necessary
    data_library = options.flowcell_name
    libs = gi.libraries.get_libraries()
    library_id = None
    for library in libs:
        if library['name'] == data_library:
            library_id = library['id']
    if not library_id:
        library = gi.libraries.create_library(data_library)
        #lib_create_data = {'name':data_library}
        print library
        library_id = library['id']
    folders = gi.libraries.get_folders(library_id, name='/')
    for f in folders:
        if f['name'] == "/":
            library_folder_id = f['id']
    if not library_id or not library_folder_id:
        print "Failure to configure library destination."
        sys.exit(1)

    # only add RunInfo.xml file to library
    data = {}
    data['folder_id'] = library_folder_id
    data['file_type'] = 'xml'
    data['dbkey'] = 'hg19'
    data['upload_option'] = 'upload_paths'
    data['filesystem_paths'] = xml_location
    data['create_type'] = 'file'
    data['link_data_only'] = 'link_to_files'

    libset = gi.libraries.upload_file_from_server(library_id, xml_location, folder_id=library_folder_id, file_type='xml', dbkey='hg19', link_data_only='link_to_files' )
    #libset = gi.libraries.upload_file_from_local_path(library_id, xml_location, folder_id=library_folder_id, file_type='xml', dbkey='hg19' )
    print libset

    ##########
    # create the batch submit file for the workflow (isaac-align-bcl-qc)
    workflow_name = "Optimized-Isaac-align-per-lane"
    wf = gi.workflows.get_workflows(name=workflow_name)

    #workflow_id = "1234567890abcd"
    workflow_id = wf[0]['id']
    project_name = options.flowcell_name

    wf_steps = gi.workflows.show_workflow(workflow_id)
    for step_id in wf_steps['steps']:
        for param_name in wf_steps['steps'][step_id]['input_steps']:
            if 'input_type_selector' in param_name:
                wf_input_id = step_id
                tool_id = wf_steps['steps'][step_id]['tool_id']

    f1=open(options.output_config, 'w+')
    f1.write("### METADATA\n")
    f1.write("#######################################\n")
    f1.write("Workflow Name\t%s\n" % workflow_name)
    f1.write("Workflow id\t%s\n" % workflow_id)
    f1.write("Project Name\t%s\n" % project_name)
    f1.write("#######################################\n\n")
    f1.write("###TABLE DATA\n")
    f1.write("#######################################\n")

    #if options.username == "mfabani":
    #    f1.write("SampleName\tSourceType::SourceName::RunInfo.xml\tParam::1768::isaac_align::input_type_selector::lane\n") ## hardcoded for Martin's workflow
    #else:
    #    f1.write("SampleName\tSourceType::SourceName::RunInfo.xml\tParam::1768::isaac_align::input_type_selector::lane\n")  ## hardcoded for Alex workflow
    f1.write("SampleName\tSourceType::SourceName::RunInfo.xml\tParam::%s::%s::input_type_selector::lane\n" % (wf_input_id, tool_id)) 
    f1.write("L001\tlibrary::%s::RunInfo.xml\tL001\n" % data_library)
    f1.write("L002\tlibrary::%s::RunInfo.xml\tL002\n" % data_library)
    f1.write("L003\tlibrary::%s::RunInfo.xml\tL003\n" % data_library)
    f1.write("L004\tlibrary::%s::RunInfo.xml\tL004\n" % data_library)
    f1.write("L005\tlibrary::%s::RunInfo.xml\tL005\n" % data_library)
    f1.write("L006\tlibrary::%s::RunInfo.xml\tL006\n" % data_library)
    f1.write("L007\tlibrary::%s::RunInfo.xml\tL007\n" % data_library)
    f1.write("L008\tlibrary::%s::RunInfo.xml\tL008\n" % data_library)
    f1.close()

    print "\n\nEnd time:"
    print time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))

if __name__=="__main__": __main__()
