# runs after the job (and after the default post-filter)
from galaxy import datatypes, jobs

def validate(incoming):
    """Validator"""
    #raise Exception, 'not quite right'
    pass

"""

def exec_before_job( app, inp_data, out_data, param_dict, tool=None):
    #Sets the name of the data
    import subprocess, os

    accessToken = param_dict.get( 'accessToken', None )
    fileId = param_dict.get( 'dataID', None )
    cmd = "curl -L -H \"x-access-token: %s\" https://api.basespace.illumina.com/v1pre3/files/%s | jq \'.Response.Name\' -r" % ( accessToken, fileId )
    print "CMD: %s" % cmd
    process = subprocess.Popen(args=cmd, stdout=subprocess.PIPE, stderr=None, shell=True)
    output = process.communicate()
    print "OUTPUT: %s" % str(output)
    filename = output[0].rstrip()
    print "FILE OUT: %s" % filename    
    #filename = filename.replace("\"", "")
    ext = os.path.basename(filename).split(".")[-1]

    #outputType = param_dict.get( 'hgta_outputType', None )
    #if isinstance(outputType, list) and len(outputType)>0: outputType = outputType[-1]
    items = out_data.items()
    print items
    for name, data in items:
        print data.name
        data.name  = filename
        print "FILENAME %s" % filename
        print "AFTER: %s" % data.name

        data = app.datatypes_registry.change_datatype(data, ext)
        out_data[name] = data
"""
        
def exec_after_process( app, inp_data, out_data, param_dict, tool, stdout, stderr):
    #"Verifies the data after the run"
    import os

    items = out_data.items()
    filename = stdout.rstrip()
    #print "FILE OUT: %s" % filename
    #filename = filename.replace("\"", "")
    ext = os.path.basename(filename).split(".")[-1]


    for name, data in items:
        #print data.name
        data.name  = filename
        #print "FILENAME %s" % filename
        #print "AFTER: %s" % data.name

        data = app.datatypes_registry.change_datatype(data, ext)
        out_data[name] = data



