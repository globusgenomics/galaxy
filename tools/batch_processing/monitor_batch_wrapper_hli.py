#!/usr/bin/python
#!/usr/bin/python
# -*- coding: utf-8 -*-

import optparse, os, shutil, subprocess, sys, glob, tempfile, psycopg2, socket, re
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
parser.add_option( '--user', dest='userid', help='The user id from the postgresql table associated with the user' )
parser.add_option( '--input', dest='table_file', help='The table file containing the parameters for the workflow to submit' )
parser.add_option( '--out', dest='log', help='Log file' )
parser.add_option( '--out-dir', dest='output_dir', help='Output directory to store BAM, FASTQ, VCF and more' )
parser.add_option( '--rootdir', dest='rootdir', help='The root directory for galaxy installation' )
parser.add_option( '--url', dest='galaxyurl', help='Base URL' )
(options, args) = parser.parse_args()

# get the database, user, host, password, url
#print os.getcwd()
universe_file = options.rootdir + "/universe_wsgi_ci.ini"
ufile = open( universe_file, "r" )
for line in ufile:
    if "database_connection = " in line:
        matches = re.search(r"database_connection = postgres:\/\/\/(.*)\?user=(.*)\&password=(.*)\&host=(.*)", line)
        database = matches.group(1)
        user = matches.group(2)
        host = matches.group(4)
        password = matches.group(3)

    #elif "display_servers = " in line: 
        #matches = re.search(r"display_servers = (.*)", line)
        #url_host = matches.group(1)
    elif "globus_endpoint = " in line:
        matches = re.search(r"globus_endpoint = (.*)#(.*)", line)
        url_host = matches.group(2)
        
ufile.close()
#print "URL_UNIVERSE: %s" % url

# comment next line if socket.gethostname is not giving correct hostname
url_host = socket.gethostname()
#print "URL_HOST: %s" % url_host
url = "https://%s" % url_host

if not os.path.exists(options.output_dir):
    os.mkdir(options.output_dir)

con = None
try:
    con = psycopg2.connect(database=database, user=user, host=host, password=password) 
    cur = con.cursor()

    # get the user key
    #key_query_sql = "SELECT  api_keys.key FROM api_keys INNER JOIN galaxy_user ON galaxy_user.id=api_keys.user_id where galaxy_user.username='%s'" % (userid)
    #key_query_sql = "SELECT  api_keys.key FROM api_keys WHERE api_keys.user_id='%s'" % (options.userid)
    key_query_sql = "SELECT  api_keys.key FROM api_keys WHERE api_keys.user_id='%s' ORDER BY api_keys.create_time DESC LIMIT 1" % (options.userid)
    #print key_query_sql
    cur.execute(key_query_sql)
    key = cur.fetchone()

    try:
        if len(key) > 0:
           api_key = key[0]
           #print api_key
 
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
           summary_output = "%s/summary_output.txt" % options.output_dir
           fh_out = open(summary_output, "w")
           while len(monitor_meta) != len(completed_meta):
               for hKey in monitor_meta:
                   print "Analyzing HKEY %s" % hKey 
                   if hKey in completed_meta.keys():
                       continue
                   elif is_complete(hKey, gi) == True:
                       #fh_out.write("%s\n" % hKey)
                       completed_meta[hKey] = {'status' : 'complete', 'sampleName': monitor_meta[hKey]['sampleName'], 'datasets' : {}, 'metrics' : {}}

                       # get metadata about the completed jobs
                       for dataset in gi.histories.show_history(hKey, contents=True):
                           meta = gi.histories.show_dataset(hKey, dataset['id'])
                           #print meta
                           if meta['visible'] is True and meta['hid'] != 1 and meta['deleted'] is False:
                               completed_meta[hKey]['datasets'][dataset['id']] = { 'data_type' : meta['data_type'], 'file_name': meta['file_name'], 'job_name' : meta['name'].replace(" ", "_")}

           #print completed_meta
           # write the final output or store the VCF or BAM files into specific directory
           # parse the output data
           tools_metrics_list = {}
           job_names = {}
           for hKey in completed_meta:
               completed_meta[hKey]['raw'] = {}
               #print "DATSETS: %s" % completed_meta[hKey]['datasets']
               for dataset_id in completed_meta[hKey]['datasets']:
                   #print "GETTING READY TO ANALYZE: %s, %s" % (completed_meta[hKey]['sampleName'], completed_meta[hKey]['datasets'][dataset_id]['job_name'])
                   data_type = completed_meta[hKey]['datasets'][dataset_id]['data_type']
                   if data_type == 'vcf' or data_type == 'bam':
                       #print "NON_HTML: %s" % dataset_id
                       file_name = "%s.%s.%s" % (completed_meta[hKey]['sampleName'], completed_meta[hKey]['datasets'][dataset_id]['job_name'], completed_meta[hKey]['datasets'][dataset_id]['data_type'])
                       link_name = "%s/%s" % (completed_meta[hKey]['sampleName'], file_name)
                       destination_name = "%s/%s" % (options.output_dir, link_name)
                       #print "%s\t%s\t%s" % (data_type, destination_name, completed_meta[hKey]['datasets'][dataset_id]['file_name'])

                       # move the file to the destination
                       #shutil.move(completed_meta[hKey]['datasets'][dataset_id]['file_name'], destination_name)
                       if data_type.upper() in completed_meta[hKey]['raw']:
                           completed_meta[hKey]['raw'][data_type.upper()].append({'path' : completed_meta[hKey]['datasets'][dataset_id]['file_name'], 'destination' : destination_name, 'sample_name' : completed_meta[hKey]['sampleName'], 'data_type' : completed_meta[hKey]['datasets'][dataset_id]['data_type'], 'link_name' : link_name, 'filename' : data_type.upper()})
                       else:
                           completed_meta[hKey]['raw'][data_type.upper()] = [{'path' : completed_meta[hKey]['datasets'][dataset_id]['file_name'], 'destination' : destination_name, 'sample_name' : completed_meta[hKey]['sampleName'], 'data_type' : completed_meta[hKey]['datasets'][dataset_id]['data_type'], 'link_name' : link_name, 'filename' : data_type.upper()}]
                       job_names[data_type.upper()] = 0

                   elif completed_meta[hKey]['datasets'][dataset_id]['data_type'] == 'html':
                       # get path of metrics file
                       metrics_path = "%s/%s_files" % (os.path.dirname(completed_meta[hKey]['datasets'][dataset_id]['file_name']), ('.').join(os.path.basename(completed_meta[hKey]['datasets'][dataset_id]['file_name']).split('.')[:-1]))
                       metrics_file = None

                       for txt_file in glob.glob("%s/*.metrics.txt" %  metrics_path):
                           metrics_file = txt_file
                       completed_meta[hKey]['raw'][completed_meta[hKey]['datasets'][dataset_id]['job_name']] = []
                       job_names[completed_meta[hKey]['datasets'][dataset_id]['job_name']] = 0
                       for txt_file in glob.glob("%s/*" %  metrics_path):
                           data_type = os.path.basename(txt_file).split('.')[-1]
                           file_name = os.path.basename(txt_file)
                           link_name = "%s/%s" % (completed_meta[hKey]['sampleName'], file_name)
                           destination_name = "%s/%s" % (options.output_dir, link_name)
                           completed_meta[hKey]['raw'][completed_meta[hKey]['datasets'][dataset_id]['job_name']].append({'path' : txt_file, 'destination' : destination_name, 'sample_name' : completed_meta[hKey]['sampleName'], 'data_type' : data_type, 'link_name' : link_name, 'filename' : file_name})

                       #print "METRICS_FILE PATH: %s" % metrics_file 
                       # get the information from the txt file
                       # assume that it's a Picard metrics file and get the data
                       f_metrics = open(metrics_file, "r")
                       metrics_dict = {}
                       for metrics_line in f_metrics:
                           metrics_line = metrics_line.rstrip("\n")
                           if metrics_line.startswith("## METRICS CLASS"):
                               classy,tool_type = metrics_line.split("\t")
                               #print "METRICS_TOOLTYPE: %s" % tool_type
                               tool_values = metrics_line.split(".")
                               tool_name = tool_values[-1]
                               #print "METRICS TOOLVALUE: %s" % tool_name
                               metrics_header = f_metrics.next().rstrip("\n")
                               #print "METRICSHEADER: %s" % metrics_header
                               metrics_stats = f_metrics.next().rstrip("\n")
                               metrics_values = metrics_stats.split("\t")
                               #print "METRICS_VALUES: %s" % metrics_values
                               index = 0
                               for col in metrics_header.split("\t"):
                                   metrics_dict[col] = metrics_values[index]
                                   index += 1
                               if tool_name not in tools_metrics_list:
                                   tools_metrics_list[tool_name] = metrics_header.split("\t")
                               completed_meta[hKey]['metrics'][tool_name] = metrics_dict
                               #print "METRICS: %s: %s" % (tool_name, completed_meta[hKey]['metrics'][tool_name])
                               break

           # Print out final statistics to stdout or file
           # header
           fh_out.write("SAMPLE")
           for tool in tools_metrics_list:
               for column_name in tools_metrics_list[tool]:
                   fh_out.write("\t%s.%s" % (tool, column_name))
           fh_out.write("\n")

           # table content - the juicy details
           #print completed_meta
           for hKey in completed_meta:
               fh_out.write(completed_meta[hKey]['sampleName'])
               #print "PRINTING %s" % completed_meta[hKey]['sampleName']
               #print completed_meta[hKey]
               for tool in tools_metrics_list:
                   for column_name in tools_metrics_list[tool]:
                       if tool in completed_meta[hKey]['metrics']:
                           fh_out.write("\t%s" % completed_meta[hKey]['metrics'][tool][column_name])
                       else:
                           fh_out.write("\tMISSING")
               fh_out.write("\n")
           fh_out.close()
           html_output = options.log

           # write the final html output
           fh_html = open(html_output, "w")
           fh_html.write("<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">\n")
           fh_html.write("<html xmlns=\"http://www.w3.org/1999/xhtml\" xml:lang=\"en\" lang=\"en\">\n")
           fh_html.write("<head>\n")
           fh_html.write("  <meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\" />\n")
           fh_html.write("  <meta name=\"generator\" content=\"Galaxy picard_wrapper tool output - see http://getgalaxy.org/\" />\n")
           fh_html.write("  <title></title>\n")
           fh_html.write("  <link rel=\"stylesheet\" href=\"/static/style/base.css\" type=\"text/css\" />\n")
           fh_html.write("<style type=\"text/css\">\n")
           fh_html.write("tr.d0 td {background-color: oldlace; color: black;}\n")
           fh_html.write("tr.d1 td {background-color: aliceblue; color: black;}\n")
           fh_html.write("</style><?xml version=\"1.0\" encoding=\"utf-8\" ?>\n")       
           fh_html.write("</head>\n")

           fh_html.write("<body>\n")
           fh_html.write("  <div class=\"document\">\n")
           fh_html.write("    <b>Batch output summary</b><br/>The following output files were created (click the filename to view/download a copy):<hr/>\n")
           fh_html.write("    <table>\n")
           fh_html.write("      <tr><td><a href=\"%s\">%s</a></td></tr>\n" % ("summary_output.txt",  "Batch Summary Output") )
           fh_html.write("    </table>\n")
           fh_html.write("  </div>\n")

           fh_html.write("  <div class=\"document\" >\n")
           fh_html.write("    <b>Batch output per Sample</b><br/>The following output files were created for each sample (click the filename to view/download a copy):<hr/>\n")
           fh_html.write("    <table class=\"manage-table colored\">\n")
           fh_html.write("      <tr class=\"header\">\n")
           fh_html.write("        <td><b>Sample Name</b></td>\n")
           for job_name in job_names:
               fh_html.write("      <td><b>%s</b></td>\n" % job_name)
           fh_html.write("      </tr>\n") 

           for hKey in completed_meta:
               #print completed_meta[hKey]['raw']
               sample_name = completed_meta[hKey]['sampleName']
               destination_sample_dir = "%s/%s" % (options.output_dir, sample_name)
               if not os.path.exists(destination_sample_dir):
                   os.makedirs(destination_sample_dir)
                
               # generate a table for the sample
               fh_html.write("      <tr>\n")
               fh_html.write("        <td>%s</td>\n" % sample_name)

               # add column for each jobname
               for job_name in job_names:
                   fh_html.write("    <td>\n")
                   fh_html.write("      <table>\n")
                   #fh_html.write("        <tr><td>%s</td></tr>\n" % job_name)
                   for raw_dataset in completed_meta[hKey]['raw'][job_name]:
                       # copy or move the files to the new directory
                       shutil.copy(raw_dataset['path'], raw_dataset['destination'])
                       if raw_dataset['data_type'] == "jpg":
                           fh_html.write("      <tr><td><a href=\"%s\"><img align=\"middle\" hspace=\"1\" src=\"%s\" height=\"75\" width=\"75\" title=\"Click image preview for a print quality PDF version\" /></a></td></tr>\n" % (raw_dataset['link_name'], raw_dataset['link_name']) )
                       else:
                           fh_html.write("    <tr><td><a href=\"%s\">%s</a></td></tr>\n" % (raw_dataset['link_name'], raw_dataset['filename']))

                   fh_html.write("      </table>\n")
                   fh_html.write("    </td>\n")
               fh_html.write("    </tr>\n")

           fh_html.write("    </table>\n")
           fh_html.write("  </div></body></html>\n")

           # Now it is safe to delete the output histories
           for hKey in completed_meta:
               print "deleting history %s: %s" % (hKey, completed_meta[hKey]['sampleName'])
               gi.histories.delete_history(hKey, purge=True)

    except Exception as e:
        print "There is some kind of issue here: %s" % (e)
        sys.exit()

except psycopg2.DatabaseError, e:
    print 'Error %s' % e    


