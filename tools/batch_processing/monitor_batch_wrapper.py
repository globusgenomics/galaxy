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
parser.add_option( '--out', dest='log', help='Log file' )
parser.add_option( '--out-dir', dest='output_dir', help='Output directory to store BAM, FASTQ, VCF and more' )
(options, args) = parser.parse_args()

url = options.url_name
if "http:" in options.url_name:
   url = options.url_name.replace("http", "https")

# get the database, user, host, password, url
#print os.getcwd()

if not os.path.exists(options.output_dir):
    os.mkdir(options.output_dir)

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
                        print meta
                        if meta['visible'] is True and meta['deleted'] is False:
                            job_name = meta['name'].replace(" ", "_")
                            job_name = job_name.replace("/", "_")
                            job_name = job_name.replace(",", "_")
                            completed_meta[hKey]['datasets'][dataset['id']] = { 'data_type' : meta['data_type'], 'file_name': meta['file_name'], 'job_name' : job_name}
            if len(monitor_meta) != len(completed_meta):
                time.sleep(600)

        #print completed_meta
        # write the final output or store the VCF or BAM files into specific directory
        # parse the output data
        tools_metrics_list = {}
        job_names = {}
        for hKey, history_dict in completed_meta.iteritems():
            completed_meta[hKey]['raw'] = {}
            print "DATSETS: %s" % history_dict['datasets']
            for dataset_id, dataset_dict in history_dict['datasets'].iteritems():
                print "GETTING READY TO ANALYZE: %s, %s" % (history_dict['sampleName'], dataset_dict['job_name'])
                data_type = dataset_dict['data_type']
                if data_type != "html":
                    if data_type == "tabular" or data_type == "vcf":
                        file_name = "%s.%s.%s" % (history_dict['sampleName'], dataset_dict['job_name'], "txt")
                    else:
                        file_name = "%s.%s.%s" % (history_dict['sampleName'], dataset_dict['job_name'], data_type)
                    link_name = "%s/%s" % (history_dict['sampleName'], file_name)
                    destination_name = "%s/%s" % (options.output_dir, link_name)

                    # skip any files in history that were used as initial inputs to the workflow. 
                    # Usually any file with a non dat extension are provided as input
                    if dataset_dict['file_name'].split(".")[-1] != "dat":
                        continue

                    # move the file to the destination
                    #shutil.move(dataset_dict['file_name'], destination_name)
                    if data_type.upper() in history_dict['raw']:
                        completed_meta[hKey]['raw'][completed_meta[hKey]['datasets'][dataset_id]['job_name']].append({'path' : dataset_dict['file_name'], 'destination' : destination_name, 'sample_name' : history_dict['sampleName'], 'data_type' : data_type, 'link_name' : link_name, 'filename' : data_type.upper()})
                    else:
                        completed_meta[hKey]['raw'][completed_meta[hKey]['datasets'][dataset_id]['job_name']] = [{'path' : dataset_dict['file_name'], 'destination' : destination_name, 'sample_name' : history_dict['sampleName'], 'data_type' : data_type, 'link_name' : link_name, 'filename' : data_type.upper()}]

                    # if output is associated with extra_files_path then copy that too
                    extra_files_path = "%s/%s_files" % (os.path.dirname(dataset_dict['file_name']), ('.').join(os.path.basename(dataset_dict['file_name']).split('.')[:-1]))
                    if os.path.exists(extra_files_path):
                        completed_meta[hKey]['raw'][completed_meta[hKey]['datasets'][dataset_id]['job_name']].append({'path' : extra_files_path, 'destination' : destination_name, 'sample_name' : history_dict['sampleName'], 'data_type' : data_type, 'link_name' : link_name, 'filename' : data_type.upper()})

                    #job_names[data_type.upper()] = 0
                    job_names[completed_meta[hKey]['datasets'][dataset_id]['job_name']] = 0

                elif data_type == 'html':
                    # get path of metrics file in case of Picard tools
                    extra_files_path = "%s/%s_files" % (os.path.dirname(dataset_dict['file_name']), ('.').join(os.path.basename(dataset_dict['file_name']).split('.')[:-1]))
                    metrics_file = None

                    for txt_file in glob.glob("%s/*.metrics.txt" %  extra_files_path):
                        metrics_file = txt_file
                    completed_meta[hKey]['raw'][dataset_dict['job_name']] = []
                    job_names[completed_meta[hKey]['datasets'][dataset_id]['job_name']] = 0
                    for txt_file in glob.glob("%s/*" % extra_files_path):
                        sub_data_type = os.path.basename(txt_file).split('.')[-1]
                        if sub_data_type is None:
                            sub_data_type = 'txt'
                        file_name = os.path.basename(txt_file)
                        link_name = "%s/%s" % (history_dict['sampleName'], file_name)
                        destination_name = "%s/%s" % (options.output_dir, link_name)
                        completed_meta[hKey]['raw'][dataset_dict['job_name']].append({'path' : txt_file, 'destination' : destination_name, 'sample_name' : history_dict['sampleName'], 'data_type' : sub_data_type, 'link_name' : link_name, 'filename' : file_name})

                    #print "METRICS_FILE PATH: %s" % metrics_file 
                    # get the information from the txt file
                    # assume that it's a Picard metrics file and get the data
                    if metrics_file is not None:
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
                    file_name = "%s.%s.%s" % (history_dict['sampleName'], dataset_dict['job_name'], data_type)
                    link_name = "%s/%s" % (history_dict['sampleName'], file_name)
                    destination_name = "%s/%s" % (options.output_dir, link_name)

                    # move the file to the destination and include the extra_files_path
                    #shutil.move(dataset_dict['file_name'], destination_name)
                    if dataset_dict['job_name'] in history_dict['raw']:
                        completed_meta[hKey]['raw'][dataset_dict['job_name']].append({'path' : dataset_dict['file_name'], 'destination' : destination_name, 'sample_name' : history_dict['sampleName'], 'data_type' : data_type, 'link_name' : link_name, 'filename' : data_type.upper(), 'extra_files_path' : extra_files_path})
                    else:
                        completed_meta[hKey]['raw'][dataset_dict['job_name']] = [{'path' : dataset_dict['file_name'], 'destination' : destination_name, 'sample_name' : history_dict['sampleName'], 'data_type' : data_type, 'link_name' : link_name, 'filename' : data_type.upper(), 'extra_files_path' : extra_files_path}]
                    job_names[dataset_dict['job_name']] = 0

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
            if ":" in job_name:
                pieces = job_name.split(":")
                job_name = "<br>".join(pieces)
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
                    if os.path.isdir(raw_dataset['path']):
                        print os.path.basename(raw_dataset['path'])
                        source_dir_name = os.path.basename(raw_dataset['path'])
                        shutil.copytree(raw_dataset['path'], "%s/%s/%s" % (options.output_dir, completed_meta[hKey]['sampleName'], source_dir_name))
                    else:
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

        # change the file permissions in the output directory to readable and executable
        for dirpath, dirnames, filenames in os.walk(options.output_dir):
            for filename in filenames:
                path = os.path.join(dirpath, filename)
                os.chmod(path, 0o755)

        # Now it is safe to delete the output histories
        for hKey in completed_meta:
            print "deleting history %s: %s" % (hKey, completed_meta[hKey]['sampleName'])
            #gi.histories.delete_history(hKey, purge=True)

except Exception as e:
    print "There is some kind of issue here: %s" % (e)
    sys.exit()



