#!/usr/bin/env python
#arodri7

"""
A wrapper script for parsing swift logs
"""

import sys, optparse, re, os, tempfile, subprocess, shutil
path = "/mnt/galaxyTools/tools/pymodules/python2.7/lib/python"
sys.path.append(path)
from datetime import datetime
import time
import numpy as np
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

    
#Parse Command Line     
parser = optparse.OptionParser()
parser.add_option( '-l', '--log', dest='log_file', type="string", action="append", help='Swift log file' )
parser.add_option( '-o', '--output', dest='output_file', type="string", help='output file' )
(options, args) = parser.parse_args()

jobs = {}
job_order = []
start_times = []
end_times = []
legend_map = []
FMT = "%Y-%m-%d %H:%M:%S"
group_number = 1
for log in options.log_file:
    start = None
    fh = open(log, "r")
    name = os.path.basename(log).split("-")[0]
    legend_map.append([len(job_order), name])
    for line in fh:
        if "DEBUG Loader arguments" in line:
            searchObj = re.search( r'(.*) (.*),.* DEBUG .*', line, re.M|re.I)
            date = searchObj.group(1)
            mytime = searchObj.group(2)
            start = "%s %s" % (date, mytime)
            start_times.append(int(time.mktime(time.strptime(start, FMT))))
        elif "THREAD_ASSOCIATION" in line:
            searchObj = re.search( r'(.*) (.*),.* DEBUG swift THREAD_ASSOCIATION jobid=(.*) thread.*', line, re.M|re.I)
            date = searchObj.group(1)
            mytime = searchObj.group(2)
            job_id = searchObj.group(3)
            jobs[job_id] = {'THREAD_ASSOCIATION' : "%s %s" % (date, mytime), 'group' : group_number }
            job_order.append(job_id)
        elif "submitting urn:" in line:
            searchObj = re.search( r'(.*) (.*),.* INFO  Cpu .*_swiftwrap (.*) -jobdir .*', line, re.M|re.I)
            date = searchObj.group(1)
            mytime = searchObj.group(2)
            job_id = searchObj.group(3)
            jobs[job_id]['job_start'] = "%s %s" % (date, mytime)
            jobs[job_id]['job_start_epoch'] = int(time.mktime(time.strptime(jobs[job_id]['job_start'], FMT)))
        elif "JOB_END" in line:
            line = line.rstrip("\n")
            searchObj = re.search( r'(.*) (.*),.* DEBUG swift JOB_END jobid=(.*)', line, re.M|re.I)
            date = searchObj.group(1)
            mytime = searchObj.group(2)
            job_id = searchObj.group(3)
            jobs[job_id]['job_end'] = "%s %s" % (date, mytime)
            jobs[job_id]['job_end_epoch'] = int(time.mktime(time.strptime(jobs[job_id]['job_end'], FMT)))

    searchObj = re.search( r'(.*) (.*),.* INFO .*', line, re.M|re.I)
    date = searchObj.group(1)
    mytime = searchObj.group(2)
    end = "%s %s" % (date, mytime)
    end_times.append(int(time.mktime(time.strptime(end, FMT))))
    group_number += 1

#start_epoch = int(time.mktime(time.strptime(start, FMT)))
start_epoch = sorted(start_times)[0]
#end_epoch = int(time.mktime(time.strptime(end, FMT)))
end_epoch = sorted(end_times)[-1]

total_cpu_seconds = 0
plot_pts = []

bottom = []
width = []
left = []
colors = []
number = 1
cmap = matplotlib.cm.get_cmap('prism')

for job_id in job_order:
    group_number = jobs[job_id]['group']
    delta = datetime.strptime(jobs[job_id]['job_end'], FMT) - datetime.strptime(jobs[job_id]['job_start'], FMT)
    jobs[job_id]['delta'] = delta.total_seconds()
    total_cpu_seconds += jobs[job_id]['delta']
    bottom.append(number)
    colors.append(cmap(group_number*5))
    number += 1
    width.append(jobs[job_id]['delta']/3600)
    left.append((float(jobs[job_id]['job_start_epoch']) - float(start_epoch))/3600)

#total_delta = datetime.strptime(end, FMT) - datetime.strptime(start, FMT)
#total_execution = total_delta.total_seconds()
total_delta = end_epoch - start_epoch
total_execution = total_delta

seconds = total_cpu_seconds
total_cpu_hours = seconds // 3600
total_cpu_minutes = (seconds % 3600) // 60
total_cpu_seconds = seconds % 60
#print ("TOTAL CPU TIME: %d hours, %d minutes, %d seconds" % (int(total_cpu_hours), int(total_cpu_minutes), int(total_cpu_seconds)))

seconds = total_execution
total_hours = seconds // 3600
total_minutes = (seconds % 3600) // 60
total_seconds = seconds % 60
#print ("TOTAL EXECUTION TIME: %d hours, %d minutes, %d seconds" % (int(total_hours), int(total_minutes), int(total_seconds)))
# plot the results
with PdfPages(options.output_file) as pdf:
    fig = plt.figure()
    ax = fig.add_axes((.1,.4,.8,.5))
    bp = ax.barh(bottom, width, left=left, color=colors, height=0.4)
    ax.set_title('Swift Log: Job Duration')
    ax.set_xlabel('CPU Hours')
    ax.set_ylabel('Jobs')
    ax.legend( (bp[x[0]] for x in legend_map), (x[1] for x in legend_map), loc='lower center', bbox_to_anchor=(0.5, -0.35), fontsize = 'x-small' )
    fig.text(0.1,0.1, "TOTAL CPU TIME: %d hours, %d minutes, %d seconds\nTOTAL EXECUTION TIME: %d hours, %d minutes, %d seconds" % (int(total_cpu_hours), int(total_cpu_minutes), 
int(total_cpu_seconds), int(total_hours), int(total_minutes), int(total_seconds)), style='italic')
    pdf.savefig(fig)
    
    d = pdf.infodict()
    d['Title'] = 'Swift Job Log'
    d['Author'] = 'GlobusGenomics'
    d['Subject'] = 'Bar plot of all jobs launched by Swift'
    d['Keywords'] = 'Swift barh log history'
    d['CreationDate'] = datetime.today()

