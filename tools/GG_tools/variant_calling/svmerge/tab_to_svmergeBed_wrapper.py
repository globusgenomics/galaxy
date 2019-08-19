#!/usr/bin/env python
"""
Converts BED, tabular, VCF file to SVmerge BED format
usage: tab_to_svmergeBed_wrapper.py [options]
"""

import glob, fnmatch, time, re, optparse, os, sys, tempfile, shutil
from subprocess import *

def get_bam_name( bamfile ):
    ## read the BAM file header to get the sample ID name
    cmd = "samtools view -H %s" % bamfile
    #print cmd
    proc = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE )
    for line in proc.stdout:
        #line = proc.stdout.readline()
        if "SM:" in line:
            # split by tabs
            for field in line.rsplit("\t"):
                if "SM:" in field:
                    matchObj = re.search( r'SM:(.*)', field)
                    if matchObj:
                        sampleName = matchObj.group(1)
                        return sampleName

def file_line_count(fname):
    with open(fname) as f:
        for i, l in enumerate(f, 1):
            pass
    return i

def __main__():
    print "\n\nStart time:"
    print time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))

    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '-b', '--bam-input', dest='bamFile', help='The input BAM dataset' )
    parser.add_option( '', '--input-calls-file', dest='inputCallsFile', help='The input Bed file with the target SV regions.' )
    parser.add_option( '', '--input-calls-dataset', dest='inputCallsDataset', help='The input Calls file stored as a dataset.' )
    parser.add_option( '', '--svtype-column', dest='svtype_column', help='SV Type column in the input file.' )
    parser.add_option( '', '--svtype', dest='svtype', help='SV Type of the input file.' )
    parser.add_option( '', '--output', dest='output_bed', help='output assembled directory or file?' )
    parser.add_option( '', '--species', dest='species', help='sample species' )
    ( options, args ) = parser.parse_args()

    sample_name = get_bam_name(options.bamFile)

    ## Create the named.merged.tab file from the BED files inputs
    calls_file = None
    if options.inputCallsFile:
        svCaller = "SVcaller"
        calls_file = options.inputCallsFile
        sv_chr_col = 1
        sv_start_col = 2
        sv_end_col = 3
        sv_header = 1
        if options.svtype_column:
            sv_type_col = options.svtype_column
            sv_type_desc = "exact"
        else:
            sv_type_col = options.svtype
            sv_type_desc = "value"

        calls_file = options.output_bed
        merged_call_file = open(calls_file,'w')
        sampleSize = file_line_count('%s' % (options.inputCallsFile))

        callFileFH = open(options.inputCallsFile, 'r')
        header = ""
        if sv_header == 1:
            header = callFileFH.readline()
        for line in iter(callFileFH):
            line = line.rstrip("\n")
            values = line.split('\t')
            sv_type = None
            if sv_type_desc == "exact":
                sv_type = values[sv_type_col-1]
            elif sv_type_desc == "value":
                sv_type = sv_type_col
            #print sv_type
            merged_call_file.write('%s\t%s\t%s\t%s_%s_%s_%s\n' % (values[sv_chr_col-1], values[sv_start_col-1], values[sv_end_col-1], sv_type, svCaller, sample_name, sampleSize ))
        callFileFH.close()    # close the input BED file
        merged_call_file.close()  # close the merged BED file input for SVmerge


    else:
        ## prepare one calls file from the dataset named merged.calls.tab
        input_calls_dir = '%s/%s_files' % (os.path.dirname(options.inputCallsDataset), ('.').join(os.path.basename(options.inputCallsDataset).split('.')[:-1]))

        ## get the information from the html file
        datasetFile = open(options.inputCallsDataset, 'r')
        callFiles = list()
        sampleSize = 0
        svCaller = None
        for line in iter(datasetFile):
            if "SV_Caller:" in line:
                matchObj = re.search( r'SV_Caller: (.*)</p>', line)
                if matchObj:
                    svCaller = matchObj.group(1)
            elif "BED" in line:
                matchObj = re.search( r'\"(.*)\"', line)
                if matchObj:
                    callFiles.append('%s/%s' % (input_calls_dir, matchObj.group(1)))
                    sampleSize += file_line_count('%s/%s' % (input_calls_dir, matchObj.group(1)))
            elif "VCF" in line:
                matchObj = re.search( r'\"(.*)\"', line)
                if matchObj:
                    callFiles.append('%s/%s' % (input_calls_dir, matchObj.group(1)))
                    sampleSize += file_line_count('%s/%s' % (input_calls_dir, matchObj.group(1)))
        datasetFile.close()

        ## foreach type of SVCALLER prepare the tab file
        print svCaller
        if svCaller == "forestSV":
            sv_type_col = 6
            sv_type_desc = "exact"
            sv_chr_col = 1
            sv_start_col = 2
            sv_end_col = 3
            sv_header = 1
        elif svCaller == "RetroSeq":
            sv_type_col = "INS"
            sv_type_desc = "value"
            sv_chr_col = 1
            sv_start_col = 2
            sv_end_col = 3
            sv_header = 1

        calls_file = options.output_bed
        merged_call_file = open(calls_file,'w')

        for callFile in callFiles:
            callFileFH = open(callFile, 'r')
            header = ""
            if sv_header == 1:
                header = callFileFH.readline()
            for line in iter(callFileFH):
                line = line.rstrip("\n")
                values = line.split('\t')
                sv_type = None
                if sv_type_desc == "exact":
                    sv_type = values[sv_type_col-1]
                elif sv_type_desc == "value":
                    sv_type = sv_type_col
                #print sv_type
                merged_call_file.write('%s\t%s\t%s\t%s_%s_%s_%s\n' % (values[sv_chr_col-1], values[sv_start_col-1], values[sv_end_col-1], sv_type, svCaller, sample_name, sampleSize ))
            callFileFH.close()    # close the input BED file
        merged_call_file.close()  # close the merged BED file input for SVmerge

    print "\n\nEnd time:"
    print time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))

if __name__=="__main__": __main__()

