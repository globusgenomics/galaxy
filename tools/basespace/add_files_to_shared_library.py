#!/usr/bin/env python
import time
from datetime import datetime, timedelta
import optparse, os, shutil, subprocess, sys, tempfile, re
from bioblend import galaxy
import requests
requests.packages.urllib3.disable_warnings()

def __main__():
    print "Start time:"
    print time.strftime('%d/%m/%Y %H:%M:%S\n', time.localtime(time.time()))

    descr = "add_files_to_shared_library.py: Add files in a directory to a Shared library"
    parser = optparse.OptionParser(description=descr)
    parser.add_option( '-p', '--input-dir', dest="input_dir", help="Path to directory data" )
    parser.add_option( '-s', '--sampleID', dest="libName", help="Library name will be your Sample ID" )
    parser.add_option( '-o', '--output', dest="output", help="output file" )
    parser.add_option( '-k', '--key', dest='user_api_key', help="The user api key" )
    parser.add_option( '-u', '--url', dest='galaxyurl', help='Base URL' )
    (options, args) = parser.parse_args()

    url = options.galaxyurl
    if "http:" in url:
        url = url.replace("http", "https")

    if os.path.exists(options.input_dir):
        filename = options.output
        directory = options.input_dir
        files_in_lib = os.listdir(directory)
        lib_name = "%s-%s" % (os.path.commonprefix(files_in_lib), options.libName)
        f = open(filename,'w')
        f.write("Adding files in directory %s to Shared Library %s" % (options.input_dir, options.libName))
        gi =  galaxy.GalaxyInstance(url=url, key=options.user_api_key)
        library = gi.libraries.create_library(lib_name, description="Datasets for basespace Sample transfer: %s" % options.libName)

        for dirName, subdirList, fileList in os.walk(directory):
            f.write('\nFound directory: %s\n' % dirName)
            for fname in fileList:
                f.write('\t%s\n' % fname)
                filepath = "%s/%s" % (dirName, fname)
                fileformat = None
                if fname.endswith("fastq.gz"):
                    fileformat = "fastqsanger"
                else:
                    fileformat = fname.split("\.")[-1]
                gi.libraries.upload_from_galaxy_filesystem(library["id"],filepath,file_type=fileformat,dbkey='hg19',link_data_only="True")
            # Remove the first entry in the list of sub-directories
            # if there are any sub-directories present
            if len(subdirList) > 0:
                del subdirList[0]
        f.close()

    else:
        f = open(filename,'w')
        f.write("Seems like nothing was transferred\n")
        f.close()

if __name__ == "__main__":
    __main__()
