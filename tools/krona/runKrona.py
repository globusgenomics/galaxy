#!/usr/bin/env python
import argparse, os, shutil, subprocess, sys, tempfile

def stop_err( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit()

def format_mothur_file (infile, sample_name):
    cmd = "awk 'BEGIN{FS=\"\t\";OFS=\"\t\"}{print $2,$3}' %s > %s;sed -i 's/;/\t/g' %s" % (infile, sample_name, sample_name)
    os.system(cmd)

def __main__():

    #Parse Command Line
    parser = argparse.ArgumentParser(description='Wrapper for the Krona visualization for multiple inputs.')
    parser.add_argument( '-m', '--mothur', action='store_true', default=False, help='Is file used from MOTHUR.' )
    parser.add_argument( '-o', '--output', dest='output' )
    parser.add_argument( '-d', '--dataset', nargs=2, action='append', help="<input_file_path> <sample_name>")
    args = parser.parse_args()

    sample_list = []
    for (path, name) in args.dataset:
        sample_name = name.replace(" ", "_")
        sample_list.append(sample_name)
        if args.mothur:
            format_mothur_file(path, sample_name)
        else:
            shutil.copy(path, sample_name)

    os.system("ktImportText %s" % " ".join(sample_list))
    
if __name__=="__main__": __main__()
