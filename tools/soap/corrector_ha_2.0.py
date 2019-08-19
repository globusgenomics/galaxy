"""
A Galaxy wrapper script for corrector
Peter Li - GigaScience and BGI-HK
"""

import optparse
import os
import shutil
import subprocess
import sys
import tempfile
import glob

def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit()

def cleanup_before_exit(tmp_dir):
    if tmp_dir and os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

def main():
    #Parse command line
    parser = optparse.OptionParser()
    #List of params
    parser.add_option("", "--filelist", type="string", dest="filelist")
    parser.add_option("", "--freq_gz", type="string", dest="freq_gz")

    parser.add_option("", "--default_full_settings_type", dest="default_full_settings_type")
    #Custom params
    parser.add_option("-k", "--kmer_size", dest="kmer_size")
    parser.add_option("-l", "--low_freq_cutoff", dest="low_freq_cutoff")
    parser.add_option("-m", "--min_length_high_freq_region", dest="min_length_high_freq_region")
    parser.add_option("-c", "--max_read_change", dest="max_read_change")
    parser.add_option("-n", "--max_node_num", dest="max_node_num")
    parser.add_option("-a", "--remove_suspicious_data", dest="remove_suspicious_data")
    parser.add_option("-Q", "--ascii_shift_quality_value", dest="ascii_shift_quality_value")
    parser.add_option("-e", "--trim_suspicious_end_regions_Q", dest="trim_suspicious_end_regions_Q")
    parser.add_option("-w", "--trim_error_bases_Q", dest="trim_error_bases_Q")
    parser.add_option("-q", "--qual_threshold_error_bases", dest="qual_threshold_error_bases")
    parser.add_option("-x", "--length_trim_low_qual_ends", dest="length_trim_low_qual_ends")
    parser.add_option("-r", "--min_length_trimmed_read", dest="min_length_trimmed_read")
    parser.add_option("-t", "--thread_num", dest="thread_num")
    parser.add_option("-j", "--convert_reads_into_paired_end_file", dest="convert_reads_into_paired_end_file")
    parser.add_option("-o", "--output_format", dest="output_format")

    #Multiple outputs; number not known before job execution
    parser.add_option("", "--output1.id", dest='output1_id')
    parser.add_option("", "--output1", dest='output1')
    parser.add_option("", "--__new_file_path__", dest='__new_file_path__')
    opts, args = parser.parse_args()

    #Temp directory for data processing
    temp_dir = tempfile.mkdtemp()

    #Files for std out and std error
    tmp_out_file = tempfile.NamedTemporaryFile(dir=temp_dir).name
    tmp_stdout = open(tmp_out_file, 'wb')
    tmp_err_file = tempfile.NamedTemporaryFile(dir=temp_dir).name
    tmp_stderr = open(tmp_err_file, 'wb')

    #Set up command line call
    if opts.default_full_settings_type == "default":
        cmd = "Corrector_HA_v2.0 %s %s" % (opts.freq_gz, opts.filelist)
    elif opts.default_full_settings_type == "full":
        cmd = "Corrector_HA_v2.0 %s %s -k %s -l %s -m %s -c %s -n %s -a %s -Q %s -e %s -w %s -q %s -r %s -t %s -j %s -o %s" % (opts.freq_gz, opts.filelist, opts.kmer_size, opts.low_freq_cutoff, opts.min_length_high_freq_region, opts.max_read_change, opts.max_node_num, opts.remove_suspicious_data, opts.ascii_shift_quality_value, opts.trim_suspicious_end_regions_Q, opts.trim_error_bases_Q, opts.qual_threshold_error_bases, opts.min_length_trimmed_read, opts.thread_num, opts.convert_reads_into_paired_end_file, opts.output_format)

        if opts.length_trim_low_qual_ends != "":
            cmd =  cmd + " -x %s" % opts.length_trim_low_qual_ends

    print "Command executed: ", cmd

    buffsize = 1048576

    #Temp directory to perform processing
    dirpath = tempfile.mkdtemp()
    print "Working directory: ", dirpath

    try:
        #Execution occurs in the directory where the input read files are
        proc = subprocess.Popen(args=cmd, shell=True, cwd=dirpath, stdout=tmp_stdout, stderr=tmp_stderr)
        returncode = proc.wait()
        #Get stdout, allowing for case where it's very large
        tmp_stdout = open(tmp_out_file, 'rb')
        stdout = ''
        try:
            while True:
                stdout += tmp_stdout.read(buffsize)
                if not stdout or len(stdout) % buffsize != 0:
                    break
        except OverflowError:
            pass

        print stdout

        #Get stderr, allowing for case where it's very large
        tmp_stderr = open(tmp_err_file, 'rb')
        stderr = ''
        try:
            while True:
                stderr += tmp_stderr.read(buffsize)
                if not stderr or len(stderr) % buffsize != 0:
                    break
        except OverflowError:
            pass

        #Close streams
        tmp_stdout.close()
        tmp_stderr.close()
        if returncode != 0:
            raise Exception, stderr

    except Exception, e:
        raise Exception, 'Problem performing Corrector process: ' + str(e)

    #Read Corrector results into outputs
    print "Unique identifier for file: " + opts.output1_id
    print "Files kept in: " + opts.__new_file_path__

    #Excel file output
    xls_out = open(opts.output1, 'w')
    xlspath = opts.filelist + ".QC.xls"
    xls_in = open(xlspath, 'r')
    data = xls_in.read()
    xls_out.write(data)
    xls_out.close()
    xls_in.close()

    #Create outputs; need to move and rename files for galaxy for display multiple files
    print "Reading filelist contents"
    file = open(opts.filelist)
    index = 1
    #Check file format
    if opts.output_format == "0":
        format = ".fa.gz"
    elif opts.output_format =="1":
        format = ".fq.gz"
    elif opts.output_format =="2":
        format = ".fa"
    elif opts.output_format =="3":
        format = ".fq"

    #Read the file paths in read.lst
    for line in file:
        print "line:", line
        #Work with cor.pair.fq.gz files
        print "Working on cor.pair files"
        #Create path to access file
        source = line.rstrip() + ".cor.pair_" + str(index) + format
        print "Renaming file: ", source
        #Create string for renaming file
        dir = line[0:line.rindex("/")]
        filename = line[line.rindex("/") + 1:].rstrip()
        filename = filename.replace("_", ".")
        dest = dir + "/primary_" + opts.output1_id + "_" + filename + ".cor.pair." + str(index) + "_visible_" + format
        print "New file name: ", dest
        #Rename and move file
        os.rename(source, dest)
        shutil.move(dest, opts.__new_file_path__)

        #Deal with cor.stat files
        print "Working on cor.stat files"
        #Create path to access file
        source = line.rstrip() + ".cor.stat"
        print "Renaming file: ", source
        #Create string for renaming file
        dir = line[0:line.rindex("/")]
        filename = line[line.rindex("/") + 1:].rstrip()
        filename = filename.replace("_", ".")
        dest = dir + "/primary_" + opts.output1_id + "_" + filename + ".cor.stat_visible_txt"
        print "New file name: ", dest
        #Rename and move file
        os.rename(source, dest)
        shutil.move(dest, opts.__new_file_path__)

        #Deal with cor single fq gz files if present
        print "Working on cor single fq gz files"
        #Create path to access file
        source = line.rstrip() + ".cor.single.fq.gz"
        print "Renaming file: ", source
        #Need to check that this file is present
        if os.path.isfile(source):
            #Create string for renaming file
            dir = line[0:line.rindex("/")]
            filename = line[line.rindex("/") + 1:].rstrip()
            filename = filename.replace("_", ".")
            dest = dir + "/primary_" + opts.output1_id + "_" + filename + ".cor.single_visible_" + format
            print "New file name: ", dest
            #Rename and move file
            os.rename(source, dest)
            shutil.move(dest, opts.__new_file_path__)

        #Deal with cor pair single stat files if present
        print "Working on cor single single stat files"
        #Create path to access file
        source = line.rstrip() + ".cor.pair.single.stat"
        print "Renaming file: ", source
        #Need to check that this file is present
        if os.path.isfile(source):
            #Create string for renaming file
            dir = line[0:line.rindex("/")]
            filename = line[line.rindex("/") + 1:].rstrip()
            filename = filename.replace("_", ".")
            dest = dir + "/primary_" + opts.output1_id + "_" + filename + ".cor.pair.single.stat_visible_txt"
            print "New file name: ", dest
            #Rename and move file
            os.rename(source, dest)
            shutil.move(dest, opts.__new_file_path__)

        index = index + 1
    file.close()

    #Clean up temp files
    cleanup_before_exit(temp_dir)
    cleanup_before_exit(dirpath)
    #Check results in output file
    if os.path.getsize(opts.output1) > 0:
        sys.stdout.write('Status complete')
    else:
        stop_err("The output is empty")

if __name__ == "__main__": main()
