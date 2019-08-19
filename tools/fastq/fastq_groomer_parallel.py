# Kenneth Sabir
# Garvan Institute
import sys
import time
import os
import math
import subprocess
import string
import shutil
import pickle
import io
from multiprocessing import Process
from galaxy_utils.sequence.fastq import fastqReader, fastqVerboseErrorReader, fastqAggregator, fastqWriter

  

def main():
    split_program = "split"
    cat_program = "cat"
    input_filename =  sys.argv[1]
    output_filename = sys.argv[3]
    number_of_processes = 1;
    if (len(sys.argv) > 7):
        number_of_processes = int(sys.argv[7])
    file_prefix = "temp_groomer_part_"
    
    t1 = time.time()
    old_path = os.getcwd()

    lines_per_process,number_of_lines = calculate_lines_per_process(input_filename, number_of_processes)                 
    temp_dir_name = move_to_temp_dir()
    sequences = number_of_lines/4;
    args = [split_program, "-l"+str(lines_per_process), input_filename, file_prefix]
#    print "The args are: " , args
    subprocess.call(args)
#    print "Finished"
    file_count = 0;
    keep_checking = True
    processes = []
    output_filenames = []
    while keep_checking: 
        
        # only need to support 26x26 different processes, so do it brute force (ie not in a loop) for 2 chars.
        lastchar = string.letters[file_count % len(string.letters)] 
        firstchar = string.letters[(file_count / len(string.letters)) % len(string.letters)] 
        temp_input_filename = "%s%c%c" % (file_prefix, firstchar, lastchar)
        
 #       print 'looking for ' + temp_input_filename
        if os.path.exists(temp_input_filename):
#            print 'found ' + temp_input_filename
            temp_output_filename = temp_input_filename + "_output"
            output_filenames.append(temp_output_filename)
            p = Process(target=partition, args=([temp_input_filename, temp_output_filename, file_count]))
            p.start()
            processes.append(p)
            file_count = file_count + 1
        else:
            break
    for p in processes :
        p.join()
    cat_params = [cat_program]
    cat_params.extend(output_filenames)
    with open(output_filename, 'w') as catOutputFile:
        subprocess.call(cat_params, stdout=catOutputFile)
    summarize_input = sys.argv[6] == 'summarize_input'
    input_type = sys.argv[2]
    output_type = sys.argv[4]
    print "Groomed %i %s reads into %s reads." % ( sequences, input_type, output_type )

    aggregators = []
    if summarize_input:
        for temp_output_filename in output_filenames :
            with open(temp_output_filename + "_summary", 'r') as summaryLogFile:
                temp_aggregator = pickle.load(summaryLogFile)
                aggregators.append(temp_aggregator)

        print_aggregators(aggregators)
    os.chdir(old_path)
    shutil.rmtree(temp_dir_name)
    time2 = time.time()
    print 'Groomer took: %0.3f ms using %d processes' % (((time2 - t1)*1000.0), number_of_processes)

def calculate_lines_per_process(input_filename, number_of_processes):
    wc_program = "wc"
    p = subprocess.Popen([wc_program, "-l", input_filename], stdout=subprocess.PIPE)
    out, err = p.communicate()
    number_of_lines = int(string.split(string.lstrip(out), ' ', 1)[0])
    exact_lines_per_process = number_of_lines * 1.0 / number_of_processes
    lines_per_process = int(math.ceil((exact_lines_per_process / 4.0))) * 4
    return lines_per_process,number_of_lines

def move_to_temp_dir():
    dirExists = False;
    dir_name = None
    
    while not dirExists:
        dir_name = "temp_groomer_part_" + str(time.time())
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)
            break;
    os.chdir(dir_name)
    return dir_name

def print_aggregators(aggregators):
    total_ascii_range = [None, None]
    total_decimal_range = [None, None]
    total_valid_formats = set()
    for aggregator in aggregators:
#        print "This aggregators valid formats are: " + str(aggregator.get_valid_formats())
        total_valid_formats = total_valid_formats.union(set(aggregator.get_valid_formats()))
        ascii_range = aggregator.get_ascii_range()
        decimal_range =  aggregator.get_decimal_range()
        
        if total_ascii_range[0] is None:
            total_ascii_range[0] = ascii_range[0]
        else:
            total_ascii_range[0] = min (total_ascii_range[0], ascii_range[0]) 

        # max of None and a value is the value
        total_ascii_range[1] = max (total_ascii_range[1], ascii_range[1]) 
        if total_decimal_range[0] is None:
            total_decimal_range[0] = decimal_range[0]
        else:
            total_decimal_range[0] = min (total_decimal_range[0], decimal_range[0]) 
        # max of None and a value is the value
        total_decimal_range[1] = max (total_decimal_range[1], decimal_range[1]) 
    print "total_valid_formats= " + str(total_valid_formats)
    print "Based upon quality and sequence, the input data is valid for: %s" % ( ", ".join( total_valid_formats )  or "None" )
    print "Input ASCII range: %s(%i) - %s(%i)" % ( repr( total_ascii_range[0] ), ord( total_ascii_range[0] ), repr( total_ascii_range[1] ), ord( total_ascii_range[1] ) ) #print using repr, since \x00 (null) causes info truncation in galaxy when printed
    print "Input decimal range: %i - %i" % ( total_decimal_range[0], total_decimal_range[1] )
                

def partition(input_filename, temp_output_filename, fileCount):
#    print 'Starting Thread: ' + str(fileCount)
#    input_filename = sys.argv[1]
    input_type = sys.argv[2]
    output_type = sys.argv[4]
    force_quality_encoding = sys.argv[5]
    summarize_input = sys.argv[6] == 'summarize_input'
    if force_quality_encoding == 'None':
        force_quality_encoding = None
    aggregator = fastqAggregator()
    temp_process_file = fastqWriter( open( temp_output_filename, 'wb'), format = output_type, force_quality_encoding = force_quality_encoding )
    read_count = None
    if summarize_input:
        reader = fastqVerboseErrorReader
    else:
        reader = fastqReader
    for read_count, fastq_read in enumerate( reader( open(input_filename, 'rb'), format = input_type, apply_galaxy_conventions = True ) ):
        if summarize_input:
            aggregator.consume_read( fastq_read )
        temp_process_file.write( fastq_read )
#        print "Just wrote (%d): " % read_count + str(fastq_read)
    temp_process_file.close()
    if read_count is not None:
        if input_type != output_type and 'solexa' in [ input_type, output_type ]:
            print "Converted between Solexa and PHRED scores."
        if summarize_input:
            with open(temp_output_filename + "_summary", 'w') as summaryLogFile :
                pickle.dump(aggregator, summaryLogFile)
    else:
        print "No valid FASTQ reads were provided."

if __name__ == "__main__": main()

