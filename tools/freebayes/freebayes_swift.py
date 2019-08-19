#!/usr/bin/env python

"""
A wrapper script for running freebayes (take directory as input).
"""

import glob, sys, optparse, os, tempfile, subprocess, shutil
from binascii import unhexlify
from string import Template
import time  

CHUNK_SIZE = 2**20 #1mb

def get_base_count(bed_file):
    base_count = 0
    with open(bed_file, "r") as inF:
        for line in inF:
            line = line.rstrip("\n")
            values = line.split("\t")
            base_count += int(values[2]) - int(values[1]) + 1
    return base_count

def split_bed_per_base(bed_file, parallel_split, bed_tmp):
    num_bases = get_base_count(bed_file)
    bases_per_file = int(float(num_bases)/float(parallel_split))+ 1
    base_count = 0
    file_count = 0
    file_list = []
    with open(bed_file, "r") as inF:
        try:
            w = open("%s/tmp_%s.bed" % (bed_tmp, file_count), 'w')
            for line in inF:
                line = line.rstrip("\n")
                values = line.split("\t")
                bases_in_line = int(values[2]) - int(values[1]) + 1
                if base_count + bases_in_line <= bases_per_file:
                    w.write(line + "\n")
                    base_count += bases_in_line
                else:
                    # get the difference between max bases per file and current base count
                    space_left = bases_per_file - base_count
                    segment_end = int(values[1]) + int(space_left) - 1
                    next_segment_start = int(segment_end) + 1
                    new_line = "%s\t%s\t%s\n" % (values[0], values[1], segment_end)
                    w.write(new_line)
                    w.close()

                    # start new bed file
                    base_count = 0
                    file_list.append("%s/tmp_%s.bed" % (bed_tmp, file_count) )
                    file_count += 1
                    line_count = 0
                    w = open("%s/tmp_%s.bed" % (bed_tmp, file_count), 'w')
                    first_line = "%s\t%s\t%s\n" % (values[0], next_segment_start, values[2])
                    w.write(first_line)
        finally:
            file_list.append("%s/tmp_%s.bed" % (bed_tmp, file_count) )
            w.close()
    return file_list

def split_bed(bed_file, parallel_split, bed_tmp):
    num_lines = sum(1 for line in open(bed_file))
    lines_per_file = int(float(num_lines)/float(parallel_split))+ 1
    line_count = 0
    file_count = 0
    file_list = []
    with open(bed_file, "r") as inF:
        try:
            w = open("%s/tmp_%s.bed" % (bed_tmp, file_count), 'w')
            for line in inF:
                if line_count > lines_per_file:
                    file_list.append("%s/tmp_%s.bed" % (bed_tmp, file_count) )
                    w.close()
                    file_count += 1
                    line_count = 0
                    w = open("%s/tmp_%s.bed" % (bed_tmp, file_count), 'w')
                w.write(line)
                line_count += 1
        finally:
            file_list.append("%s/tmp_%s.bed" % (bed_tmp, file_count) )
            w.close()
    return file_list



def cleanup_before_exit( tmp_dir ):
    if tmp_dir and os.path.exists( tmp_dir ):
        shutil.rmtree( tmp_dir )

def open_file_from_option( filename, mode = 'rb' ):
    if filename:
        return open( filename, mode = mode )
    return None

def index_bam_files( bam_filenames, tmp_dir ):
    for bam_filename in bam_filenames:
        bam_index_filename = "%s.bai" % bam_filename
        if not os.path.exists( bam_index_filename ):
            #need to index this bam file
            stderr_name = tempfile.NamedTemporaryFile( prefix = "bam_index_stderr" ).name
            command = 'samtools index %s %s' % ( bam_filename, bam_index_filename )
            proc = subprocess.Popen( args=command, shell=True, stderr=open( stderr_name, 'wb' ) )
            return_code = proc.wait()
            if return_code:
                for line in open( stderr_name ):
                    print >> sys.stderr, line
                os.unlink( stderr_name ) #clean up
                cleanup_before_exit( tmp_dir )
                raise Exception( "Error indexing BAM file" )
            os.unlink( stderr_name ) #clean up

def __main__():
    #liubo added, print current time for calculating execution time of GATK
    print "Start time:"
    print time.strftime('%d/%m/%Y %H:%M:%S\n', time.localtime(time.time()))

    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '-p', '--pass_through', dest='pass_through_options', action='append', type="string", help='These options are passed through directly, without any modification.' )
    parser.add_option( '', '--input_dir', dest='input_dir', action='store', type="string", help='Input directory path of BAM files' )
    parser.add_option( '', '--input_dir_file', dest='input_dir_file', action='store', type="string", help='Input directory File containing path of BAM files' )
    parser.add_option( '', '--targets', dest='targets_file', action='store', type="string", help='BED file of regions to look at' )
    parser.add_option( '', '--parallel-split', dest='parallel_split', action='store', type="string", help='number of processes to split to' )
    parser.add_option( '', '--out-dir', dest='output_dir', action='store', type="string", default=None, help='output directory.' )
    parser.add_option( '', '--out-vcf', dest='out_vcf', action='store', type="string", default=None, help='vcf output.' )
    parser.add_option( '', '--out-log', dest='swift_log', action='store', type="string", default=None, help='log output.' )
    parser.add_option( '', '--trace', dest='trace_output', action='store', type="string", default=None, help='trace output.' )
    parser.add_option( '', '--failed-alleles', dest='failed_alleles_output', action='store', type="string", default=None, help='failed alleles output.' )
    parser.add_option( '', '--fasta-reference', dest='fasta_reference', action='store', type="string", default=None, help='Fasta reference should have a fai file, otherwise create one' )   
    (options, args) = parser.parse_args()

    swift_bin = 'swift'
    sites_file = '/opt/galaxy/tools/freebayes/freebayes_sites.xml'
    tc_file = '/opt/galaxy/tools/swift/tc.data'
    swift_file = '/opt/galaxy/tools/freebayes/freebayes.swift'
    samtools_bin = "/mnt/galaxyTools/tools/samtools/1.2/bin/samtools"
    freebayes_bin = "/mnt/galaxyTools/tools/freebayes/v0.9.15-1-g076a2a2"

    #set output directories workspace
    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)
    bed_tmp = "%s/beds" % options.output_dir
    if not os.path.exists(bed_tmp):
        os.mkdir(bed_tmp)
    inputs_tmp = "%s/inputs" % options.output_dir
    if not os.path.exists(inputs_tmp):
        os.mkdir(inputs_tmp)
    output_dir = "%s/output" % options.output_dir
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    tmp_dir = tempfile.mkdtemp( dir=options.output_dir , prefix='tmp-FREEBAYES-' )

    #input directory
    if options.input_dir_file:
        infile = open(options.input_dir_file, 'r')
        inputDirLine = infile.readline()
        inputDirectory = inputDirLine.rstrip('\n\r ')
    else:
        inputDirectory = options.input_dir

    # set up the reference file
    reference_path = None
    if os.path.exists("%s.fai" % options.fasta_reference):
        reference_path = options.fasta_reference
    else:
        # create a fasta index file
        reference_path = "%s/local_fasta.fa" % os.getcwd()
        os.symlink(options.fasta_reference, reference_path)
        os.system("samtools faidx %s" % reference_path)

    # setup the parallelization
    targets_flag = ""
    if options.targets_file != "chromosome":
        # split the input bed file
        parallel_list = split_bed_per_base(options.targets_file, options.parallel_split, bed_tmp)
        targets_flag = "--targets BED"

    #set up stdout and stderr output options
    stdout = tempfile.NamedTemporaryFile( prefix="freebayes-stdout-", dir=tmp_dir )
    stderr = tempfile.NamedTemporaryFile( prefix="freebayes-stderr-", dir=tmp_dir )

    #traverse all the bam files in the input directory and add them to a list file
    #options.input_dir is the input directory's full path
    inputs_list_file = "%s/inputs.txt" % inputs_tmp
    bam_list_fh = open(inputs_list_file, "w")
    if options.input_dir_file:
	infile = open(options.input_dir_file, 'r')
        inputDirLine = infile.readline()
        inputDirectory = inputDirLine.rstrip('\n\r ')
    else:
        inputDirectory = options.input_dir

    if inputDirectory:
        BAMlist = sorted(glob.glob("%s/*.bam" % inputDirectory))
        bam_list_fh.write("\n".join(BAMlist))
    bam_list_fh.close()

    if options.pass_through_options:
        pass_through_options = ' '.join( options.pass_through_options )

    tool_cmd = "source /mnt/galaxyTools/tools/pymodules/python2.7/env.sh; %s/freebayes -f %s -L %s %s -v OUTPUT "  % (freebayes_bin, reference_path, inputs_list_file, pass_through_options)
    if len(targets_flag) > 0:
        tool_cmd += targets_flag
    if options.trace_output:
        tool_cmd += " --trace TRACE " 
    if options.failed_alleles_output:
        tool_cmd += " --failed-alleles FAILED"

    #####

    # prepare command line
    swift_params = list()
    swift_params.append('-outputdir=' + output_dir)
    swift_params.append('-parallel_list=\"%s\"' % ",".join(parallel_list))
    if options.targets_file != "chromosome":
        swift_params.append('-parallel=bed')
    else:
        swift_params.append('-parallel=chromosome')

    ## construct the swift command
    swift_cmd = "%s -sites.file %s -tc.file %s %s " %   (swift_bin, sites_file, tc_file, swift_file)
    cmd = "%s %s %s" % (swift_cmd, ' '.join(swift_params), '-tool_cmd=\"'+tool_cmd+'\"')
    print cmd
    return_code = None

    if return_code is None or not return_code:
        proc = subprocess.Popen( args=cmd, stdout=stdout, stderr=stderr, shell=True, cwd=tmp_dir )
        return_code = proc.wait()
    if return_code:
        stderr_target = sys.stderr
    else:
        if stdout:
            stderr_target = stdout
        else:
            stderr_target = sys.stdout
    stderr.flush()
    stderr.seek(0)
    while True:
        chunk = stderr.read( CHUNK_SIZE )
        if chunk:
            stderr_target.write( chunk )
        else:
            break
    stderr.close()

    # concatenate all VCF files into one output
    vcf_outputs = glob.glob("%s/*.vcf" % output_dir)
    vcf_string = ""
    for i in range(len(vcf_outputs)):
        vcf_string += "%s/tmp_%s.vcf " % (output_dir,i)
    vcf_cmd = "vcf-concat %s > %s " % (vcf_string, options.out_vcf)
    #set up stdout and stderr output options
    #if no stderr file is specified, we'll use our own
    stderr = tempfile.NamedTemporaryFile( prefix="VCFtools-stderr-", dir=tmp_dir )
    stdout = tempfile.NamedTemporaryFile( prefix="VCFtools-stderr-", dir=tmp_dir )

    # run the cmd
    proc = subprocess.Popen( args=vcf_cmd, stdout=stdout, stderr=stderr, shell=True, cwd=tmp_dir )
    return_code = proc.wait()
    if return_code:
        stderr_target = sys.stderr
    else:
        if stdout:
            stderr_target = stdout
        else:
            stderr_target = sys.stdout
    stderr.flush()
    stderr.seek(0)
    while True:
        chunk = stderr.read( CHUNK_SIZE )
        if chunk:
            stderr_target.write( chunk )
        else:
            break
    stderr.close()

    swift_log_files = glob.glob("%s/*.log" % tmp_dir)
    cmdSummary = "export PYTHONPATH=/mnt/galaxyTools/tools/pymodules/python2.7/lib/python:/opt/galaxy/lib:\$PYTHONPATH;. /mnt/galaxyTools/tools/pymodules/python2.7/env.sh;/opt/galaxy/tools/swift/parse_swift_log.py "
    for logF in swift_log_files:
        if "swift.log" in logF:
            continue
        cmdSummary += " -l %s " % logF
    cmdSummary += " -o %s" % options.swift_log

    return_code = None
    stderr = tempfile.NamedTemporaryFile( prefix="TOOL-stderr-", dir=tmp_dir )
    stdout = tempfile.NamedTemporaryFile( prefix="TOOL-stdout-", dir=tmp_dir )
    if return_code is None or not return_code:
        proc = subprocess.Popen( args=cmdSummary, stdout=stdout, stderr=stderr, shell=True, cwd=tmp_dir )
        return_code = proc.wait()
        if return_code:
            stderr_target = sys.stderr
        else:
            if stdout:
                stderr_target = stdout
            else:
                stderr_target = sys.stdout
        stderr.flush()
        stderr.seek(0)
        while True:
            chunk = stderr.read( CHUNK_SIZE )
            if chunk:
                stderr_target.write( chunk )
            else:
                break
        stderr.close()

    #cleanup_before_exit( options.output_dir )


    #print current time for calculating execution time of GATK
    print "\n\nEnd time:"
    print time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))

if __name__=="__main__": __main__()
