#!/usr/bin/env python


"""
A wrapper script for running SAMTools commands.
"""

import sys, optparse, os, glob, tempfile, subprocess, shutil
from string import Template

GALAXY_EXT_TO_SAMTOOLS_EXT = { 'bam_index':'bam.bai', } #items not listed here will use the galaxy extension as-is
GALAXY_EXT_TO_SAMTOOLS_FILE_TYPE = GALAXY_EXT_TO_SAMTOOLS_EXT #for now, these are the same, but could be different if needed
DEFAULT_SAMTOOLS_PREFIX = "SAMTools_file"
CHUNK_SIZE = 2**20 #1mb


def cleanup_before_exit( tmp_dir ):
    if tmp_dir and os.path.exists( tmp_dir ):
        shutil.rmtree( tmp_dir )

def SAMTOOLS_filename_from_galaxy( galaxy_filename, galaxy_ext, target_dir = None, prefix = None ):
    suffix = GALAXY_EXT_TO_SAMTOOLS_EXT.get( galaxy_ext, galaxy_ext )
    if prefix is None:
        prefix = DEFAULT_SAMTOOLS_PREFIX
    if target_dir is None:
        target_dir = os.getcwd()
    SAMTools_filename = os.path.join( target_dir, "%s.%s" % ( prefix, suffix ) )
    os.symlink( galaxy_filename, SAMTools_filename )
    return SAMTools_filename

def SAMTOOLS_filetype_argument_substitution( argument, galaxy_ext ):
    return argument % dict( file_type = GALAXY_EXT_TO_SAMTOOLS_FILE_TYPE.get( galaxy_ext, galaxy_ext ) )

def open_file_from_option( filename, mode = 'rb' ):
    if filename:
        return open( filename, mode = mode )
    return None

def html_report_from_directory( html_out, dir ):
    html_out.write( '<html>\n<head>\n<title>Galaxy - SAMTOOLS Output</title>\n</head>\n<body>\n<p/>\n<ul>\n' )
    for fname in sorted( os.listdir( dir ) ):
        html_out.write(  '<li><a href="%s">%s</a></li>\n' % ( fname, fname ) )
    html_out.write( '</ul>\n</body>\n</html>\n' )

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

def __main__():
    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '-p', '--pass-through', dest='pass_through_options', action='append', type="string", help='These options are passed through directly to SAMTOOLS, without any modification.' )
    parser.add_option( '-d', '--dataset', dest='datasets', action='append', type="string", nargs=4, help='"-argument" "original_filename" "galaxy_filetype" "name_prefix"' )
    parser.add_option( '', '--stdout', dest='stdout', action='store', type="string", default=None, help='If specified, the output of stdout will be written to this file.' )
    parser.add_option( '', '--out-mpileup', dest='out_mpileup', action='store', type="string", default=None, help='mpileup output.' )
    parser.add_option( '', '--stderr', dest='stderr', action='store', type="string", default=None, help='If specified, the output of stderr will be written to this file.' )
    parser.add_option( '', '--html_report_from_directory', dest='html_report_from_directory', action='append', type="string", nargs=2, help='"Target HTML File" "Directory"')

    ##alex added (3/17/2014)
    parser.add_option( '', '--input_dir', dest='input_dir', action='store', type="string", help='Input directory path of BAM files' )
    parser.add_option( '', '--input_dir_file', dest='input_dir_file', action='store', type="string", help='Input directory File containing path of BAM files' )
    parser.add_option( '', '--parallel-bed', dest='bed_file', action='store', type="string", help='BED file of regions to look at' )
    parser.add_option( '', '--parallel-split', dest='parallel_split', action='store', type="string", help='number of processes to split to' )
    parser.add_option( '', '--out-dir', dest='output_dir', action='store', type="string", default=None, help='If specified, the output directory for extra files.' )
    (options, args) = parser.parse_args()
    
    swift_bin = 'swift'
    sites_file = '/opt/galaxy/tools/samtools/samtools_mpileup_sites.xml'
    tc_file = '/opt/galaxy/tools/swift/tc.data'
    swift_file = '/opt/galaxy/tools/samtools/samtools_mpileup.swift'
    samtools_bin = "/mnt/galaxyTools/tools/samtools/1.1/bin/samtools"
    bcftools_bin = "/mnt/galaxyTools/tools/samtools/1.1/bin/bcftools"

    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)
    bed_tmp = "%s/beds" % options.output_dir
    if not os.path.exists(bed_tmp):
        os.mkdir(bed_tmp)
    output_dir = "%s/output" % options.output_dir
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    tmp_dir = tempfile.mkdtemp( dir=options.output_dir , prefix='tmp-SAMTOOLS-' )

    if options.input_dir_file:
        infile = open(options.input_dir_file, 'r')
        inputDirLine = infile.readline()
        inputDirectory = inputDirLine.rstrip('\n\r ')
    else:
        inputDirectory = options.input_dir

    tool_cmd = "%s  mpileup --output OUTPUT" % samtools_bin
    if options.pass_through_options:
        tool_cmd += ' '.join( options.pass_through_options )

    # setup the parallelization
    if options.parallel_split != "chromosome":
        # split the input bed file
        parallel_list = split_bed_per_base(options.bed_file, options.parallel_split, bed_tmp)
        tool_cmd += " --positions BED_FILE"
    else:
        # split the job per chromosome
        # determine the number of chromosomes by looking at the BAM header
        tool_cmd += "-r CHR"
        sample_bam_file = glob.glob("%s/*.bam" % inputDirectory)[0]
        process = subprocess.Popen("samtools view -H %s | grep @SQ", stdout=subprocess.PIPE)
        parallel_list = []
        for line in iter(process.stdout.readline, ""):
            parallel_list.append(line.rstrip("\n").split("\t")[1].split(":")[1])

    #set up stdout and stderr output options
    stdout = open_file_from_option( options.stdout, mode = 'w' )
    stderr = open_file_from_option( options.stderr, mode = 'w' )
    #if no stderr file is specified, we'll use our own
    if stderr is None:
        stderr = tempfile.NamedTemporaryFile( prefix="SAMTOOLS-stderr-", dir=tmp_dir )
 
    if options.datasets:
        for ( dataset_arg, filename, galaxy_ext, prefix ) in options.datasets:
            SAMTools_filename = SAMTOOLS_filename_from_galaxy( filename, galaxy_ext, target_dir = tmp_dir, prefix = prefix )
            if dataset_arg:
                tool_cmd = '%s %s %s ' % ( tool_cmd, SAMTOOLS_filetype_argument_substitution( dataset_arg, galaxy_ext ), SAMTools_filename )

            #auto index fasta files:
            if galaxy_ext == 'fa':
                index_cmd = 'samtools faidx %s' % ( SAMTools_filename )
                proc = subprocess.Popen( args=index_cmd, stdout=stdout, stderr=stderr, shell=True, cwd=tmp_dir )
                return_code = proc.wait()
                if return_code:
                    break

    #####  Alex added this section to include directory as an input for mpileup
    if inputDirectory:
        BAMlist = sorted(glob.glob("%s/*.bam" % inputDirectory))
        tool_cmd += " ".join(BAMlist)

    # create a BCF index file for the output
    tool_cmd += "; %s index BCF_INPUT" % (bcftools_bin)
    #####

    # prepare command line
    swift_params = list()
    swift_params.append('-outputdir=' + output_dir)
    swift_params.append('-parallel_list=\"%s\"' % ",".join(parallel_list))
    if options.parallel_split != "chromosome":
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
    
    # concatenate all BCF files into one output
    bcf_outputs = glob.glob("%s/*.bcf" % output_dir)
    bcf_cmd = "bcftools concat -a -O u -o %s " % (options.out_mpileup)
    for i in range(len(bcf_outputs)):
        bcf_cmd += "%s/tmp_%s.bcf " % (output_dir,i)

    #set up stdout and stderr output options
    stdout = open_file_from_option( options.stdout, mode = 'w' )
    stderr = open_file_from_option( options.stderr, mode = 'w' )
    #if no stderr file is specified, we'll use our own
    if stderr is None:
        stderr = tempfile.NamedTemporaryFile( prefix="SAMTOOLS-stderr-", dir=tmp_dir )

    # run the cmd
    proc = subprocess.Popen( args=bcf_cmd, stdout=stdout, stderr=stderr, shell=True, cwd=tmp_dir )
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
    cmdSummary = "/opt/galaxy/tools/swift/parse_swift_log.py "
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

    #generate html reports
    #if options.html_report_from_directory:
    #    for ( html_filename, html_dir ) in options.html_report_from_directory:
    #        html_report_from_directory( open( html_filename, 'wb' ), html_dir )
    
    #cleanup_before_exit( tmp_dir )

if __name__=="__main__": __main__()
