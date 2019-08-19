"""
soapdenovo2_scaff.py
A wrapper script for SOAPdenovo2 scaff module
Copyright   Peter Li - GigaScience and BGI-HK
"""

import optparse, os, shutil, subprocess, sys, tempfile

def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit()

def cleanup_before_exit(tmp_dir):
    if tmp_dir and os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

def main():
    #Parse command line
    parser = optparse.OptionParser()
    #Inputs
    parser.add_option('', '--arc', dest='arc')
    parser.add_option('', '--pegrads', dest='pegrads')
    parser.add_option('', '--pregraph_basic', dest='pregraph_basic')
    parser.add_option('', '--updated_edge', dest='updated_edge')
    parser.add_option('', '--contig', dest='contig')
    parser.add_option('', '--read_in_gap', dest='read_in_gap')
    parser.add_option('', '--read_on_contig', dest='read_on_contig')

    parser.add_option("", "--analysis_settings_type", dest="analysis_settings_type")
    parser.add_option("", "--default_full_settings_type", dest="default_full_settings_type")
    parser.add_option("", "--fill_gaps", dest="fill_gaps")
    parser.add_option("", "--compatible_mode", dest="compatible_mode")
    parser.add_option("", "--unmask_contigs", dest="unmask_contigs")
    parser.add_option("", "--keep_contigs_connected", dest="keep_contigs_connected")
    parser.add_option("", "--ass_visual", dest="ass_visual")
    parser.add_option("", "--gap_len_diff", dest="gap_len_diff")
    parser.add_option("", "--min_contig_len", dest="min_contig_len")
    parser.add_option("", "--min_contig_cvg", dest="min_contig_cvg")
    parser.add_option("", "--max_contig_cvg", dest="max_contig_cvg")
    parser.add_option("", "--insert_size_upper_bound", dest="insert_size_upper_bound")
    parser.add_option("", "--bubble_coverage", dest="bubble_coverage")
    parser.add_option("", "--genome_size", dest="genome_size")
    parser.add_option("", "--ncpu", dest="ncpu")

    #Outputs
    parser.add_option("", "--new_contig_index", dest='new_contig_index')
    parser.add_option("", "--links", dest='links')
    parser.add_option("", "--scaf_gap", dest='scaf_gap')
    parser.add_option("", "--gap_seq", dest='gap_seq')
    parser.add_option("", "--scaf", dest='scaf')
    parser.add_option("", "--scaf_seq", dest='scaf_seq')
    parser.add_option("", "--contig_positions_scaff", dest='contig_positions_scaff')
    parser.add_option("", "--bubble_in_scaff", dest='bubble_in_scaff')
    parser.add_option("", "--scaf_stats", dest='scaf_stats')

    opts, args = parser.parse_args()

    #Need to write inputs to a temporary directory
    dirpath = tempfile.mkdtemp()
    arc_data = open(opts.arc, 'r')
    arc_file = open(dirpath + "/out.Arc", "w")
    for line in arc_data:
        arc_file.write(line)
    arc_data.close()
    arc_file.close()

    pegrads_data = open(opts.pegrads, 'r')
    pegrads_file = open(dirpath + "/out.peGrads", "w")
    for line in pegrads_data:
        pegrads_file.write(line)
    pegrads_data.close()
    pegrads_file.close()

    pregraph_basic_data = open(opts.pregraph_basic, 'r')
    pregraph_basic_file = open(dirpath + "/out.preGraphBasic", "w")
    for line in pregraph_basic_data:
        pregraph_basic_file.write(line)
    pregraph_basic_data.close()
    pregraph_basic_file.close()

    updated_edge_data = open(opts.updated_edge, 'r')
    updated_edge_file = open(dirpath + "/out.updated.edge", "w")
    for line in updated_edge_data:
        updated_edge_file.write(line)
    updated_edge_data.close()
    updated_edge_file.close()

    contig_data = open(opts.contig, 'r')
    contig_file = open(dirpath + "/out.contig", "w")
    for line in contig_data:
        contig_file.write(line)
    contig_data.close()
    contig_file.close()

    read_in_gap_out = open(dirpath + "/out.readInGap.gz", "wb")
    with open(opts.read_in_gap, mode='rb') as file: # b is important -> binary
        fileContent = file.read()
        read_in_gap_out.write(fileContent)
    read_in_gap_out.close()
    file.close()

    read_on_contig_out = open(dirpath + "/out.readOnContig.gz", "wb")
    with open(opts.read_on_contig, mode='rb') as file: # b is important -> binary
        fileContent = file.read()
        read_on_contig_out.write(fileContent)
    read_on_contig_out.close()
    file.close()

#    f=open(dirpath + "/out.readOnContig.gz", "wb")
#    f.write(opts.read_on_contig)
#    f.close()

#    read_in_gap_out = open(opts.read_in_gap, 'wb')
#    with open(dirpath + "/out.readInGap.gz", mode='rb') as file: # b is important -> binary
#        fileContent = file.read()
#        read_in_gap_out.write(fileContent)
#    read_in_gap_out.close()
#    file.close()
#
#    read_on_contig_out = open(opts.read_on_contig, 'wb')
#    with open(dirpath + "/out.readOnContig.gz", mode='rb') as file: # b is important -> binary
#        fileContent = file.read()
#        read_on_contig_out.write(fileContent)
#    read_on_contig_out.close()
#    file.close()

    #Set up command line call
    #TODO - remove hard coded path
    #Code for adding directory path to other file required as output
    if opts.default_full_settings_type == "default":
        cmd = "SOAPdenovo-63mer_v2.0 scaff -g %s -F" % (dirpath + "/out")
    elif opts.default_full_settings_type == "full":
        cmd = "SOAPdenovo-63mer_v2.0 scaff -g %s -F %s -z %s -u %s -w %s -v %s -G %s -L %s -c %s -C %s -b %s -B %s -N %s -p %s" % (dirpath + "/out",
            opts.fill_gaps, opts.compatible_mode, opts.unmask_contigs, opts.keep_contigs_connected, opts.ass_visual, opts.gap_len_diff,
            opts.min_contig_len, opts.min_contig_cvg, opts.max_contig_cvg, opts.insert_size_upper_bound, opts.bubble_coverage,
            opts.genome_size, opts.ncpu)
        #Check
        
    print cmd


    #Perform SOAPdenovo2_config analysis
    buffsize = 1048576
    #Create temp directory for standard error and out
    tmp_dir = tempfile.mkdtemp()
    try:

        tmp_out_file = tempfile.NamedTemporaryFile(dir=tmp_dir).name
        tmp_stdout = open(tmp_out_file, 'wb')
        tmp_err_file = tempfile.NamedTemporaryFile(dir=tmp_dir).name
        tmp_stderr = open(tmp_err_file, 'wb')

        #Call SOAPdenovo2
        #New additional datasets must be placed in the directory provided by $__new_file_path__
        proc = subprocess.Popen(args=cmd, shell=True, cwd=dirpath, stdout=tmp_stdout, stderr=tmp_stderr.fileno())
        returncode = proc.wait()
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
        raise Exception, 'Problem performing scaff process ' + str(e)

    #Read soap config file into its output
    new_contig_index_out = open(opts.new_contig_index, 'wb')
    file = open(dirpath + "/out.newContigIndex")
    for line in file:
        new_contig_index_out.write(line)
    new_contig_index_out.close()
    file.close()

    links_out = open(opts.links, 'wb')
    file = open(dirpath + "/out.links")
    for line in file:
        links_out.write(line)
    links_out.close()
    file.close()

    scaf_gap_out = open(opts.scaf_gap, 'wb')
    file = open(dirpath + "/out.scaf_gap")
    for line in file:
        scaf_gap_out.write(line)
    scaf_gap_out.close()
    file.close()

    gap_seq_out = open(opts.gap_seq, 'wb')
    file = open(dirpath + "/out.gapSeq")
    for line in file:
        gap_seq_out.write(line)
    gap_seq_out.close()
    file.close()

    scaf_out = open(opts.scaf, 'wb')
    file = open(dirpath + "/out.scaf")
    for line in file:
        scaf_out.write(line)
    scaf_out.close()
    file.close()

    scaf_seq_out = open(opts.scaf_seq, 'wb')
    file = open(dirpath + "/out.scafSeq")
    for line in file:
        scaf_seq_out.write(line)
    scaf_seq_out.close()
    file.close()

    contig_positions_scaff_out = open(opts.contig_positions_scaff, 'wb')
    file = open(dirpath + "/out.contigPosInscaff")
    for line in file:
        contig_positions_scaff_out.write(line)
    contig_positions_scaff_out.close()
    file.close()

    bubble_in_scaff_out = open(opts.bubble_in_scaff, 'wb')
    file = open(dirpath + "/out.bubbleInScaff")
    for line in file:
        bubble_in_scaff_out.write(line)
    bubble_in_scaff_out.close()
    file.close()

    scaf_stats_out = open(opts.scaf_stats, 'wb')
    file = open(dirpath + "/out.scafStatistics")
    for line in file:
        scaf_stats_out.write(line)
    scaf_stats_out.close()
    file.close()

    #Clean up temp files
    cleanup_before_exit(tmp_dir)
    #Check results in output file
    if os.path.getsize(opts.scaf_stats) > 0:
        sys.stdout.write('Status complete')
    else:
        stop_err("The output is empty")

if __name__ == "__main__": main()
