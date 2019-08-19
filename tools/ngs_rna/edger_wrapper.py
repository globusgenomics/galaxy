#!/usr/bin/env python

import argparse
import tempfile
import subprocess
import sys
import os

def append_dge_file(gene_counts_files, group, sample_name):
    for i, f in enumerate(gene_counts_files):
        tmp_dge_fh.write("\t".join([f, group, sample_name[i]]) + "\n")
    return

#
# Parse command line
#
parser = argparse.ArgumentParser(description='EdgeR wrapper')

# Input files
parser.add_argument(
                    '--group1-name',
                    '-g1',
                    required=True,
                    help='Name for group #1',
                   )
parser.add_argument(
                    '--group2-name',
                    '-g2',
                    required=True,
                    help='Name for group #2',
                   )
parser.add_argument(
                    '--group1-gene-counts',
                    '-gc1',
                    required=True,
                    action='append',
                    help='Gene counts file for group #1',
                   )
parser.add_argument(
                    '--group2-gene-counts',
                    '-gc2',
                    required=True,
                    action='append',
                    help='Gene counts file for group #2',
                   )
parser.add_argument(
                    '--group1-sample-name',
                    '-gsn1',
                    required=True,
                    action='append',
                    help='Name of group 1 sample name',
                   )
parser.add_argument(
                    '--group2-sample-name',
                    '-gsn2',
                    required=True,
                    action='append',
                    help='Name of group 2 sample name',
                   )

# Output files
parser.add_argument(
                    '--mds-plot-output',
                    '-m',
                    required=True,
                    help='Name for MDS plot output file (PNG format)',
                   )
parser.add_argument(
                    '--tagwise-disp-plot-output',
                    '-t',
                    required=True,
                    help='Name for tagwise dispersion plot output file (PNG format)',
                   )
parser.add_argument(
                    '--ma-plot-output',
                    '-a',
                    required=True,
                    help='Name for MA plot output file (PNG format)',
                   )
parser.add_argument(
                    '--diff-expr-output',
                    '-d',
                    required=True,
                    help='Name for differential expression results output file (TSV format)',
                   )
parser.add_argument(
                    '--sig-diff-expr-output',
                    '-s',
                    required=True,
                    help='Name for significant differential expression results output file (TSV format)',
                   )

args = parser.parse_args()
# Create DGE file
# TODO not make tmp, but output file?
tmp_dge_obj = tempfile.NamedTemporaryFile(delete=False)
tmp_dge = tmp_dge_obj.name
tmp_dge_fh = open(tmp_dge, 'w')

tmp_dge_fh.write("\t".join(['files','group','description']) + "\n")
print args.group1_sample_name
append_dge_file(args.group1_gene_counts, args.group1_name, args.group1_sample_name)
append_dge_file(args.group2_gene_counts, args.group2_name, args.group2_sample_name)
tmp_dge_fh.close()



# Build EdgeR command
cmd = 'Rscript /opt/galaxy/tools/ngs_rna/edger_wrapper.R %s %s %s %s %s %s' % (tmp_dge, args.mds_plot_output, args.tagwise_disp_plot_output, args.ma_plot_output, args.diff_expr_output, args.sig_diff_expr_output)
print cmd
proc = subprocess.Popen(cmd, shell=True)
returncode = proc.wait()

# Clean up
#os.remove(tmp_dge)


