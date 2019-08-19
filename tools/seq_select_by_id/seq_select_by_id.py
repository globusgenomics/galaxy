#!/usr/bin/env python
"""Select FASTA, QUAL, FASTQ or SSF sequences by IDs from a tabular file.

Takes five command line options, tabular filename, ID column number (using
one based counting), input filename, input type (e.g. FASTA or SFF) and the
output filename (same format as input sequence file).

When selecting from an SFF file, any Roche XML manifest in the input file is
preserved in both output files.

This tool is a short Python script which requires Biopython 1.54 or later
for SFF file support. If you use this tool in scientific work leading to a
publication, please cite the Biopython application note:

Cock et al 2009. Biopython: freely available Python tools for computational
molecular biology and bioinformatics. Bioinformatics 25(11) 1422-3.
http://dx.doi.org/10.1093/bioinformatics/btp163 pmid:19304878.

This script is copyright 2011-2013 by Peter Cock, The James Hutton Institute UK.
All rights reserved. See accompanying text file for licence details (MIT
license).

This is version 0.0.6 of the script.
"""
import sys, os

def stop_err(msg, err=1):
    sys.stderr.write(msg.rstrip() + "\n")
    sys.exit(err)

if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.6"
    sys.exit(0)

galhtmlprefix = """<?xml version="1.0" encoding="utf-8" ?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<meta name="generator" content="Galaxy tool output - see http://getgalaxy.org/" />
<title></title>
<link rel="stylesheet" href="/static/style/base.css" type="text/css" />
</head>
<body>
<div class="document">
"""
galhtmlpostfix = """</div></body></html>\n"""

#Parse Command Line
try:
    tabular_file, col_arg, in_file, seq_format, out_file, group_col, extra_files_path = sys.argv[1:]
except ValueError:
    stop_err("Expected seven arguments, got %i:\n%s" % (len(sys.argv)-1, " ".join(sys.argv)))
try:
    if col_arg.startswith("c"):
        column = int(col_arg[1:])-1
    else:
        column = int(col_arg)-1
except ValueError:
    stop_err("Expected column number, got %s" % col_arg)

try:
    if group_col.startswith("c"):
        group_column = int(group_col[1:])-1
    else:
        group_column = int(group_col)-1
except ValueError:
    stop_err("Expected column number, got %s" % col_arg)


if seq_format == "fastqcssanger":
    stop_err("Colorspace FASTQ not supported.")
elif seq_format.lower() in ["sff", "fastq", "qual", "fasta"]:
    seq_format = seq_format.lower()
elif seq_format.lower().startswith("fastq"):
    #We don't care how the qualities are encoded    
    seq_format = "fastq"
elif seq_format.lower().startswith("qual"):
    #We don't care what the scores are
    seq_format = "qual"
else:
    stop_err("Unrecognised file format %r" % seq_format)

if group_col != "None":
    os.mkdir(extra_files_path)

try:
    from Bio import SeqIO
except ImportError:
    stop_err("Biopython 1.54 or later is required")


def parse_ids(tabular_file, col, group_col="None"):
    """Read tabular file and record all specified identifiers."""
    handle = open(tabular_file, "rU")
    for line in handle:
        if line.strip() and not line.startswith("#"):
            if group_col=="None":
                yield line.rstrip("\n").split("\t")[col].strip()
            else:
                yield line.rstrip("\n").split("\t")[col].strip(), line.rstrip("\n").split("\t")[group_col].strip()
    handle.close()

#Index the sequence file.
#If very big, could use SeqIO.index_db() to avoid memory bottleneck...
records = SeqIO.index(in_file, seq_format)
print "Indexed %i sequences" % len(records)

if seq_format.lower()=="sff":
    #Special case to try to preserve the XML manifest
    try:
        from Bio.SeqIO.SffIO import SffIterator, SffWriter
    except ImportError:
        stop_err("Requires Biopython 1.54 or later")

    try:
        from Bio.SeqIO.SffIO import ReadRocheXmlManifest
    except ImportError:
        #Prior to Biopython 1.56 this was a private function
        from Bio.SeqIO.SffIO import _sff_read_roche_index_xml as ReadRocheXmlManifest

    in_handle = open(in_file, "rb") #must be binary mode!
    try:
        manifest = ReadRocheXmlManifest(in_handle)
    except ValueError:
        manifest = None
    in_handle.close()

    out_handle = open(out_file, "wb")
    writer = SffWriter(out_handle, xml=manifest)
    count = 0
    #This does have the overhead of parsing into SeqRecord objects,
    #but doing the header and index at the low level is too fidly.
    iterator = (records[name] for name in parse_ids(tabular_file, column))
    try:
        count = writer.write_file(iterator)
    except KeyError, err:
        out_handle.close()
        if name not in records:
            stop_err("Identifier %r not found in sequence file" % name)
        else:
            raise err
    out_handle.close()
else:
    #Avoid overhead of parsing into SeqRecord objects,
    #just re-use the original formatting from the input file.
    out_handle = open(out_file, "w")
    count = 0
    if group_col == "None":
        for name in parse_ids(tabular_file, column):
            try:
                out_handle.write(records.get_raw(name))
            except KeyError:
                out_handle.close()
                stop_err("Identifier %r not found in sequence file" % name)
            count += 1
    else:
        seqs_by_group = dict()
        out_handle.write(galhtmlprefix)
        for name, group in parse_ids(tabular_file, column, group_column):
            # push reads to dictionary of lists
            if group in seqs_by_group:
                seqs_by_group[group].append(records.get_raw(name))
            else:
                seqs_by_group[group] = list(records.get_raw(name))
            count += 1
        ## write output for each group to a file
        for group, seq_list in seqs_by_group.iteritems():
            filename1 = group.replace (" ", "_")
            filename = filename1.replace ("/", "_")
            group_filename = "%s/%s.fasta" % (extra_files_path, filename)
            group_handle = open(group_filename, "w")
            group_handle.write("".join(seq_list))
            group_handle.close()
            out_handle.write('<a href="%s.fasta">%s.fasta</a><br>\n' % (filename, filename))
        out_handle.write(galhtmlpostfix)
    out_handle.close()

print "Selected %i sequences by ID" % count
