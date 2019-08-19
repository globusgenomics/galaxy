#!/usr/bin/env python
"""Select FASTA sequences by IDs from a BLAST output tabular file.

Takes five command line options, tabular filename, ID column number (using
one based counting), input filename, input type (e.g. FASTA) and the
output filename (same format as input sequence file).

This tool is a short Python script which requires Biopython 1.54 or later
for SFF file support. If you use this tool in scientific work leading to a
publication, please cite the Biopython application note:

Cock et al 2009. Biopython: freely available Python tools for computational
molecular biology and bioinformatics. Bioinformatics 25(11) 1422-3.
http://dx.doi.org/10.1093/bioinformatics/btp163 pmid:19304878.

"""
import sys

def stop_err(msg, err=1):
    sys.stderr.write(msg.rstrip() + "\n")
    sys.exit(err)

if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.6"
    sys.exit(0)

#Parse Command Line
try:
    tabular_file, col_arg, in_file, seq_format, out_file = sys.argv[1:]
except ValueError:
    stop_err("Expected five arguments, got %i:\n%s" % (len(sys.argv)-1, " ".join(sys.argv)))
try:
    if col_arg.startswith("c"):
        column = int(col_arg[1:])-1
    else:
        column = int(col_arg)-1
except ValueError:
    stop_err("Expected column number, got %s" % col_arg)

if seq_format.lower() in ["sff", "fastq", "qual", "fasta"]:
    seq_format = seq_format.lower()
else:
    stop_err("Unrecognised file format %r" % seq_format)


try:
    from Bio import SeqIO
except ImportError:
    stop_err("Biopython 1.54 or later is required")


def parse_ids(tabular_file, col):
    """Read tabular file and record all specified identifiers."""
    handle = open(tabular_file, "rU")
    for line in handle:
        if line.strip() and not line.startswith("#"):
            yield line.rstrip("\n").split("\t")[col].strip(), line.rstrip("\n").split("\t")[col-1].strip()
    handle.close()

#Index the sequence file.
#If very big, could use SeqIO.index_db() to avoid memory bottleneck...
records = SeqIO.index(in_file, seq_format)
print "Indexed %i sequences" % len(records)

#Avoid overhead of parsing into SeqRecord objects,
#just re-use the original formatting from the input file.
out_handle = open(out_file, "w")
count = 0
for name, query in parse_ids(tabular_file, column):
    try:
        description = records.get_raw(name).split('\n')[0]
        values = description.split(" ")
        id_dict = dict()
        for value in values:
            if ":" in value:
                new_value = value.replace("[", "")
                id_sections = new_value.split(":")
                db_name = id_sections[0]
                db_value = id_sections[-1]
                if db_name in id_dict:
                    id_dict[db_name].append(db_value)
                else:
                    id_dict[db_name] = [db_value]
        ## print the output
        for db, value_list in id_dict.iteritems():
            value_list.sort()
            uniq_list = list(set(value_list))
            for id_value in uniq_list:
                out_handle.write("%s\t%s\t%s\t%s\n" % (query, name, db, id_value) )
    except KeyError:
        out_handle.close()
        stop_err("Identifier %r not found in sequence file" % name)
    count += 1
out_handle.close()

#print "Selected %i sequences by ID" % count
