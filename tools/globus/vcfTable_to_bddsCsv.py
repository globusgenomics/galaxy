#!/usr/bin/python

import sys, os, csv
from datetime import datetime
from optparse import OptionParser

def __main__():

    usage = "--input <input_table> --output <output_file>"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", "--input", dest="input_file", help="input file")
    parser.add_option("-o", "--output", dest="output_file", help="output file")
    options, args = parser.parse_args()

    csv_in = open(options.input_file, 'r')
    csv_out = open(options.output_file, 'w')
    reader = csv.DictReader(csv_in, delimiter='\t')
    header = list(['RID', 'PID'])
    header.extend(list(fn for fn in reader.fieldnames))
    old_header =  list(fn for fn in reader.fieldnames)
    pid = header[-1].split(".")[0].split("_")[-1]
    i = 0
    for col in header:
        if pid in col:
            col = col.split(".")[-1]
        header[i] = col
        i += 1

    #del(header[-1])  # remove patient specific last column name
    writer = csv.writer(csv_out, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(header)  # for 2.6 support since no writeheader
    rownum = 1
    for row in reader:
        dt = datetime.now()
        rowid = "%s_%s_%s" % (pid, rownum, dt.isoformat("T"))
        newrow = list()
        newrow.append(rowid)
        newrow.append(pid)  # put patient ID here
        for col in old_header:
            newrow.append(row[col])
        #newrow.extend(row.values())
        writer.writerow(newrow)
        rownum += 1
    csv_in.close()
    csv_out.close()

if __name__=="__main__": __main__()

