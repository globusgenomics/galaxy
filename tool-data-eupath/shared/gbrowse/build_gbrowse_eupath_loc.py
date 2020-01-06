#!/usr/bin/python
import sys, os, shutil

fh = open(sys.argv[1], "r")
items = {}
for line in fh:
    (db, dbkey, chrname) = line.split("\t")
    if db in items:
        items[db][dbkey] = chrname.split(",")[0]
    else:
        items[db] = { dbkey:chrname.split(",")[0]}

for db in items:
    current_line = "%s\t" % db
    dbkey_field = items[db].keys()
    chrom_field = items[db].values()
    print "%s\t%s\t%s" % (db, ",".join(dbkey_field), ",".join(chrom_field))
        
