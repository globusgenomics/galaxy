#!/usr/bin/python
import os, sys, shutil, glob

def get_first_fasta_id(fasta_path):
    fh = open(fasta_path, "r")
    chrm = None
    for line in fh:
        ##sequence-region GG745328 1 3145054
        if line.startswith(">"):
            values = line.split(" ")
            chrm = values[0]
            values = chrm.split(">")
            chrm = values[1]
            return chrm
    fh.close()

def get_first_chr(gff_path):
    fh = open(gff_path, "r")
    chrm = None
    for line in fh:
        ##sequence-region GG745328 1 3145054
        if line.startswith("##sequence-region"):
            values = line.split(" ")
            chrm = values[1]
            return chrm
    fh.close()

    if chrm is None:
        fh = open(gff_path, "r")
        for line in fh:
            if not line.startswith("#"):
                values = line.split("\t")
                return values[0]

path = sys.argv[1]
items = {}
display_names = {}
for dbkey in os.listdir(path):
    ffs = glob.glob("%s/%s/seq/*.fasta" % (path, dbkey))
    if len(ffs) > 0:
        ff =  ffs[0]
        chrm = get_first_fasta_id(ff)
        values = dbkey.split("-")
        dbname = values[0]
        db = dbname.lower()
        if db in items:
            items[db][dbkey] = chrm
        else:
            items[db] = {dbkey:chrm}
            display_names[db] = "%s GBrowse" % dbname

for db in items:
    current_line = "%s\t" % db
    display_name = display_names[db]
    dbkey_field = items[db].keys()
    chrom_field = items[db].values()
    #print db
    #print dbkey_field
    #print chrom_field
    print "%s\t%s\t%s\t%s" % (db, display_name, ",".join(dbkey_field), ",".join(chrom_field))
