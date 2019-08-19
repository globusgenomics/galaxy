#!/usr/bin/python
import sys, argparse

def __main__():
    parser = argparse.ArgumentParser(description='Filter expression profiling')
    parser.add_argument('-p', '--peptides', type=argparse.FileType('r'), required=True,
                    help='File containing the peptides list')
    parser.add_argument('-t', '--threshold', type=float, help='FPKM threshold')
    parser.add_argument('-n', '--index', type=argparse.FileType('r'), required=True,
                    help='Peptides Fasta index information')
    parser.add_argument('-c', '--cufflinks', type=argparse.FileType('r'), required=True,
                    help='File containing cufflinks FPKM values')
    args = parser.parse_args()

    # use the index file of the peptide fasta file to check where the variant is
    index = {}
    locs = []
    for line in args.index:
        values = line.rstrip("\n").split("\t")
        name = values[1].replace(".", "_")[0:15]
        index[name] = {'loc': values[3]}
        locs.append(values[3])

    fpkm = {}
    for line in args.cufflinks:
        values = line.rstrip("\n").split("\t")
        if values[2] != "exon":
            continue
        chrm = values[0]
        start = values[3]
        end = values[4]
        fpkm_val = values[8].split(";")[3].split("\"")[-2]
        if float(fpkm_val) >= 1:
            # make sure this region is asked for
            for loc in locs:
                if loc in fpkm:
                    continue
                loc_chrm, loc_pos = loc.split(":")
                if loc_chrm != chrm:
                    continue
                if int(loc_pos) >= int(start) and int(loc_pos) <= int(end):
                    fpkm[loc] = fpkm_val
                    break

    peptides = {}
    for line in args.peptides:
        values = line.rstrip("\n").split("\t")
        name = values[10]
        if index[name]['loc'] in fpkm:
            values.append(fpkm[index[name]['loc']])
            values.append(index[name]['loc'])
            print "\t".join(values) 

if __name__=="__main__": __main__()
