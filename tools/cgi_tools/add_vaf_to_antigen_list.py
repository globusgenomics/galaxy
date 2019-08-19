#!/usr/bin/python
import sys, argparse

def __main__():
    parser = argparse.ArgumentParser(description='Filter expression profiling')
    parser.add_argument('-p', '--peptides', type=argparse.FileType('r'), required=True,
                    help='File containing the peptides list')
    parser.add_argument('-v', '--vcf', type=argparse.FileType('r'), required=True,
                    help='Consensus VCF file')
    args = parser.parse_args()

    # parse list of final neoantigen candidates
    peptides = {}
    order = []
    header = args.peptides.readline()
    for line in args.peptides:
        values = line.rstrip("\n").split("\t")
        peptides[values[-1]] = values
        order.append(values[-1])

    # use the vcf file to get the AVF value for each neoantigen candidate
    vafs = {}
    for line in args.vcf:
        if line.startswith("#"):
            continue
        values = line.rstrip("\n").split("\t")
        location = "%s:%s" % (values[0], values[1])
        if location in peptides:
            format_values = values[-1].split(":")
            #ad = format_values[5]  # freebayes
            #total = format_values[3] # freebayes
            #vaf = float(float(ad)/float(total)) # freebayes
            ad = format_values[1] # mutect2
            vaf = format_values[2].split(",")[1] # mutect2
            vafs[location] = [str(vaf), str(ad)]
            continue

    # print final information
    print "SEQ_NAME\tNAME\tENSP\tNEOANTIGEN\tMUTANT AFFINITY (NM)\tWT AFFINITY (NM)\tNM_FOLD_CHANGE\tFPKM\tLOCATION\tVAF\tAD"
    for location in order:
        print "%s\t%s" % ("\t".join(peptides[location]), "\t".join(vafs[location]))

if __name__=="__main__": __main__()
