#!/usr/bin/python
import sys, argparse

def __main__():
    parser = argparse.ArgumentParser(description='Filter expression profiling')
    parser.add_argument('-p', '--peptides', type=argparse.FileType('r'), required=True,
                    help='File containing the peptides list')
    parser.add_argument('-t', '--threshold', type=float, help='FPKM threshold')
    parser.add_argument('-n', '--index', type=argparse.FileType('r'), required=True,
                    help='Peptides Fasta index information')
    parser.add_argument('-m', '--netmhc', required=True,
                    help='NetMHC raw output') 
    #parser.add_argument('-s', '--sites', required=True, type=argparse.FileType('w'),
    #                help='sites output file')
    args = parser.parse_args()

    # parse list of final neoantigen candidates
    peptides = {}
    seen = {}
    for line in args.peptides:
        values = line.rstrip("\n").split("\t")
        if values[15] in seen:
            continue
        seen[values[15]] = 1 
        #sites.write("%s\t%s\t%s" % (values[15].split(":")[0], values[15].split(":")[1], values[15].split(":")[1]))
        peptides[values[10]] = {'pos_index' : values[0], 'hla_type' : values[1],
                                'mt_peptide' : values[2], 'core' : values[3], 'offset' : values[4],
                                'ipos' : values[5], 'ilen' : values[6], 'dpos' : values[7],
                                'dlen' : values[8], 'icore' : values[9], 'aff_log' : values[11],
                                'aff_nm_mutant' : values[12], 'fpkm' : values[14], 'loc' : values[15]}

    # use the index file of the peptide fasta file to check where the variant is
    for line in args.index:
        values = line.rstrip("\n").split("\t")
        name = values[1].replace(".", "_")[0:15]
        if name not in peptides:
            continue
        peptides[name]['wt_name'] = name.replace("MT_", "WT_")
        peptides[name]['ensp'] = values[0]
        peptides[name]['fasta_seq_name'] = values[1]
        peptides[name]['aa_pos_change'] = values[2]
        peptides[name]['gene_name'] = values[4]
        (mt_wt, gene, variant) = peptides[name]['fasta_seq_name'].split(".")
        peptides[name]['aa_wt'] = variant[0]
        peptides[name]['aa_mt'] = variant[-1]
        mt_sequence = values[-1]
        wt_sequence = mt_sequence[:int(values[2])-1] + peptides[name]['aa_wt'] + mt_sequence[int(values[2]):]
        peptides[name]['wt_peptide'] = wt_sequence[int(peptides[name]['pos_index']):int(peptides[name]['pos_index'])+len(peptides[name]['mt_peptide'])]

        # get the wildtype affinity value by looking at the raw NetMHC output
        for mhc_line in open(args.netmhc,"r"):
            if peptides[name]['wt_peptide'] in mhc_line and \
               peptides[name]['wt_name'] in mhc_line and peptides[name]['hla_type'] in mhc_line:
                mhc_values = mhc_line.rstrip("\n").split()
                if len(mhc_values[2]) == len(peptides[name]['wt_peptide']):
                    peptides[name]['aff_nm_wildtype'] = mhc_values[12]
                    peptides[name]['fold_change'] = float(peptides[name]['aff_nm_wildtype']) / float(peptides[name]['aff_nm_mutant'])
                    break

    # print final information
    print "SEQ_NAME\tNAME\tENSP\tNEOANTIGEN\tMUTANT AFFINITY (NM)\tWT AFFINITY (NM)\tAFFINITY FOLD CHANGE (WT/MT)\tFPKM\tLOCATION"
    for name, value in sorted(peptides.iteritems(), key=lambda (x, y): y['fold_change'], reverse=True):
        print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (value['fasta_seq_name'], name, value['ensp'], value['mt_peptide'], value['aff_nm_mutant'], value['aff_nm_wildtype'], value['fold_change'], value['fpkm'], value['loc'])

if __name__=="__main__": __main__()
