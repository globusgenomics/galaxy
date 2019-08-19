#!/usr/bin/python
import sys, argparse

def __main__():
    parser = argparse.ArgumentParser(description='Parse netMHC output')
    parser.add_argument('-i', '--input', type=argparse.FileType('r'), required=True,
                    help='File containing the netMHC hits')
    parser.add_argument('-t', '--threshold', type=float, help='Affinity threshold')
    parser.add_argument('-n', '--index', type=argparse.FileType('r'), required=True,
                    help='Peptides Fasta index information')
    args = parser.parse_args()

    # use the index file of the peptide fasta file to check where the variant is
    index = {}
    for line in args.index:
        values = line.rstrip("\n").split("\t")
        name = values[1].replace(".", "_")[0:15]
        index[name] = {'mut_loc': int(values[2])-1}

    candidates = {}
    seen = {}
    for line in args.input:
        if line.startswith("#") or line.startswith("\n") or \
           line.startswith("--") or line.startswith("  pos") or \
           line.startswith("Protein"):
            continue
        if line.rstrip("\n") in seen:
            continue
        seen[line.rstrip("\n")] = 0
        values = line.rstrip("\n").split()

        # keep only hits below the threshold (500)
        if float(values[12]) <= args.threshold:
            #print line
            #print values[12]
            if values[10] in candidates:
                candidates[values[10]].append(values)
            else:
                candidates[values[10]] = [values]

    # pick the best binder for each variant that contains the aa variant
    final_candidates = {}
    for candy in candidates:
        if "WT_" in candy:
            continue
        lowest = 500.1
        for hit in candidates[candy]:
            # make sure aa variant is in the hit
            mut_loc = int(index[hit[10]]['mut_loc'])
            pep_start = int(hit[0])
            pep_end = int(int(pep_start) + len(hit[2]))
            if mut_loc >= pep_start and mut_loc <= pep_end:
                if float(hit[12]) <= lowest:
                    final_candidates[candy] = hit
                    lowest = hit[12]

    for candy in final_candidates:
        print "\t".join(final_candidates[candy][0:14])

if __name__=="__main__": __main__()
