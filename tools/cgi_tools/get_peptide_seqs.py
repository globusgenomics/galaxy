#!/usr/bin/python
import requests, sys, argparse
requests.packages.urllib3.disable_warnings()

def get_list_from_annotations(annotations):
    locs = {}
    for line in annotations:
        if line.startswith("#"):
            continue
        if "missense_variant" not in line:
            continue
        if "SYMBOL_SOURCE=HGNC" not in line:
            continue
        values = line.rstrip("\n").split("\t")
        extras = values[13].split(";")
        meta = {}
        for i in extras:
            key,value = i.split("=")
            meta[key] = value
        # only get the value for SYMBOL_SOURCE=HGNC
        # get the ENSP accession value

        ref_aa, allele_aa = values[10].split("/")
        if values[1] not in locs:
            locs[values[1]] = [{'gene_name' : meta['SYMBOL'], 'accesion' : meta['ENSP'],
                                'ref_aa' : ref_aa, 'meta' : meta,
                                'allele_aa' : allele_aa, 'aa_pos' : values[9]}]
        else:
            locs[values[1]].append({'gene_name' : meta['SYMBOL'], 'accesion' : meta['ENSP'],
                                    'ref_aa' : ref_aa, 'meta' : meta,
                                    'allele_aa' : allele_aa, 'aa_pos' : values[9]})

    return locs

def get_list_from_targets(targets):
    locs = {}
    for line in targets:
        (chrm, nuc_pos, prot_pos, change, gene_name) = line.rstrip("\n").split("\t")
        location = "%s:%s" % (chrm, nuc_pos)
        ref_aa, allele_aa = change.split("/")
        if location not in locs:
            locs[locations] = [{'gene_name' : gene_name, 'accession' : '', 'ref_aa' : ref_aa,
                               'allele_aa' : allele_aa, 'aa_pos' : prot_pos, 'meta' : {}}]
    return locs

def __main__():
    parser = argparse.ArgumentParser(description='Parse variant list and get peptide sequences by making API calls to the uniprot service.')
    parser.add_argument('-t', '--targets', type=argparse.FileType('r'), required=False,
                    help='File containing the positions and protein name to get')
    parser.add_argument('-a', '--annotations', type=argparse.FileType('r'), required=False,
                    help='File containing the annotated positions of the final VCF')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), required=True,
                    help='Output file to write the peptide sequences')
    args = parser.parse_args()

    #server = "http://rest.ensembl.org/sequence/id/"
    server = "http://grch37.rest.ensembl.org/sequence/id/"
    target_locations = {}
    if args.annotations:
        target_locations = get_list_from_annotations(args.annotations)
    else:
        target_locations = get_list_from_targets(args.targets)

    seen = {}
    for loc in target_locations:
        for target in target_locations[loc]:
            sequence = None
            accession = target['accesion']
            masked_prot_loc = int(target['aa_pos']) - 1
            wt_name = "%s.%s.%s%s%s" % ("WT", target['meta']['SYMBOL'], target['ref_aa'], target['aa_pos'], target['allele_aa'])
            mt_name = "%s.%s.%s%s%s" % ("MT", target['meta']['SYMBOL'], target['ref_aa'], target['aa_pos'], target['allele_aa'])
            if mt_name in seen:
                continue
            seen[mt_name] = 0
            url = "%s%s?" % (server, target['accesion'])
            r = requests.get(url, headers={ "Accept" : "application/json"})
            if not r.ok:
                accession = None
                if 'SWISSPROT' in target['meta']:
                    accession = target['meta']['SWISSPROT'].split(",")[0]
                elif 'TREMBL' in target['meta']:
                    accession = target['meta']['TREMBL'].split(",")[0]
                
                if accession is not None:
                    url = "https://www.ebi.ac.uk/proteins/api/proteins/%s" % accession
                    r = requests.get(url, headers={ "Accept" : "application/json"})
                    if not r.ok:
                        #args.output.write( ">" + wt_name + "\n")
                        #args.output.write( "NONE-%s" % target['accesion'] + "\n")
                        #args.output.write( ">" + mt_name + "\n" )
                        #args.output.write( "NONE-%s" % target['accesion'] + "\n")
                        continue
                    sequence = r.json()['sequence']
                else:
                    continue
            else:
                sequence = r.json()['seq']

            if len(sequence) > masked_prot_loc and sequence[masked_prot_loc] == target['ref_aa']:
                wt_seq = sequence
                mt_seq = sequence[:masked_prot_loc] + target['allele_aa'] + sequence[masked_prot_loc+1:]
                aa_change_loc = None
                if len(sequence) - (masked_prot_loc) < 9:
                    extra_length = 10 + (11 - (len(sequence) - (masked_prot_loc)))
                    wt_pep = wt_seq[masked_prot_loc-extra_length:len(sequence)]
                    mt_pep = mt_seq[masked_prot_loc-extra_length:len(sequence)]
                    aa_change_loc = len(wt_pep) - (len(sequence)-masked_prot_loc) + 1
                elif masked_prot_loc < 9:
                    extra_length = 11 + (10 - masked_prot_loc)
                    wt_pep = wt_seq[0:masked_prot_loc+extra_length]
                    mt_pep = mt_seq[0:masked_prot_loc+extra_length]
                    aa_change_loc = masked_prot_loc+1
                else:
                    wt_pep = wt_seq[masked_prot_loc-10:masked_prot_loc+11]
                    mt_pep = mt_seq[masked_prot_loc-10:masked_prot_loc+11]
                    aa_change_loc = 11

                print "%s\t%s\t%s\t%s\t%s\t%s" % (accession, wt_name, aa_change_loc, loc, target['gene_name'], wt_pep)
                print "%s\t%s\t%s\t%s\t%s\t%s" % (accession, mt_name, aa_change_loc, loc, target['gene_name'], mt_pep)
                args.output.write( ">" + wt_name +"\n")
                args.output.write(wt_pep +"\n")
                args.output.write( ">" + mt_name+"\n")
                args.output.write(mt_pep+"\n")

if __name__=="__main__": __main__()
