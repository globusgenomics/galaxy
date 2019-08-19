import os, sys
infile = sys.argv[1]

fh = open(infile, "r")
header = fh.readline()
fusion_pairs = {}
fusions =  []
header_vals = ['FUSION_ID', 'SPANNING_READS', 'ENCOMPASSING_READS', 'GENOMIC', '5_FPG_GENE_NAME', '5_IN_CDS', '5_SEGMENT_TYPE', '3_FPG_GENE_NAME', '3_IN_CDS', '3_SEGMENT_TYPE', 'P_VAL_CORR', 'DRIVER_PROB', 'EXPRESSION_GAIN']
print "\t".join(header_vals)
for line in fh:
    values = line.split("\t")
    fusion_id = "%s-%s" % (values[6],values[13])
    if fusion_id in fusion_pairs:
        #print line
        fusion_pairs[fusion_id]['SPANNING_READS'] += int(values[3])
        fusion_pairs[fusion_id]['ENCOMPASSING_READS'] += int(values[4])
    else:
        fusions.append(fusion_id)
        fusion_pairs[fusion_id] = {}
        fusion_pairs[fusion_id]['FUSION_ID'] = values[1]
        fusion_pairs[fusion_id]['SPANNING_READS'] = int(values[3])
        fusion_pairs[fusion_id]['ENCOMPASSING_READS'] = int(values[4])
        fusion_pairs[fusion_id]['GENOMIC'] = values[5]
        fusion_pairs[fusion_id]['5_FPG_GENE_NAME'] = values[6]
        fusion_pairs[fusion_id]['5_IN_CDS'] = values[7]
        fusion_pairs[fusion_id]['5_SEGMENT_TYPE'] = values[8]
        fusion_pairs[fusion_id]['3_FPG_GENE_NAME'] = values[13]
        fusion_pairs[fusion_id]['3_IN_CDS'] = values[14]
        fusion_pairs[fusion_id]['3_SEGMENT_TYPE'] = values[15]
        fusion_pairs[fusion_id]['P_VAL_CORR'] = values[21]
        fusion_pairs[fusion_id]['DRIVER_PROB'] = values[22]
        fusion_pairs[fusion_id]['EXPRESSION_GAIN'] = values[23]

for i in fusions:
    print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (i, str(fusion_pairs[i]['SPANNING_READS']), str(fusion_pairs[i]['ENCOMPASSING_READS']), fusion_pairs[i]['GENOMIC'], fusion_pairs[i]['5_FPG_GENE_NAME'], fusion_pairs[i]['5_IN_CDS'], fusion_pairs[i]['5_SEGMENT_TYPE'], fusion_pairs[i]['3_FPG_GENE_NAME'], fusion_pairs[i]['3_IN_CDS'], fusion_pairs[i]['3_SEGMENT_TYPE'], fusion_pairs[i]['P_VAL_CORR'], fusion_pairs[i]['DRIVER_PROB'], fusion_pairs[i]['EXPRESSION_GAIN'])
