import sys
import re

def extract_vcf_variant (strand, pos, ref, alt):
    pos = int(pos)
    reflen = len(ref)
    altlen = len(alt)
    minlen = min(reflen, altlen)
    new_ref = ref
    new_alt = alt

    if reflen == 1 and altlen == 1 and ref == alt:
        return pos, ref, alt

    for nt_pos in xrange(0, minlen):
        if ref[reflen - nt_pos - 1] == alt[altlen - nt_pos - 1]:
            new_ref = ref[:reflen - nt_pos - 1]
            new_alt = alt[:altlen - nt_pos - 1]
        else:
            break
    newreflen = len(new_ref)
    newaltlen = len(new_alt)

    minlen = min(newreflen, newaltlen)
    new_pos = pos
    new_ref2 = new_ref
    new_alt2 = new_alt

    for nt_pos in xrange(minlen):
        if new_ref[nt_pos] == new_alt[nt_pos]:
            if strand == '+':
                new_pos += 1
            elif strand == '-':
                new_pos -= 1
            new_ref2 = new_ref[nt_pos + 1:]
            new_alt2 = new_alt[nt_pos + 1:]
        else:
            new_ref2 = new_ref[nt_pos:]
            new_alt2 = new_alt[nt_pos:]
            break
    if new_ref == '':
        new_ref2 = '-'
    if new_alt2 == '':
        new_alt2 = '-'

    return new_pos, new_ref2, new_alt2

input_filename = sys.argv[1]
output_filename = sys.argv[2]

f = open(input_filename)
wf = open(output_filename, 'w')

vcf_line_no = 0

for line in f:
    vcf_line_no += 1
    if len(line) < 6:
        continue
    if line[:6] == '#CHROM':
        toks = re.split('\s+', line.rstrip())
        if len(toks) > 8:
            samples = toks[9]
        break
no_samples = len(samples)

for line in f:
    vcf_line_no += 1

    if line[0] == '#':
        continue

    toks = re.split('\s+', line.rstrip())

    if len(toks) < 8:
        continue

    [chrom, pos, uidbase, ref, alts, dummy, dummy, dummy] = toks[:8]
    reflen = len(ref)
    if uidbase == '.':
        uidbase = 'VAR' + str(vcf_line_no)
    if 'chr' not in chrom[:3]:
        chrom = 'chr' + chrom
    alts = alts.split(',')

    len_alts = len(alts)
    if len(toks) == 8:
        for altno in xrange(len_alts):
            alt = alts[altno]
            if len_alts == 1:
                uid = uidbase
            else:
                uid = uidbase + '_' + str(altno + 1)
            newpos, newref, newalt = extract_vcf_variant('+', pos, ref, alt)
            cravat_line = '\t'.join([uid, chrom, str(newpos), '+', newref, newalt, 'Unknown'])
            wf.write(cravat_line + '\n')
    elif len(toks) > 8:
        sample_datas = toks[9:]

        genotype_fields = {}
        genotype_field_no = 0
        for genotype_field in toks[8].split(':'):
            genotype_fields[genotype_field] = genotype_field_no

        if not ('GT' in genotype_fields):
            print 'No GT Field at line ' + str(vcf_line_no) + ' [' + line.strip() + ']'
            continue
        
        gt_field_no = genotype_fields['GT']

        for sample_no in xrange(len(sample_datas)):
            sample = samples[sample_no]
            sample_data = sample_datas[sample_no].split(':')
            gts = {}
            for gt in sample_data[gt_field_no].replace('/', '|').split('|'):
                if gt == '.':
                    continue
                else:
                    gts[int(gt)] = True
            for gt in gts.keys():
                if gt == 0:
                    continue
                else:
                    alt = alts[gt - 1]
                    if len_alts == 1:
                        uid = uidbase
                    else:
                        uid = uidbase + ':' + str(gt)
                    newpos, newref, newalt = extract_vcf_variant('+', pos, ref, alt)
                    cravat_line = '\t'.join([uid + '_' + sample, chrom, str(newpos), '+', newref, newalt, sample])
                    wf.write(cravat_line + '\n')
f.close()
wf.close()
