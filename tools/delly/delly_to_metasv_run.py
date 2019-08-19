#!/usr/bin/python

import sys, re, glob
from optparse import OptionParser
#import contra_to_metasv_reformat as con2metasv

def header(line):
    """
    add '>' if the symbol is missing at the end of line
    """
    line=line.strip()
    if re.match("##.*<.*", line) and line[-1:] != '>':
        line=line+'>'
    return line

def body(line):
    """
    modify the INFO column by adding SNVTYPE and SNVLEN
    if the last field each line is gain/loss, add SNVTYPE as DUP/DEL
    calcuate SNVLEN by POS - END
    """
    line=line.strip()
    m=line.split('\t')
    EndPosition=re.match(r".*;END=(\d+).*",m[7])
    SVLen=int(EndPosition.group(1))-int(m[1])
    m[7] = 'SVLEN='+str(SVLen)+';'+EndPosition.group(0)
    return '\t'.join(m)

def main():
    """
       usage: ./contra_vcf_reformat_metasv.py [file_path]
       This will get all the vcf files under the directory
       and generate reformated files under the current directory
    """
    parser = OptionParser()
    parser.add_option('--fin', dest='fin', nargs=1, default=None, 
                      action='store', type='string', help="vcf input file name.")
    parser.add_option('--fout', dest='fout', nargs=1, action='store', 
                       type='string', default=None, help="modified vcf output file name.")
    (options, args) = parser.parse_args()
    try:
        if options.fin:
            f1 = open(options.fin, 'r')
            f2 = open(options.fout, 'w')
            for line in f1:
                if line.lstrip().startswith('##'):
                    f2.write(header(line)+"\n")
                elif line.lstrip().startswith('#'):
                    f2.write(line)
                else:
                    f2.write(body(line)+"\n")
    except:
        sys.stdout.write( 'Could not find input file\n' )        
        
if __name__ == '__main__':
    main()
