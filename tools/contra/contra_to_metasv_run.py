#!/usr/bin/python

import sys, re, glob
#import contra_to_metasv_reformat as con2metasv

def header(line):
    """
    add '>' if the symbol is missing at the end of line
    """
    line=line.strip()
    if re.match("##.*<.*", line) and line[-1:] != '>':
        line=line+'>'
    return line

vcf_data={}
def body(line):
    """
    modify the INFO column by adding SNVTYPE and SNVLEN
    if the last field each line is gain/loss, add SNVTYPE as DUP/DEL
    calcuate SNVLEN by POS - END
    """
    line=line.strip()
    m=line.split('\t')
    EndPosition=re.match(r"SVTYPE=CNV;END=(\d+)",m[7])
    SVType=re.match(r"(...):(\w+)",m[9])

    if SVType.group(2) == 'gain':
        vcf_data['SVTYPE']='DUP'
        SVLen=int(EndPosition.group(1))-int(m[1])
    elif SVType.group(2) == 'loss':
        vcf_data['SVTYPE']='DEL'
        SVLen=int(m[1])-int(EndPosition.group(1))
    else:
        vcf_data['SVTYPE']='None'
    vcf_data['SVLEN']=SVLen
    m[7]=';'.join(['%s=%s' % (key, value) for (key, value) in vcf_data.items()])
    m[7]=m[7]+';END='+EndPosition.group(1)
    return '\t'.join(m)

def main():
    """
       usage: ./contra_vcf_reformat_metasv.py [file_path]
       This will get all the vcf files under the directory
       and generate reformated files under the current directory
    """
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('--dir', dest='dirPath', nargs=1, default=None,
                       help="Directory path for vcf files.")
    parser.add_option('--fin', dest='fin', nargs=1, default=None, 
                      action='store', type='string', help="vcf input file name.")
    parser.add_option('--fout', dest='fout', nargs=1, action='store', 
                       type='string', default=None, help="modified vcf output file name.")
    (options, args) = parser.parse_args()
    print options 
#    print options
    try:
        if options.dirPath:
            dir_path=options.dirPath
            flist=glob.glob(dir_path+"/*.vcf")

            for i in xrange(0, len(flist)):
                fname=re.search('.*/(.*).vcf',flist[i]).group(1)
                fout=options.fout+'.vcf'
                f1 = open(flist[i], 'r')
                f2 = open(fout, 'w')
                for line in f1:
                    if line.lstrip().startswith('##'):
                        f2.write(header(line)+"\n")
                    elif line.lstrip().startswith('#'):
                        f2.write(line)
                    else:
                        f2.write(body(line)+"\n")
    except:
        sys.stdout.write( 'Could not determine directory path\n' )
    
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
