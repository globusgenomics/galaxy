#!/usr/bin/python

import sys, re, glob, os, shutil

def main():
    """
    """
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('--input', dest='input', type='string', help="vcf input file directory path.")
    parser.add_option('--output', dest='output', type='string', help="output file name.")
    parser.add_option('--output-dir', dest='outputdir', type='string', help="output directore.")
    (opts, args) = parser.parse_args()

    if not os.path.isdir(opts.outputdir):
        os.mkdir(opts.outputdir)

    if opts.input:
        files=glob.glob(opts.input + "/*.vcf")
    else:
        sys.exit()
    
    with open("%s/all.csv" % opts.outputdir, "w") as output:
        for fname in files:
            samplename=os.path.basename(fname).split('.')[0]  # sample name (last column)
            with open(fname, 'r') as f:
                for line in f:
                    line=line.strip()
                    if re.match(r'^chr', line):
                        match=re.match(r'(\w+)\t(\d+)\t..*\t<(\w+)>.*END=(\d+);.*SVLEN.*', line)
                        newline='\t'.join([match.group(1), match.group(2), match.group(4), match.group(3), samplename])
                        output.write(newline); output.write("\n")
    shutil.copy("%s/all.csv" % opts.outputdir, opts.output)
        
if __name__ == '__main__':
    main()
