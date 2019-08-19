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
        input_dir=opts.input+'/output/'
        files=glob.glob(input_dir + "*/table/*.LargeDeletion.txt")
    else:
        sys.exit()
    
    with open("%s/allyears.csv" % opts.outputdir, "w") as output:
        for fname in files:	
            samplename=os.path.basename(fname).split('.')[0]  # sample name (last column)
            with open(fname, 'r') as f:
                next(f)
                for line in f:
                    chrname=line.split('\t')[0]; start=line.split('\t')[4]; stop=line.split('\t')[5]
                    newline='\t'.join([chrname, start, stop, 'DEL', samplename])
                    if "NA" not in newline:
                        output.write(newline); output.write("\n")
    shutil.copy("%s/allyears.csv" % opts.outputdir, opts.output)
        
if __name__ == '__main__':
    main()
