#!/usr/bin/python

import tempfile, os, sys, optparse, subprocess

def main():

    # Parse the command line options
    usage = "Usage: vcfconvert.py [options]"
    parser = optparse.OptionParser(usage = usage)
    parser.add_option("-i", "--in",
                      action="store", type="string",
                      dest="vcfFile", help="input vcf file (stdin for piped vcf)")
    parser.add_option("-o", "--out",
                      action="store", type="string",
                      dest="output", help="output validation file")
    (options, args) = parser.parse_args()

    # Check that a vcf file is given.
    if options.vcfFile == None:
        parser.print_help()
        print >> sys.stderr, "\nInput vcf file (--in, -i) is required."
        exit(1)

    # Check that the version of the VCF file is VCFv4.2, else print the same file as output
    fh = open(options.vcfFile, "r")
    cmd = ""
    for line in fh:
        if "##fileformat=VCFv4.2" in line:
            cmd = "cat %s | sed 's/VCFv4.2/VCFv4.1/' | sed 's/,Version=3>/>/' | sed 's/,Version=\"3\">/>/' | sed 's/Number=R/Number=./' > %s " % (options.vcfFile, options.output)
            break
        elif "##fileformat=VCFv4.1" in line or "##fileformat=VCFv4.0" in line:
            cmd = "cat %s > %s" % (options.vcfFile, options.output)
        break
    fh.close()
    stderr_reformat = tempfile.NamedTemporaryFile( prefix = "reformat_stderr" ).name
    process = subprocess.Popen(cmd, shell=True, stderr=open( stderr_reformat, 'wb' ))
    rval = process.wait()

if __name__ == "__main__":
  main()

