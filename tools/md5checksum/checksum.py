#!/usr/bin/python

import sys, os, re
import hashlib, optparse
'''
Usage:
Usage:
python checksum.py -d [file1] -d [file2] ... -m [md5checksum file] -f [output file]
last updated: March 24, 2016
'''

(HASH, FILE) = (0, 2)

def generate_file_md5(fin, blocksize=2**20):
'''
performs md5 checksum and returns a hashcode
'''

    m = hashlib.md5()
    with open(fin, "rb" ) as f:
        while True:
            buf = f.read(blocksize)
            if not buf:
                break
            m.update( buf )
    return m.hexdigest()

def main():
    check=[]; checksum={}
    parser = optparse.OptionParser()
    parser.add_option( '-d', '--dataset', dest='datasets', action='append', type="string", nargs=1, help='input file for checksum' )
    parser.add_option( '-m', '--md5checksum', dest='md5', help='If specified, the md5 will be used to check' )
    parser.add_option("-f", "--file", dest="filename", help="write report to FILE", metavar="FILE")

    (options, args) = parser.parse_args()

# parse md5 file
    if options.md5:
        with open(options.md5, 'r') as f1:
            for line in f1:
                hashList = line.strip().split(' ')[HASH]
                fileList = line.strip().split(' ')[FILE]
                checksum[hashList]=fileList

# run md5 checksum for each file
    if options.datasets:
        for f2 in options.datasets:
            check.append(generate_file_md5(f2))

# print the final result
    with open(options.filename,'w') as fout:
        for i in check:
            try:
                fout.write("%s %s pass\n" % (i, checksum[i]))
            except KeyError as e:
                fout.write("%s failed\n" % e)

if __name__ == '__main__':
    main()
