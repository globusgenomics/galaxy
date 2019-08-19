#!/usr/bin/env python

import sys, os

def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit()
    
def main():
    outfile = sys.argv[1]
    infile = sys.argv[2]

    if len(sys.argv) >3:
        columnList = sys.argv[3]
    
    try:
        fout = open(sys.argv[1],'w')
    except:
        stop_err("Output file cannot be opened for writing.")
        
    try:
        fin = open(sys.argv[2],'r')
    except:
        stop_err("Input file cannot be opened for reading.")

    if len(sys.argv) >3:
        cList=columnList.split(',')
        selColumn=[ "-k %s,%s" % (cList[i], cList[i]) for i in range(0,len(cList))] 
        par=(' ').join(selColumn)
        cmdline = "sort -u %s %s  > %s" % (par,infile, outfile) 
    else:
        cmdline = "sort %s | uniq > %s" % (infile, outfile)

    try:
        os.system(cmdline)
    except:
        stop_err("Error encountered with cat.")
        
if __name__ == "__main__": main()
