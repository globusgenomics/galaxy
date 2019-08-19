#!/usr/bin/python

CHUNK_SIZE = 2**20 #1mb

import argparse, os, shutil, subprocess, sys, tempfile, shlex, vcf, pysam

parser = argparse.ArgumentParser(description='')
parser.add_argument ( '--bed', dest='bed', help='the bed file')
parser.add_argument ( '--sam', dest='sam', help='the sam file')
parser.add_argument ( '--bam', dest='bam', help='the bam file')	
#parser.add_argument ( '-p', nargs = '*', dest='pass_through_options', help='These options are passed through directly to contra, without any modification')
parser.add_argument ( '-o', dest='output', help='output log file' )
parser.add_argument ( '--output-dir', dest='output_dir', help='output dir' )

def execute( cmd, output="" ):
    tmp_dir = tempfile.mkdtemp()
    try:
            err = open(tmp_dir+"/errorLog", 'a')
            if output != "":
                    out = open(output, 'w')
            else:
                    out = subprocess.PIPE
            process = subprocess.Popen( args=shlex.split(cmd), stdout=out, stderr=err )
            process.wait()
            err.close()
            if out != subprocess.PIPE:
                    out.close()
    except Exception, e:
            sys.stderr.write("problem doing : %s\n" %(cmd))
            sys.stderr.write( '%s\n\n' % str(e) )

def __main__():
    args = parser.parse_args()

    if args.output_dir is not None:
        output_dir = args.output_dir
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
    else:
        output_dir = tempfile.mkdtemp();

#    if args.pass_through_options is not None:
#        ptc = ' '.join(args.pass_through_options )
        ## In our implementation we are assuming BAI files are already generated. 
        ## Will need to modify for future implementations
    input_dir = "%s/input" % output_dir
    if os.path.exists(input_dir):
        shutil.rmtree(input_dir)
    os.mkdir(input_dir)

    try:
        if args.bam is not None:
            bamfiles = args.bam.split()
            linked_bamfiles = []
            for bamfile in bamfiles:
                name = os.path.basename(bamfile)
                linkname = "%s/%s" % (input_dir, name)
                os.symlink(bamfile, linkname)
                linked_bamfiles.append(linkname)
                execute("samtools index " + linkname)

            inputfiles = " ".join(linked_bamfiles)
        elif args.sam is not None:
            samfiles = args.sam.split()
            linked_samfiles = []
            for samfile in samfiles:
                name = os.path.basename(samfile)
                linkname = "%s/%s" % (input_dir, name)
                os.symlink(samfile, linkname)
                linked_samfiles.append(linkname)

            inputfiles = " ".join(linked_samfiles)
        else: 
            bedfiles = args.bed.split()
            linked_bedfiles = []
            for bedfile in bedfiles:
                name = os.path.basename(bedfile)
                linkname = "%s/%s" % (input_dir, name)
                os.symlink(bedfile, linkname)
                linked_bedfiles.append(linkname)

            inputfiles = " ".join(linked_bedfiles)
    except Exception, e:
        print("problem with input files " + str(e))

#    if args.pass_through_options is not None:
#        cmd="makeTagDirectory %s %s %s" %(output_dir, ptc, inputfiles) 
#    else:
    cmd="makeTagDirectory %s %s" %(output_dir, inputfiles)

    print("cmd: %s" % cmd)
    stdout = tempfile.NamedTemporaryFile( prefix="maketagdir-stdout-", dir=output_dir )
    stderr = tempfile.NamedTemporaryFile( prefix="maketagdir-stderr-", dir=output_dir )  

    proc = subprocess.Popen( args=cmd, stdout=stdout, stderr=stderr, shell=True)
    return_code = proc.wait()

    if return_code:
        stderr_target = sys.stderr
    else:
        stderr_target = sys.stdout
    stderr.flush()
    stderr.seek(0)

    while True:
        chunk = stderr.read( CHUNK_SIZE )
        if chunk:
            stderr_target.write( chunk )
        else:
            break
    ##shutil.copy("%s/
    stderr.close()
    stdout.close() 

    shutil.copy('%s/tagInfo.txt' % output_dir, args.output)

if __name__=="__main__":
	__main__()
