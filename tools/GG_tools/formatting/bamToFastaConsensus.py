#!/usr/bin/env python
"""
Converts BAM data to consensus FASTA files.
usage: bamToFastaConsensus.py [options]
   -b BAM file
   -r Fasta Reference
   -o output File name
   --no-split: do not split the output file into chromosomes
   --files_path output for chr consensus files
"""

import optparse, os, sys, subprocess, tempfile, shutil

galhtmlprefix = """<?xml version="1.0" encoding="utf-8" ?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<meta name="generator" content="Galaxy tool output - see http://getgalaxy.org/" />
<title></title>
<link rel="stylesheet" href="/static/style/base.css" type="text/css" />
</head>
<body>
<div class="document">
"""
galhtmlpostfix = """</div></body></html>\n"""

def stop_err( msg ):
    sys.stderr.write( '%s\n' % msg )
    sys.exit()

def cleanup_before_exit( tmp_dir ):
    if tmp_dir and os.path.exists( tmp_dir ):
        shutil.rmtree( tmp_dir )

def create_dict_idx_from_fasta( fasta_filename, target_dir = None ):
    if target_dir is None:
        target_dir = os.getcwd()
    basename = os.path.basename(fasta_filename)
    link_filename = os.path.join( target_dir, "%s.fasta" % ( basename ) )
    os.symlink( fasta_filename, link_filename )

    # if input reference is a fasta file, create the dict and index file needed 
    stderr_picard = tempfile.NamedTemporaryFile( prefix = "picard_stderr" ).name
    #print stderr_picard
    dict_filename = os.path.join( target_dir, "%s.dict" % ( basename ) )
    picard_jar = os.path.join( os.environ['PICARD_JAR'], "CreateSequenceDictionary.jar")

    cmd1 = "java -jar %s R= %s O= %s" % (picard_jar, link_filename, dict_filename)
    cmd2 = "samtools faidx %s" % (link_filename)
    process = subprocess.Popen(cmd1, shell=True, stderr=open( stderr_picard, 'wb' ), cwd=target_dir)
    rval = process.wait()
    if rval:
        for line in open( stderr_picard ):
            print >> sys.stderr, line
            os.unlink( stderr_picard ) #clean up
            cleanup_before_exit( tmp_dir )
            raise Exception( "Error creating reference dict file" )
        os.unlink( stderr_picard ) #clean up

    stderr_sam = tempfile.NamedTemporaryFile( prefix = "sam_stderr" ).name
    process2 = subprocess.Popen(cmd2, shell=True, stderr=open( stderr_sam, 'wb' ), cwd=target_dir)
    rval2 = process2.wait()
    if rval2:
        for line in open( stderr_sam ):
            print >> sys.stderr, line
            os.unlink( stderr_sam ) #clean up
            cleanup_before_exit( tmp_dir )
            raise Exception( "Error creating reference index file" )
        os.unlink( stderr_picard ) #clean up        

def create_bam_idx(bamFile, target_dir = None):
    if target_dir is None:
        target_dir = os.getcwd()

    index_name = "%s/%s.bai" % (target_dir, os.path.basename(bamFile))
    stderr_name = tempfile.NamedTemporaryFile( prefix = "bam_stderr" ).name
    command = 'samtools index %s %s' % ( bamFile, index_name )
    print command
    proc = subprocess.Popen( args=command, shell=True, stderr=open( stderr_name, 'wb' ) )
    exit_code = proc.wait()
    if exit_code:
        for line in open( stderr_name ):
            print >> sys.stderr, line
        os.unlink( stderr_name ) #clean up
        #cleanup_before_exit( tmp_dir )
        raise Exception( "Error indexing BAM file" )
    os.unlink( stderr_name ) #clean up


def __main__():
    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '-b', '--bam-input', dest='bamFile', help='The input BAM dataset' )
    parser.add_option( '-i', '--indexBam', dest='bamFileIndex', help='The BAM index file' )
    parser.add_option( '-r', '--ref-input', dest='refFile', help='The reference FASTA dataset' )
    parser.add_option( '-o', '--output', dest='outputFile', help='Output file will either be a log file pointing to the location of each of the files or the combined FASTA Consensus file' )
    parser.add_option( '', '--no-split', dest='splitChrs', action="store_true", help="Combine all chromosome consensus FASTA files into one" ) 
    parser.add_option( '', '--files_path', dest='extra_files_path', help="Location for chr consensus output" )
    ( options, args ) = parser.parse_args()

    # exit if input file empty
    if os.path.getsize( options.bamFile ) == 0:
        raise Exception, 'Initial BAM file empty'

    #print "Creating directory: " + args.files_path
    if not os.path.exists(options.extra_files_path):
        os.makedirs(options.extra_files_path)

    tmp_dir = tempfile.mkdtemp(prefix="bamIndex", dir=options.extra_files_path)

    # create link of input BAM file to the output directory
    basename = os.path.basename(options.bamFile)
    link_filename = os.path.join( tmp_dir, "%s.bam" % ( basename ) )
    os.symlink( options.bamFile, link_filename )

    # create a BAM index file
    if options.bamFileIndex:
        # create a link to the index file in the same directory as the bam file
        bai_link_filename = os.path.join(tmp_dir, "%s.bam.bai" % ( basename ) )
        os.symlink( options.bamFileIndex, bai_link_filename)
    else:
        create_bam_idx(link_filename, target_dir = tmp_dir)

    # does dict and fasta index files exist? if not, create link to a temp directory and create index, dict files
    fasta_path_basename = ('.').join(options.refFile.split('.')[:-1])
    base = os.path.basename(options.refFile)
    if not os.path.isfile("%s.dict" % fasta_path_basename):
        create_dict_idx_from_fasta( options.refFile, target_dir = tmp_dir)
        index_dir = tmp_dir
        fastapath = "%s/%s.fasta" % (index_dir, base)
    else:
        index_dir = os.path.realpath(options.refFile)
        fastapath = options.refFile

    #print "Creating directory: " + args.files_path
    if not os.path.exists(options.extra_files_path):
        os.makedirs(options.extra_files_path)

    ## run the consensus fasta tool
    #samtools view -H $inputBam | grep "\@SQ" | sed 's/^.*SN://g' | cut -f 1 | grep -v GL | grep -v NC_ |xargs -I {} -n 1 -P 32 sh -c "samtools mpileup -BQ0 -d 100000 -uf $inputref -r {} $inputBam | bcftools view -vcg - | vcfutils.pl vcf2fq - | seqtk fq2fa - 20 > tmp.{}.fasta"
    cmd = "samtools view -H %s | grep \"\\@SQ\" | sed \'s/^.*SN://g\' | cut -f 1 | grep -v GL | grep -v NC_  |xargs -I {} -n 1 -P 30 sh -c \"samtools mpileup -BQ0 -d 100000 -uf %s -r {} %s | bcftools view -vcg - | vcfutils.pl vcf2fq - | seqtk fq2fa - 20 > %s/consensus.{}.fasta\"" % (link_filename, fastapath, link_filename, options.extra_files_path)
    print cmd

    stderr_name = tempfile.NamedTemporaryFile( prefix = "consensus_stderr" ).name
    proc = subprocess.Popen( args=cmd, shell=True, stderr=open( stderr_name, 'wb' ) )
    exit_code = proc.wait()
    if exit_code:
        for line in open( stderr_name ):
            print >> sys.stderr, line
            os.unlink( stderr_name ) #clean up
            raise Exception( "Error creating consensus sequences " )
        os.unlink( stderr_name ) #clean up

    if options.splitChrs:
        # concatenate all chr consensus files
        cmd2 = "cat %s/consensus.*.fasta > %s" % (options.extra_files_path, options.outputFile)
        print cmd2

        stderr_name = tempfile.NamedTemporaryFile( prefix = "cat_stderr" ).name
        proc = subprocess.Popen( args=cmd2, shell=True, stderr=open( stderr_name, 'wb' ) )
        exit_code = proc.wait()
        if exit_code:
            for line in open( stderr_name ):
                print >> sys.stderr, line
                os.unlink( stderr_name ) #clean up
                raise Exception( "Error concatenating consensus sequences " )
            os.unlink( stderr_name ) #clean up


    else:
        # create an output file with the name of all files created and stored in the library, and name of library
        #Create HTML file
        f = open(options.outputFile,'w')
        f.write(galhtmlprefix)
        #f.write("<p>Consensus files generated from BAM file: %s</p>\n" % (os.path.basename(options.bamFile))
        for fastaF in os.listdir("%s" % (options.extra_files_path)):
            f.write('<a href="%s">%s</a><br>\n' % (fastaF, fastaF))
        f.write(galhtmlpostfix)
        f.close()




if __name__=="__main__": __main__()
