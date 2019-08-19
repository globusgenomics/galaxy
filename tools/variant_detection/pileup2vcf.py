#!/usr/bin/env python
"""
Converts pileup format to VCF format.
usage: 
perl /nfs/software/galaxy/tools/variant_detection/pileup2vcf.pl -r reference_genome < input > outputvcf
"""

import optparse, os, sys, subprocess, tempfile, shutil
from galaxy import eggs
from galaxy import util

galaxyhome=os.environ.get('GALAXY_HOME')

def stop_err( msg ):
    sys.stderr.write( '%s\n' % msg )
    sys.exit()

def check_seq_file( dbkey, cached_seqs_pointer_file ):
    seq_path = ''
    for line in open( cached_seqs_pointer_file ):
        line = line.rstrip( '\r\n' )
        if line and not line.startswith( '#' ) and line.startswith( 'index' ):
            fields = line.split( '\t' )
            if len( fields ) < 3:
                continue
            if fields[1] == dbkey:
                seq_path = fields[2].strip()
                break
    return seq_path

def __main__():
    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '', '--input1', dest='input1', help='The input pileup file' )
    parser.add_option( '', '--dbkey', dest='dbkey', help='The build of the reference dataset' )
    parser.add_option( '', '--ref_file', dest='ref_file', help='The reference dataset from the history' )
    parser.add_option( '', '--output1', dest='output1', help='The output VCF file' )
    parser.add_option( '', '--index_dir', dest='index_dir', help='GALAXY_DATA_INDEX_DIR' )
    ( options, args ) = parser.parse_args()

    cached_seqs_pointer_file = '%s/sam_fa_indices.loc' % options.index_dir
    if not os.path.exists( cached_seqs_pointer_file ):
        stop_err( 'The required file (%s) does not exist.' % cached_seqs_pointer_file )

    seq_path = check_seq_file( options.dbkey, cached_seqs_pointer_file )

    cmd = "perl " + galaxyhome + "/tools/variant_detection/pileup2vcf.pl -r " + seq_path + " < " + options.input1 + " > " + options.output1

    try:
	tmp_dir = tempfile.mkdtemp()
        tmp = tempfile.NamedTemporaryFile( dir=tmp_dir ).name
        tmp_stderr = open( tmp, 'wb' )
	proc = subprocess.Popen( args=cmd, shell=True, cwd=tmp_dir, stderr=tmp_stderr.fileno() )
        returncode = proc.wait()
        tmp_stderr.close()
        # get stderr, allowing for case where it's very large
        tmp_stderr = open( tmp, 'rb' )
        stderr = ''
        buffsize = 1048576
        try:
            while True:
                stderr += tmp_stderr.read( buffsize )
                if not stderr or len( stderr ) % buffsize != 0:
                    break
        except OverflowError:
            pass
        tmp_stderr.close()
        if returncode != 0:
            raise Exception, stderr 
    except Exception, e:
        #clean up temp files
        if os.path.exists( tmp_dir ):
            shutil.rmtree( tmp_dir )
        stop_err( 'Error with reference (%s), %s' % ( options.ref_file, str( e ) ) )

if __name__=="__main__": __main__()
