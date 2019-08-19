#!/usr/bin/env python
import optparse, os, shutil, subprocess, sys, tempfile, gzip


def __main__():
    descr = "sortmerna_wrapper.py"
    parser = optparse.OptionParser(description=descr)
    parser.add_option( '--ref_fasta', help='The reference fasta file' )
    parser.add_option( '--ref_index', help='The reference index files path' )
    parser.add_option( '--reads', help='FASTA/FASTQ read file. SortMeRNA accepts one Fastq file only. If running paired-end reads you will need to run the "Prepare paired-end reads" under the "NGS Assembly: Velvet" tool section.' )
    parser.add_option( '--fastx_aligned', help='Output FASTA/FASTQ file' )
    parser.add_option( '--fastx_rejected', help='Output FASTA/FASTQ file' )
    parser.add_option( '--sam', help='Output SAM alignment' )
    parser.add_option( '--blast', help='Output BLAST-like alignment' )
    parser.add_option( '--log', help='Output overall statistics' )
    parser.add_option( '--fast', action='store_true', help='Use --fast option' )
    parser.add_option( '--sensitive', action='store_true', help='Use --sensitive option' )
    parser.add_option( '--max_pos', help='Use --max_pos option' )
    parser.add_option( '--sq', help='Use -SQ option' )
    parser.add_option( '--feeling_lucky', help='Use --feeling-lucky option' )
    parser.add_option( '--num_alignments', help='Use --num_alignments option. Cannot be used together with --best option' )
    parser.add_option( '--best', help='Use --best option. Cannot be used together with --num_alignments option' )
    parser.add_option( '--paired_in', help='Use --paired_in option' )
    parser.add_option( '--paired_out', help='Use --paired_out option' )
    parser.add_option( '--match', help='Use --match option' )
    parser.add_option( '--mismatch', help='Use --mismatch option' )
    parser.add_option( '--gap_open', help='Use --gap_open option' )
    parser.add_option( '--gap_ext', help='Use --gap_ext option' )
    parser.add_option( '-N', help='Use -N option' )
    parser.add_option( '-F', action='store_true', help='Use -F option' )
    parser.add_option( '-R', action='store_true', help='Use -R option' )
    parser.add_option( '-a', help='Use -a option' )
    parser.add_option( '-e', help='Use -e option' )
    parser.add_option( '--passes', help='Use --passes option' )
    parser.add_option( '--edges', help='Use --edges option' )
    parser.add_option( '--num_seeds', help='Use --num_seeds option' )
    parser.add_option( '--full_search', help='Use --full_search option' )
    parser.add_option( '--fileSource', help='Whether to use a previously indexed reference sequence or one form history (indexed or history)' )
    parser.add_option( '--dbkey', help='Dbkey for reference genome' )
    (options, args) = parser.parse_args()
    if len(args) > 0:
        parser.error('Wrong number of arguments')

    # make temp directory for placement of indices
    tmp_index_dir = tempfile.mkdtemp()
    tmp_dir = tempfile.mkdtemp()
    # modify reads file to be in proper format
    (prefix, sep, suffix) = options.reads.rpartition('.')
    if suffix == 'gz':
        read_file_name = "%s/my_reads.fastq.gz" % ( tmp_index_dir )
        new_read_file_name = "%s/my_reads.fastq" % ( tmp_index_dir )

        # if file is gz, the file needs to be copied to a temp directory and unzipped. Linking won't work
        shutil.copy( options.reads, read_file_name )
        print "link is at: %s" % read_file_name
    else:
        new_read_file_name = "%s/my_reads.fastq" % ( tmp_index_dir )
        os.symlink( options.reads, new_read_file_name )
    if suffix == 'gz':
        #unzip_cmd = "inF = gzip.GzipFile("access_log.1.gz", 'rb'); s = inF.read(); inF.close(); outF = file("access_log.1", 'wb'); outF.write(s); outF.close()"  
        #reads_update_cmd = "inF = gzip.GzipFile(%s, 'rb'); s = inF.read(); inF.close(); outF = file(%s, 'wb'); outF.write(s); outF.close()" % ( read_file_name, new_read_file_name)
        read_file_unzip_cmd = "gunzip %s" % ( read_file_name )
        print "unzip cmd: %s" % read_file_unzip_cmd
        subprocess.Popen( args=read_file_unzip_cmd, shell=True ) 
        #os.symlink( options.reads, new_read_file_name )
    # index if necessary
    if options.fileSource == 'history':
        ref_file = tempfile.NamedTemporaryFile( dir=tmp_index_dir )
        ref_file_name = ref_file.name
        ref_file.close()
        os.symlink( options.ref_fasta, ref_file_name )
        #cmd1 = 'indexdb_rna --ref %s,%s/index_file' % ( options.ref, ref_file_name )
        #cmd1 = 'indexdb_rna --ref %s,%s/index_file %s' % ( options.ref, tmp_index_dir, options.fast )
        cmd1 = "indexdb_rna --ref %s,%s/index_file" % ( options.ref_fasta, tmp_index_dir )
        if options.fast is not None:
            cmd1 += ' --fast'
        if options.sensitive is not None:
            cmd1 += ' --sensitive'
        if options.max_pos is not None:
            cmd1 += ' --max_pos %s' % options.max_pos 
        print cmd1
        try:
            tmp = tempfile.NamedTemporaryFile( dir=tmp_index_dir ).name
            tmp_stderr = open( tmp, 'wb' )
            #proc = subprocess.Popen( args=cmd1, shell=True, cwd=tmp_index_dir, stderr=tmp_stderr.fileno() )
            proc = subprocess.Popen( args=cmd1, shell=True, stderr=tmp_stderr.fileno() )
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
            # clean up temp dirs
            if os.path.exists( tmp_index_dir ):
                shutil.rmtree( tmp_index_dir )
            if os.path.exists( tmp_dir ):
                shutil.rmtree( tmp_dir )
            raise Exception, 'Error indexing reference sequence. ' + str( e )
    else:
        ref_file_name = options.ref_fasta
        ref_index_path = options.ref_index
        #ref_file_name = os.path.basename(ref_file_name).rsplit( ".", 1 )[ 0 ]
    
    #cmd = "sortmerna --ref %s,%s/index_file --reads %s --log --blast --fastx --sam --aligned aligned --other rejected" % ( options.ref, ref_file_name, options.reads ) 
    #cmd = "sortmerna --ref %s,%s/index_file --reads %s %s %s --log --blast --fastx --sam --aligned aligned --other rejected" % ( options.ref, tmp_index_dir, options.reads, options.sq, options.feeling_lucky )
    if options.fileSource == 'history': 
        cmd = "sortmerna --ref %s,%s/index_file --reads %s %s" % ( options.ref_fasta, tmp_index_dir, new_read_file_name, options.sq )
    else:
        cmd = "sortmerna --ref %s,%s --reads %s %s" % ( options.ref_fasta, options.ref_index, new_read_file_name, options.sq ) 
    if options.num_alignments is not None:
        cmd += " --num_alignments %s" % options.num_alignments
    if options.feeling_lucky is not None:
        cmd += " %s" % options.feeling_lucky
    if options.best is not None:
        cmd += " --best %s" % options.best
    if options.paired_in is not None:
        cmd += " %s" % options.paired_in
    if options.paired_out is not None:
        cmd += " %s" % options.paired_out
    if options.match is not None:
        cmd += " --match %s" % options.match
    if options.mismatch is not None:
        cmd += " --mismatch %s" % options.mismatch
    if options.gap_open is not None:
        cmd += " --gap_open %s" % options.gap_open
    if options.gap_ext is not None:
        cmd += " --gap_ext %s" % options.gap_ext
    if options.N is not None:
        cmd += " -N %s" % options.N
    if options.F is not None:
        cmd += ' -F'
    if options.R is not None:
        cmd += ' -R'
    if options.a is not None:
        cmd += " -a %s" % options.a
    if options.e is not None:
        cmd += " -e %s" % options.e
    if options.passes is not None:
        cmd += " --passes %s" % options.passes
    if options.edges is not None:
        cmd += " --edges %s" % options.edges
    if options.num_seeds is not None:
        cmd += " --num_seeds %s" % options.num_seeds
    if options.full_search is not None:
        cmd += " %s" % options.full_search
    cmd += " --log --blast --fastx --sam --aligned aligned --other rejected"
    # perform alignments
    buffsize = 1048576
    try:
        # need to nest try-except in try-finally to handle 2.4
        try:
            #tmp = tempfile.NamedTemporaryFile( dir=tmp_dir ).name
            tmp = tempfile.NamedTemporaryFile( dir=tmp_index_dir ).name
            tmp_stderr = open( tmp, 'wb' )
            print "The cmd is %s" % cmd
            #proc = subprocess.Popen( args=cmd, shell=True, cwd=tmp_dir, stderr=tmp_stderr.fileno() )
            #proc = subprocess.Popen( args=cmd, shell=True, cwd=tmp_index_dir, stderr=tmp_stderr.fileno() )
            proc = subprocess.Popen( args=cmd, shell=True, stderr=tmp_stderr.fileno() )
            returncode = proc.wait()
            tmp_stderr.close()
            # get stderr, allowing for case where it's very large
            tmp_stderr = open( tmp, 'rb' )
            stderr = ''
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
            raise Exception, 'Error generating alignments. ' + str( e )
        
    finally:
        # clean up temp dir
        if os.path.exists( tmp_index_dir ):
            shutil.rmtree( tmp_index_dir )
        #if os.path.exists( tmp_dir ):
        #    shutil.rmtree( tmp_dir )

if __name__ == "__main__":
    __main__()


