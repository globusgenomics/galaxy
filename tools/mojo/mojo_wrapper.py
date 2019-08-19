#!/usr/bin/env python


import optparse, os, sys, subprocess, tempfile, shutil
#from galaxy import eggs
#import pkg_resources; pkg_resources.require( "bx-python" )
#from bx.cookbook import doc_optparse
#from galaxy import util

def stop_err( msg ):
    sys.stderr.write( '%s\n' % msg )
    sys.exit()

def __main__():
    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '-1', '--fq1', dest='fq1', action='append', type="string", help='The input forward fastqs' )
    parser.add_option( '-2', '--fq2', dest='fq2', action='append', type="string", help='The input reverse fastqs' )
    parser.add_option( '', '--output-fusions', dest='output_fusions', help='The output fusions calls' )
    parser.add_option( '', '--output-pileup', dest='output_pileup', help='The output fusions pileup' )
    parser.add_option( '-R', '--ref-path', dest='ref_path', help='The path to the reference files for the genome' )
    parser.add_option( '-p', '--pass_through', dest='pass_through_options', action='append', type="string", help='These options are passed through directly to GATK, without any modification.' )
    parser.add_option( '', '--output-tmp', dest='output_tmp', help='The output tmp directory' ) 
    ( options, args ) = parser.parse_args()

    if options.output_tmp is not None:
        tmp_dir = options.output_tmp
        os.mkdir(tmp_dir)
    else:    
        tmp_dir = tempfile.mkdtemp()

    # modify the config file
    path_to_mojo_bin = os.path.dirname(os.popen("which %s" % "MOJO").read().strip())
    path_to_mojo = "%s/.." % (path_to_mojo_bin)
    path_to_mojo_tools = "%s/../external" % (path_to_mojo_bin)
    path_to_mojo_ref = options.ref_path
    config_file = "%s/Sample.configfile.txt" % path_to_mojo
    run_config_file = "%s/configfile.txt" % tmp_dir

    fw = open(run_config_file, "w")
    fw.write("mojo_install_dir = %s/\n" % path_to_mojo_bin)
    fw.write("mojo_reference_dir = %s/\n" % path_to_mojo_ref)
    fw.write("mojo_tools_dir = %s/\n" % path_to_mojo_tools)
    fw.write("samtools_binary =   %s/samtools/samtools\n" % path_to_mojo_tools)
    fw.write("bwa_binary =   %s/bwa/bwa\n" % path_to_mojo_tools)
    fw.write("bowtie2_binary =  %s/bowtie2/bowtie2\n" % path_to_mojo_tools)
    fw.write("bowtie2_build_binary = %s/bowtie2/bowtie2-build\n" % path_to_mojo_tools)
    fw.write("blat_binary = %s/blat/blat\n" % path_to_mojo_tools)
    fw.write("bwa_transcriptome_index        =   %s/transcriptome/transcriptome\n" % path_to_mojo_ref)
    fw.write("bowtie2_all_isoforms_index     =   %s/all_isoforms/all.isoforms\n" % path_to_mojo_ref)
    fw.write("bowtie2_genome_index           =   %s/genome/genome\n" % path_to_mojo_ref)
    fw.write("blat_genome_2bit               =   %s/blat_filter_refs/genome.2bit\n" % path_to_mojo_ref)
    fw.write("blat_reference_junctions_2bit  =   %s/blat_filter_refs/allReferenceJunctsDB.2bit\n" % path_to_mojo_ref)
    fw.write("blat_filter_chroms_dir         =   %s/blat_filter_refs/\n" % path_to_mojo_ref)
    fw.write("master_gene_file               =   %s/gene_model/Gene.txt\n" % path_to_mojo_ref)
    fw.write("master_exon_file               =   %s/gene_model/Exon.txt\n" % path_to_mojo_ref)
    fw.write("master_isoform_file            =   %s/gene_model/Isoform.txt\n" % path_to_mojo_ref)
    fw.write("megablast_output_file          =   %s/gene_model/gene2gene.megablast.txt\n" % path_to_mojo_ref)
    fw.write("repeat_masker_file             =   %s/gene_model/rmsk.regions.txt\n" % path_to_mojo_ref)
    fw.write("max_bwa_mem                    = 6\n")
    fw.write("min_span                       = 2,2,80000000\n")
    fw.write("read_through                   = 200000\n")
    fw.write("junct_mismatch                 = 0.03\n")
    fw.write("split_fastq_binary             = %s/StreamNthFastqSplit\n" % path_to_mojo_bin)
    fw.write("filter_junct_output_binary     = %s/FilterJunctAlignOutput\n" % path_to_mojo_bin)
    fw.close()

    pass_through = None
    if options.pass_through_options:
        pass_through = ' '.join( options.pass_through_options )
    else:
        pass_through = ""

    try:
        command = 'MOJO %s --config %s --cores 32 --sample_name mojo_run --output_dir %s --fq1 %s --fq2 %s' % (pass_through, run_config_file, tmp_dir, ",".join(options.fq1), ",".join(options.fq2) )

        print command
        tmpdir = tempfile.mkdtemp()
        tmp = tempfile.NamedTemporaryFile( dir=tmp_dir ).name
        tmp_stderr = open( tmp, 'wb' )
        proc = subprocess.Popen( args=command, shell=True, cwd=tmpdir, stderr=tmp_stderr.fileno() )
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

        # move the output files to its final destination
        fusion_calls = "%s/mojo_run/mojo_run.fusions" % (tmp_dir)
        fusion_pileup = "%s/mojo_run/mojo_run.fusions.pileup" % (tmp_dir)
        shutil.move(fusion_calls, options.output_fusions)
        shutil.move(fusion_pileup, options.output_pileup)
    except Exception, e:
        raise Exception, 'Error generating alignments. ' + str( e ) 
    

if __name__=="__main__": __main__()

