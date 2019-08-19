#!/usr/bin/python
import sys, optparse, os, tempfile, subprocess, shutil
CHUNK_SIZE = 2**20 #1mb

def __main__():

    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '-i', '--input', dest="input_file", action="store", type="string", help="Input VCF" )
    parser.add_option( '-o', '--output', dest="output_file", action="store", type="string", help="Output VCF" )
    parser.add_option( '--assembly', dest="assembly", action="store", type="string", help="Cached human assembly" )
    parser.add_option( '--options', dest="db_options", action="store", type="string", help="db options" )
    parser.add_option( '--vcf', dest="output_format", action="store_true", help="Format output in VCF" )
    (options, args) = parser.parse_args()

    fasta_path = None
    vep_index = "/mnt/galaxyIndices/genomes/vep"
    db_path = "%s/database/" % vep_index
    assembly = options.assembly
    fork_count = "32"
    version = "95"
    if assembly == "GRCh37":
        fasta_path = "%s/fasta/Homo_sapiens.GRCh37.dna.toplevel.fa" % vep_index
    elif assembly == "GRCh38":
        fasta_path = "%s/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa" % vep_index

    db_string = ""
    if "None" not in options.db_options:
        db_string = "--" + options.db_options.replace(",", " --")
    if "everything" in options.db_options:
        db_string = "--everything"
    if options.output_format is True:
        vcf_output_option = "--vcf "
    else:
        vcf_output_option = " "

    cmd = "vep --buffer_size 1000 --offline -i %s -o %s --cache --dir %s --force_overwrite --merged --cache_version %s --assembly %s --fasta %s --fork %s %s --everything %s" % (options.input_file, options.output_file, db_path, version, assembly, fasta_path, fork_count, db_string, vcf_output_option)
    print cmd
    stderr = tempfile.NamedTemporaryFile( prefix = "cmd_stderr" )
    proc = subprocess.Popen( args=cmd, stderr=stderr, shell=True )

    exit_code = proc.wait()

    if exit_code:
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
    stderr.close()

if __name__=="__main__": __main__()
