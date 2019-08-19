# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""
Runs the FermiKit Assembly program 
See below for options
"""

import time, optparse, os, shutil, subprocess, sys, tempfile, multiprocessing, gzip
CHUNK_SIZE = 2**20 #1mb

def time_stamp():
    return time.strftime('%d/%m/%Y %H:%M:%S\n', time.localtime(time.time()))

def get_ncores():
    return multiprocessing.cpu_count()

def get_linuxRAM():
    # return the available memory on the instance in MB
    totalMemory = os.popen("free -m").readlines()[1].split()[1]
    return int(totalMemory)

def get_readLength(input_fastq):
    if file_type(input_fastq) == "gz":
        fh_zip = gzip.open(input_fastq)
        fh_zip.readline()
        return len(fh_zip.readline().rstrip()) 
    else:
        return len(os.popen("head -n 2 %s" % input_fastq).readlines()[1].rstrip())

def run_cmd ( cmd , wd_tmpdir, descriptor):
    #print cmd
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

def __main__():
    descr = ""
    parser = optparse.OptionParser(description=descr)
    parser.add_option( '-s', '--fileSource', dest="fileSource", help='is the reference file in the history or indexed' )
    parser.add_option( '-r', '--ref', dest="reference_path", help='path to reference' )
    parser.add_option( '-i', '--input-assembly', dest="input_assembly", help='The assembly file' )
    parser.add_option( '', '--output-snp', dest="output_snp", help="The file to save the output snps" )
    parser.add_option( '', '--output-sv', dest="output_sv", help="The file to save the output sv" )
    (options, args) = parser.parse_args()

    if len(args) > 0:
        parser.error('Wrong number of arguments')

    # make temp directory 
    output_dir = tempfile.mkdtemp(prefix="optimized-")
    tmp_dir = tempfile.mkdtemp(prefix="optimized-tmp-")

    # index if necessary
    if options.fileSource == 'history':
        ref_file = tempfile.NamedTemporaryFile( dir=tmp_dir )
        ref_file_name = ref_file.name
        ref_file.close()
        os.symlink( options.ref, ref_file_name )
        # determine which indexing algorithm to use, based on size
        try:
            size = os.stat( options.ref ).st_size
            if size <= 2**30:
                indexingAlg = 'is'
            else:
                indexingAlg = 'bwtsw'
        except:
            indexingAlg = 'is'
        indexing_cmds = '-a %s' % indexingAlg
        cmd1 = 'bwa index %s %s' % ( indexing_cmds, ref_file_name )
        try:
            tmp = tempfile.NamedTemporaryFile( dir=tmp_dir ).name
            tmp_stderr = open( tmp, 'wb' )
            proc = subprocess.Popen( args=cmd1, shell=True, cwd=tmp_dir, stderr=tmp_stderr.fileno() )
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
            if os.path.exists( tmp_dir ):
                shutil.rmtree( tmp_dir )
            raise Exception, 'Error indexing reference sequence. ' + str( e )
    else:
        ref_file_name = options.reference_path

    # create the meta dictionary of the commands to run
    # run the command from the output_dir path
    os.chdir(output_dir) 

    # create a link to the input file
    input_link_name = "%s/prefix.mag.gz" % output_dir
    os.symlink(options.input_assembly, input_link_name)
    wf_meta = {}
    output_prefix = "%s/prefix" % output_dir
    make_file = "%s.mak" % output_prefix
    cmd = "run-calling -t%s %s %s | sh" % (get_ncores(), ref_file_name, input_link_name)
    wf_meta['1'] = [{'cl' : cmd}]

    output_gz = "%s/prefix.flt.vcf.gz" % (output_dir)
    cmd = "gunzip %s" % (output_gz)
    wf_meta['2'] = [{'cl' : cmd}]

    output_vcf = "%s/prefix.flt.vcf" % (output_dir)
    cmd = "mv %s %s" % (output_vcf, options.output_snp)
    wf_meta['3'] = [{'cl' : cmd}]
  
    output_gz = "%s/prefix.sv.vcf.gz" % (output_dir)
    cmd = "gunzip %s" % (output_gz)
    wf_meta['4'] = [{'cl' : cmd}]

    output_vcf = "%s/prefix.sv.vcf" % (output_dir)
    cmd = "mv %s %s" % (output_vcf, options.output_sv)
    wf_meta['5'] = [{'cl' : cmd}]

    #print wf_meta
    for step_id, step_list in sorted(wf_meta.iteritems()):
        print "START STEP %s: %s" % (step_id, time_stamp())
        step_input_files = []
        for job in step_list:
            print "run step %s:\n%s" % (step_id, job['cl'])
            run_cmd ( job['cl'], tmp_dir, "running some job")
        print "END STEP %s: %s" % (step_id, time_stamp())

    #clean up
    shutil.rmtree(output_dir)

if __name__ == "__main__":
    __main__()
