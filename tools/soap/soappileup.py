"""
soappileup.py
A wrapper script for SOAPpileup
Copyright Peter Li peter@gigasciencejournal.com
"""

import sys, optparse, os, tempfile, shutil, subprocess

def stop_err(msg):
    sys.stderr.write('%s\n' % msg)
    sys.exit()


def cleanup_before_exit(tmp_dir):
    if tmp_dir and os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

def __main__():
    #Parse command line
    parser = optparse.OptionParser()
    #Generic input params
    parser.add_option("-r", "--ref", dest="ref", help="A reference genome in FASTA format")
    parser.add_option("-s", "--sorted_soap", dest="sorted_soap", help="A sorted SOAP results file")
    parser.add_option("-l", "--regions", dest="regions", help="Regions to undergo the pileup process")

    parser.add_option("", "--default_full_settings_type", dest="default_full_settings_type")
    #Custom params
    parser.add_option("-d", "--discard_uncovered_loci", dest="discard_uncovered_loci")
    parser.add_option("-c", "--regions_sorted", dest="regions_sorted", help="Regions are already sorted and combined (not overlapped with others)")
    parser.add_option("-e", "--read_edge_length", dest="read_edge_length", help="x bp towards 5' or 3' end of reads would be defined as edge of reads")
    parser.add_option("-m", "--region_limited", dest="region_limited", help="region limited (region will be split into smaller ones)")

    #Outputs
    parser.add_option("-o", "--output", dest='output', help="Results from pileup process")
    opts, args = parser.parse_args()

    #Create temp directory
    tmp_dir = tempfile.mkdtemp()
    #tmp_output = tempfile.NamedTemporaryFile(dir=tmp_dir).name+".gz"

    #Set up command line call
    if opts.default_full_settings_type == "default":
        cmd = "SOAPpileup -r %s -l %s -s %s -o %s -d" % (opts.ref, opts.regions, opts.sorted_soap, opts.output)
    elif opts.default_full_settings_type == "full":
        cmd = "SOAPpileup -r %s -l %s -s %s -o %s -c %s -e %s -m %s -d %s" % (opts.ref, opts.regions, opts.sorted_soap, opts.output, opts.regions_sorted, opts.read_edge_length, opts.region_limited, opts.discard_uncovered_loci)

    print cmd

    #Run
    try:
        tmp_out = tempfile.NamedTemporaryFile(dir=tmp_dir).name
        tmp_stdout = open(tmp_out, 'wb')
        tmp_err = tempfile.NamedTemporaryFile().name
        tmp_stderr = open(tmp_err, 'wb')
        proc = subprocess.Popen( args=cmd, shell=True, cwd=".", stdout=tmp_stdout, stderr=tmp_stderr )
        returncode = proc.wait()
        tmp_stderr.close()
        # get stderr, allowing for case where it's very large
        tmp_stderr = open(tmp_err, 'rb')
        stderr = ''
        buffsize = 1048576
        try:
            while True:
                stderr += tmp_stderr.read(buffsize)
                if not stderr or len(stderr) % buffsize != 0:
                    break
        except OverflowError:
            pass
        tmp_stdout.close()
        tmp_stderr.close()
        if returncode != 0:
            raise Exception, stderr

            # TODO: look for errors in program output.
    except Exception, e:
        #Clean up temp files
        cleanup_before_exit(tmp_dir)
        stop_err('Error in running SOAPpileup from (%s), %s' % (tmp_output, str(e)))
       
    #pass the tmp_output to opts.output
#    gz_output = open(opts.output, 'wb')
#    with open(tmp_output, mode='rb') as file: # b is important -> binary
#        fileContent = file.read()
#        gz_output.write(fileContent)
#    gz_output.close()
#    file.close()
    
    #Clean up temp files
    cleanup_before_exit(tmp_dir)
    #Check results in output file
    if os.path.getsize(opts.output) > 0:
        sys.stdout.write('Status complete')
    else:
        stop_err("The output is empty")

if __name__ == "__main__": __main__()
