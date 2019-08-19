#/usr/bin/python

import optparse, sys, subprocess, tempfile
import shutil, os

CHUNK_SIZE = 2**20 #1mbi

def cleanup_before_exit( tmp_dir ):
    if tmp_dir and os.path.exists( tmp_dir ):
        shutil.rmtree( tmp_dir )

def __main__():
    parser = optparse.OptionParser()
    parser.add_option( '-p', '--pass_through', dest='pass_through_options', action='append', type="string", help='These options are passed through directly to GATK, without any modification.' )
    parser.add_option('', '--genotype', dest='genotype', action='store', type="string", help='output format')
    parser.add_option('', '--known-haps-g', dest='known_haps_g', action='store', type="string", help='output format')
    parser.add_option('', '--haplotype', dest='haplotype', action='store', type="string", help='output format')
    parser.add_option('', '--legend', dest='legend', action='store', type="string", help='output format')
    parser.add_option('', '--map', dest='maps', action='store', type="string", help='output format')
    parser.add_option( '', '--output-dir', dest='output_dir', action='store', type="string", default=None, help='If specified, the output directory for extra files.' )
    parser.add_option('', '--output', dest='output_impute', action='store', type="string", help='output format')
    parser.add_option('', '--output1', dest='output_genofile', action='store', type="string", help='output format')
    parser.add_option('', '--output2', dest='output_info', action='store', type="string", help='output format')
    parser.add_option('', '--output3', dest='output_summary', action='store', type="string", help='output format')
    (options, args) = parser.parse_args()

    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)
    output_dir = "%s/output" % options.output_dir
    os.mkdir(output_dir)

    pass_through = ""
    if options.pass_through_options:
        pass_through = ' '.join( options.pass_through_options )
    
    #tmp_dir = tempfile.mkdtemp( prefix='tmp-tool-' )
    
    if options.genotype:
        cmd = "impute2 -m %s -g %s %s -o prephasing.impute2" % (options.maps, options.genotype, pass_through)
        print cmd
        stdout = tempfile.NamedTemporaryFile( prefix="unphased-stdout-", dir=output_dir )
        stderr = tempfile.NamedTemporaryFile( prefix="unphased-stderr-", dir=output_dir )
        
    if options.known_haps_g:
        cmd = "impute2 -m %s -h %s -l %s -known_haps_g %s %s -o phased.impute2" % (options.maps, options.haplotype, options.legend, options.known_haps_g, pass_through)
        print cmd
        stdout = tempfile.NamedTemporaryFile( prefix="phased-stdout-", dir=output_dir )
        stderr = tempfile.NamedTemporaryFile( prefix="phased-stderr-", dir=output_dir )


    proc = subprocess.Popen( args=cmd, stdout=stdout, stderr=stderr, shell=True, cwd=output_dir )
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
    stderr.close()
    stdout.close()

    if options.genotype:
        shutil.copy("%s/prephasing.impute2_haps" % output_dir, options.output_impute)
        shutil.copy("%s/prephasing.impute2_info" % output_dir, options.output_info)
        shutil.copy("%s/prephasing.impute2_summary" % output_dir, options.output_summary)   
    if options.known_haps_g:
        shutil.copy("%s/phased.impute2" % output_dir, options.output_impute)
        shutil.copy("%s/phased.impute2_info" % output_dir, options.output_info)
        shutil.copy("%s/phased.impute2_summary" % output_dir, options.output_summary)    

if __name__ == "__main__": __main__()    
