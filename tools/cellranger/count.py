#/usr/bin/python

import optparse, sys, subprocess, tempfile
import shutil, os

CHUNK_SIZE = 2**20 #1mbi

def cleanup_before_exit( tmp_dir ):
    if tmp_dir and os.path.exists( tmp_dir ):
        shutil.rmtree( tmp_dir )

def __main__():
    parser = optparse.OptionParser()
    parser.add_option('', '--id', dest='id', action='store', type="string", help='')
    parser.add_option('', '--fastqs', dest='fastqs_dir', action='store', type="string", help='')
    parser.add_option('', '--transcriptome', dest='transcriptome_dir', action='store', type="string", help='')
    parser.add_option('', '--sample', dest='sample', action='store', type="string", help='')
    parser.add_option('', '--expect-cells', dest='expect_cells', action='store', type="string", help='')
    parser.add_option('', '--chemistry', dest='chemistry', action='store', type="string", help='')

    parser.add_option( '', '--output-dir', dest='output_dir', action='store', type="string", default=None, help='If specified, the output directory for extra files.' )
    parser.add_option('', '--metric', dest='metric', action='store', type="string", help='')
    parser.add_option('', '--hdf5', dest='hdf5', action='store', type="string", help='')
    parser.add_option('', '--mhdf5', dest='mhdf5', action='store', type="string", help='')
    (options, args) = parser.parse_args()

    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)
    output_dir = options.output_dir
    cmd = "cellranger count  --id=%s --transcriptome=%s --fastqs=%s --sample=%s --expect-cells=%s --chemistry=%s" % (options.id, options.transcriptome_dir, options.fastqs_dir, options.sample, options.expect_cells, options.chemistry)

    print cmd
    stdout = tempfile.NamedTemporaryFile( prefix="cellranger-count-stdout-", dir=output_dir )
    stderr = tempfile.NamedTemporaryFile( prefix="cellranger-count-stderr", dir=output_dir )


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

    shutil.copy('%s/%s/outs/metrics_summary.csv' % (output_dir, options.id), options.metric)
    shutil.copy('%s/%s/outs/filtered_gene_bc_matrices_h5.h5' % (output_dir, options.id), options.hdf5)
    shutil.copy('%s/%s/outs/molecule_info.h5' % (output_dir, options.id), options.mhdf5)

if __name__ == "__main__": __main__()
