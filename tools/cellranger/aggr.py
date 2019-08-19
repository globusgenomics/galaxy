#/usr/bin/python

import optparse, sys, subprocess, tempfile
import shutil, os

CHUNK_SIZE = 2**20 #1mbi

def __main__():
    parser = optparse.OptionParser()
    parser.add_option('', '--id', dest='id', action='store', type="string", help='')
    parser.add_option('', '--normalize', dest='normalize', action='store', type="string", help='')
    parser.add_option('', '--aggregation', dest='aggregation', action='append', nargs=2, type="string", help='')
    parser.add_option( '-o', '--output', dest='library_csv', action='store', type="string", help='output_batch' )
    parser.add_option('', '--output-dir', dest='output_dir', action='store', type="string", default=None, help='If specified, the output directory for extra files.' )

    (options, args) = parser.parse_args()

    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)
    output_dir = options.output_dir

    f1=open(options.library_csv, 'w+')
    f1.write("library_id,molecule_h5\n")
    for library in options.aggregation:
        f1.write("%s\n" %  ",".join(library))
    f1.close()

    cmd = "cellranger aggr  --id=%s  --normalize=%s --csv=%s" % (options.id, options.normalize, options.library_csv)

    print cmd
    stdout = tempfile.NamedTemporaryFile( prefix="cellranger-aggr-stdout-", dir=output_dir )
    stderr = tempfile.NamedTemporaryFile( prefix="cellranger-aggr-stderr", dir=output_dir )


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

if __name__ == "__main__": __main__()
