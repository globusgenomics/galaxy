#/usr/bin/python

import optparse, sys, subprocess, tempfile
import shutil, os

CHUNK_SIZE = 2**20 #1mbi

def cleanup_before_exit( tmp_dir ):
    if tmp_dir and os.path.exists( tmp_dir ):
        shutil.rmtree( tmp_dir )

def __main__():
    parser = optparse.OptionParser()
    parser.add_option('', '--input-dir', dest='input_dir', action='store', type="string", help='')
    parser.add_option('', '--samplesheet', dest='samplesheet', action='store', type="string", help='')
    parser.add_option( '', '--output-dir', dest='output_dir', action='store', type="string", default=None, help='If specified, the output directory for extra files.' )

    (options, args) = parser.parse_args()

    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)
    output_dir = options.output_dir

    cmd = "cellranger mkfastq  --run=%s --samplesheet=%s --output-dir=%s " % (options.input_dir, options.samplesheet, output_dir)
    print cmd
    stdout = tempfile.NamedTemporaryFile( prefix="cellranger-mkfastq-stdout-", dir=output_dir )
    stderr = tempfile.NamedTemporaryFile( prefix="cellranger-mkfastq", dir=output_dir )


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
