#!/usr/bin/env python


import optparse, os, sys, glob, subprocess, tempfile, shutil

def stop_err( msg ):
    sys.stderr.write( '%s\n' % msg )
    sys.exit()

def create_input_rscript(input_dir):
    fh_R = open("%s/mmbgx_create_inputs.R" % os.getcwd(), "w")
    fh_R.write("library(mmbgx)\n")
    fh_R.write("m<- list.files(path = \"%s\", pattern = \".CEL\", full.names = TRUE)\n" % input_dir)
    fh_R.write("d <- read.celfiles(m,pkgname=\"pd.huex.1.0.st.v2\")\n") 
    fh_R.write("mmbgx(d, arrayType=\"HuEx-1_0-st-v2\", geneLevel=FALSE, inputDirs=paste(\"%s/transcriptInput\", c(1:ncol(d)), sep=\".\"), standalone=TRUE)\n" % (input_dir))
    fh_R.close()
    return("%s/mmbgx_create_inputs.R" % os.getcwd())


def __main__():
    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '', '--input', dest='inputF', help='The input path to CEL file' )
    parser.add_option( '', '--output', dest='output', help='The output file' )
    parser.add_option( '', '--output-directory', dest='output_directory', help='The output directory' )
    ( options, args ) = parser.parse_args()

    # create the output directory
    if not os.path.isdir(options.output_directory):
        os.mkdir(options.output_directory)

    # get temp dir
    #tmp_dir = tempfile.mkdtemp(prefix="mmbgx-tmp-")
    tmp_dir = tempfile.mkdtemp(prefix="mmbgx-tmp-")
    input_dir = "%s/input" % tmp_dir
    output_tmp_dir = "%s/transcriptRuns" % tmp_dir
    os.mkdir(input_dir)
    os.mkdir(output_tmp_dir)

    # get the path of the input CEL from the input file
    fh = open(options.inputF, "r")
    cel_file = None
    for line in fh:
        line = line.rstrip("\n")
        cel_file = line
        break

    if cel_file is None:
        print cel_file
        message = "There is no CEL file to analyze. Enter a valid path."
        sys.exit(message)

    # link the cel file to the tmp input directory
    link_cel = "%s/%s" % (input_dir, os.path.basename(cel_file))
    print link_cel
    print cel_file
    print os.path.basename(cel_file)
    os.symlink(cel_file, link_cel)

    # create Rscript and run
    rscript_path = create_input_rscript(input_dir)
    print rscript_path

    # output version # of tool
    try:
        tmp = tempfile.NamedTemporaryFile().name
        tmp_stdout = open( tmp, 'wb' )
        proc = subprocess.Popen( args='Rscript %s' % rscript_path, shell=True, stdout=tmp_stdout )
        tmp_stdout.close()
        returncode = proc.wait()
        stdout = None
    except:
        sys.stdout.write( 'Could not run directory setup\n' )

    fh_out = open(options.output, "r")
    for idir in glob.glob("%s/transcriptInput.*" % input_dir):
        try:
            command = 'mmbgx-icc %s/infile.txt %s' % (idir,output_tmp_dir) 
            print command
            tmp_cmd_dir = tempfile.mkdtemp()
            tmp = tempfile.NamedTemporaryFile( dir=tmp_cmd_dir ).name
            tmp_stderr = open( tmp, 'wb' )
            proc = subprocess.Popen( args=command, shell=True, cwd=tmp_cmd_dir, stderr=tmp_stderr.fileno() )
            returncode = proc.wait()
            tmp_stderr.close()

            for outdir in glob.glob("%s/run.*" % output_tmp_dir):
                shutil.move(outdir, options.output_directory)
                fh_out.write("OUTPUT: %s" % outdir)
            
        except:
            sys.stdout.write( 'Could not run mmbgx\n' )
    fh_out.close()



if __name__=="__main__": __main__()

