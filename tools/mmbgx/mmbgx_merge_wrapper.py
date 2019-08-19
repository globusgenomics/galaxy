#!/usr/bin/env python


import optparse, os, sys, glob, subprocess, tempfile, shutil

def stop_err( msg ):
    sys.stderr.write( '%s\n' % msg )
    sys.exit()

def create_input_rscript(inputF, output_directory, outputR, outputMuave):
    inputs_string = "\"" + "\", \"".join(inputF) + "\""
    inputs_len = len(inputF)

    fh_R = open("%s/mmbgx_merge_outputs.R" % os.getcwd(), "w")
    fh_R.write("library(mmbgx)\n")
    fh_R.write("dir.names <- dir(c(%s), full.names=TRUE)\n" % (inputs_string))
    fh_R.write("combineRuns.mmbgx(dir.names, sampleSets=rep(1,%s), outputDir=\"%s\")\n" % (inputs_len, output_directory))
    fh_R.write("resT <- readSingle.mmbgx(\"%s\")\n" % (output_directory))
    fh_R.write("Tr<-resT$muave\n")
    fh_R.write("rownames(Tr)<- resT$Psids\n")
    fh_R.write("save(resT, file=\"%s\")\n" % outputR)
    fh_R.write("write(Tr, file=\"%s\")\n" % outputMuave)
    fh_R.close()
    return("%s/mmbgx_merge_outputs.R" % os.getcwd())


def __main__():
    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '-i', '--input', dest='inputF', action="append", help='The input directories' )
    parser.add_option( '-o', '--output', dest='outputR', help='The output R file' )
    parser.add_option( '-d', '--output-directory', dest='output_directory', help='The output directory' )
    parser.add_option( '', '--output-muave', dest='outputMuave', help='The output Muave file' )
    ( options, args ) = parser.parse_args()

    # create the output directory
    if not os.path.isdir(options.output_directory):
        os.mkdir(options.output_directory)

    # create Rscript and run
    outputDirectory = "%s/combinedRuns" % options.output_directory
    rscript_path = create_input_rscript(options.inputF, outputDirectory, options.outputR, options.outputMuave)
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

if __name__=="__main__": __main__()
