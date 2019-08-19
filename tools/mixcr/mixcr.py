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
    parser.add_option('', '--input', dest='inputtype', action='store', type="string", help='output format')
    parser.add_option('','--inputfastq', dest='inputfastq', action='store', type="string", help='output format')
    parser.add_option('','--parameters', dest='parameters', action='store', type="string", help='output format')
    parser.add_option('', '--OassemblingFeatures', dest='OassemblingFeatures', action='store', type="string", help='output format')
    parser.add_option('', '--report', dest='report', action='store', type="string", help='output format')
    parser.add_option( '', '--output-dir', dest='output_dir', action='store', type="string", default=None, help='If specified, the output directory for extra files.' )
    parser.add_option('', '--output', dest='output', action='store', type="string", help='output format')
    parser.add_option('', '--output1', dest='clones_txt', action='store', type="string", help='output format')
    parser.add_option('', '--output2', dest='alignments_txt', action='store', type="string", help='output format')

    (options, args) = parser.parse_args()

    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)
    input_dir = "%s/input" % options.output_dir
    if not os.path.exists(input_dir):
        os.mkdir(input_dir)
    output_dir = "%s/output" % options.output_dir
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    pass_through = ""
    if options.pass_through_options:
        pass_through = ' '.join( options.pass_through_options )
    print "options.inputfastq: %s" % options.inputfastq
    
    if options.inputtype == "standard":
        cmd = "mixcr align --report %s %s --parameters %s %s %s/alignments.vdjca; mixcr assemble -OassemblingFeatures=%s --report %s %s/alignments.vdjca %s/clones.clns; mixcr exportClones %s/clones.clns %s/clones.txt; mixcr exportAlignmentsPretty %s/alignments.vdjca %s/alignments.txt" % (options.report, pass_through, options.parameters, options.inputfastq, output_dir, options.OassemblingFeatures, options.report, output_dir, output_dir, output_dir, output_dir, output_dir, output_dir)

    if options.inputtype == "rnaseq":
        cmd = "mixcr align --parameters %s -OallowPartialAlignments=true --report %s %s %s/alignments.vdjca; mixcr assemblePartial --report %s %s/alignments.vdjca %s/alignmentsRescued.vdjca; mixcr assemble -OassemblingFeatures=%s --report %s %s/alignments.vdjca %s/clones.clns; mixcr exportClones %s/clones.clns %s/clones.txt; mixcr exportAlignmentsPretty %s/alignmentsRescued.vdjca %s/alignments.txt" % (options.parameters, options.report, options.inputfastq, output_dir, options.report, output_dir, output_dir, options.OassemblingFeatures, options.report, output_dir, output_dir, output_dir, output_dir, output_dir, output_dir)
        #cmd = "mixcr align --parameters %s -OallowPartialAlignments=true --report %s %s/%s %s/alignments.vdjca; mixcr assemblePartial --report %s %s/alignments.vdjca %s/alignmentsRescued.vdjca; mixcr assemble -OassemblingFeatures=%s --report %s %s/alignments.vdjca %s/clones.clns; mixcr exportClones %s/clones.clns %s/clones.txt; mixcr exportAlignmentsPretty %s/alignmentsRescued.vdjca %s/alignments.txt" % (options.parameters, options.report, input_dir, options.inputfastq, output_dir, options.report, output_dir, output_dir, options.OassemblingFeatures, options.report, output_dir, output_dir, output_dir, output_dir, output_dir, output_dir)

    stdout = tempfile.NamedTemporaryFile( prefix="mixcr-stdout-", dir=output_dir )
    stderr = tempfile.NamedTemporaryFile( prefix="mixcr-stderr-", dir=output_dir )
        
    print "cmd: %s" % cmd

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

    shutil.copy("%s/mixcrReport.log" % output_dir, options.output)
    shutil.copy("%s/clones.txt" % output_dir, options.clones_txt)
    shutil.copy("%s/alignments.txt" % output_dir, options.alignments_txt)   

if __name__ == "__main__": __main__()    
