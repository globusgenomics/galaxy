#!/usr/bin/env python

import optparse, sys, subprocess, tempfile
import shutil, os, glob, re
from threading import Timer

CHUNK_SIZE = 2**20 #1mb


def cleanup_before_exit( tmp_dir ):
    if tmp_dir and os.path.exists( tmp_dir ):
        shutil.rmtree( tmp_dir )

def generate_wiggle(out_dir, out_name):
    with open("%s/%s" % (out_dir,out_name)) as f1:
        out_wig=("%s/%s_mean_coverage.wig" % (out_dir, out_name))
        with open(out_wig, 'w') as out:
            out.write('track type=wiggle_0 name="All Samples Mean Coverage"\n')
            for line in f1:
                if re.match(r'^chr', line):
                    match=re.match(r'(\w+):(\d+)\t\d+\t([-+]?\d*\.\d+|\d+).*', line)
                    outname=' '.join([match.group(1),match.group(2),match.group(2),match.group(3)])
                    out.write(outname)
                    out.write("\n")

def createOutputHTML (outputF,outdir):
    ofh = open(outputF, "w")
    ofh.write( '<html>\n<head>\n<title>IGV snapshot</title>\n</head>\n<body>\n<p/>\n<ul>\n' )
    fileNames=glob.glob("%s/*.tdf" % outdir)
    #outputDir='%s_files' % ''.join(outputF.split('.')[:-1])
    #tdf_fileNames=glob.glob("%s/igv_output/*.tdf" % outputDir)
    #png_fileNames=glob.glob("%s/igv_output/*.png" % outputDir)
    for name in fileNames:
        values = name.split("/")
        sn = values[len(values)-1]
        #outVCF = "%s/%s.vcf" % (outputDir, sn)
        #print "\noutVCF: %s\n" % outVCF
#        if os.path.exists(outVCF):
        ofh.write('<li><a href="%s">%s</a></li>\n' % (name, sn ) )
    ofh.write( '</ul>\n</body>\n</html>\n' )
    ofh.close()

def __main__():
    gatk_jar_path = "/mnt/galaxyTools/tools/gatk3/3.3/GenomeAnalysisTK.jar"
    igvtools_path = "/mnt/galaxyTools/tools/igvtools/12.19.2016/igvtools"
    igv_jar_path = "/mnt/galaxyTools/tools/igv/2.3.88/igv.jar"
    xvfb_run = "/mnt/galaxyTools/tools/xvfb/08-22-2018/xvfb-run"
    parser = optparse.OptionParser()
    parser.add_option('', '--inputbampath', dest='input_bam_path', action='store', type="string", help='output format')
    parser.add_option('', '--inputbed', dest='input_bed', action='store', type="string", help='output format')
    parser.add_option('', '--reference', dest='reference', action='store', type="string", help='output format')
    parser.add_option('', '--output', dest='output', action='store', type="string", help='output format')
    parser.add_option( '', '--output-dir', dest='output_dir', action='store', type="string", default=None, help='If specified, the output directory for extra files.' )
    (options, args) = parser.parse_args()

    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)
    input_dir = "%s/input" % options.output_dir
    gatk_dir = "%s/gatk_output" % options.output_dir
    igv_dir = "%s/igv_output" % options.output_dir
    if not os.path.exists(input_dir):
        os.mkdir(input_dir)
    if not os.path.exists(gatk_dir):
        os.mkdir(gatk_dir)
    if not os.path.exists(igv_dir):
        os.mkdir(igv_dir)

    if options.input_bam_path:
        inputbams = glob.glob(options.input_bam_path + "/*.bam")
        with open("%s/bam.list" % input_dir,"w") as fout:
            for bam in inputbams:
                fout.write("%s\n" % bam)
    with open(options.input_bed) as f, open("%s/IGV_snapshot_batch.txt" % igv_dir, "w") as f2:
        for line in f:
            line=line.strip()
            if re.match(r'^chr', line):
                match=re.match(r'(.*)\t(.*)\t(.*)\t(.*)\t(.*)', line)
                #match=re.match(r'(\w+)\t(\d+)\t(\d+)\t\w+\t.*', line)
                outname='_'.join([match.group(1),match.group(2),match.group(3)])
                region=' '.join([match.group(1),match.group(2),match.group(3)])
                svtype=match.group(4)
                bamname=match.group(5)
                snapshot='_'.join([bamname,outname,svtype,"IGV_Snapshot.png"])
                tmp='-'.join([match.group(2), match.group(3)])
                interval=':'.join([match.group(1), tmp]) 
                logname="%s_depthofcoverage.log" % outname
                #print "new\ngenome hg19\nload %s/%s_post.bam,%s/%s_mean_coverage.tdf\nsnapshotDirectory %s\ngoto %s\nregion %s\nviewaspairs\nsort base\nsnapshot %s\n" % (options.input_bam_path, bamname, output_dir, outname, output_dir, interval, region , snapshot)
          
                f2.write("new\ngenome hg19\nload %s/%s.bam,%s/%s_mean_coverage.tdf\nsnapshotDirectory %s\ngoto %s\nregion %s\nviewaspairs\nsort base\nsnapshot %s\n" % (options.input_bam_path,bamname,igv_dir, outname, igv_dir, interval, region , snapshot))
                cmd = "java -Xmx2g -jar %s -T DepthOfCoverage -R %s -o %s/%s -I %s/bam.list -L %s --log_to_file %s" % (gatk_jar_path,options.reference, gatk_dir, outname, input_dir, interval, logname)  
                stdout = tempfile.NamedTemporaryFile( prefix="tool-stdout-", dir=gatk_dir )
                stderr = tempfile.NamedTemporaryFile( prefix="tool-stderr-", dir=gatk_dir )
               
                print "cmd: %s" % cmd

                proc = subprocess.Popen( args=cmd, stdout=stdout, stderr=stderr, shell=True, cwd=gatk_dir )
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

                # step 2: generate wiggle file
                generate_wiggle(gatk_dir, outname)
                cmd = "%s toTDF %s/%s_mean_coverage.wig %s/%s_mean_coverage.tdf hg19" % (igvtools_path, gatk_dir, outname, igv_dir, outname)
                stdout = tempfile.NamedTemporaryFile( prefix="tool-stdout-", dir=igv_dir )
                stderr = tempfile.NamedTemporaryFile( prefix="tool-stderr-", dir=igv_dir )

                print "cmd: %s" % cmd

                proc = subprocess.Popen( args=cmd, stdout=stdout, stderr=stderr, shell=True, cwd=igv_dir )
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
 
    kill = lambda process: process.kill()
    cmd = "%s --server-args='-screen 0 1920x1080x24' java -Xmx2g -jar %s -b %s/IGV_snapshot_batch.txt -g hg19" % (xvfb_run, igv_jar_path, igv_dir)
    stdout = tempfile.NamedTemporaryFile( prefix="tool-stdout-", dir=igv_dir )
    stderr = tempfile.NamedTemporaryFile( prefix="tool-stderr-", dir=igv_dir )

    print "cmd: %s" % cmd

    proc = subprocess.Popen( args=cmd, stdout=stdout, stderr=stderr, shell=True, cwd=igv_dir )
    
    #return_code = proc.wait()

    #if return_code:
    #    stderr_target = sys.stderr
    #else:
    #    stderr_target = sys.stdout
    #stderr.flush()
    #stderr.seek(0)

    #while True:
    #    chunk = stderr.read( CHUNK_SIZE )
    #    if chunk:
    #        stderr_target.write( chunk )
    #    else:
    #        break
    #stderr.close()
    #stdout.close()
    
    my_timer = Timer(240, kill, [proc])

    try:
        my_timer.start()
        stdout, stderr = proc.communicate()
    finally:
        my_timer.cancel()

    createOutputHTML(options.output, igv_dir)                
                #with open("%s/IGV_snapshot_batch.txt" % output_dir, "w") as f2:
                #print "new\ngenome hg19\nload %s/%s_post.bam,%s/%s_mean_coverage.tdf\nsnapshotDirectory %s\ngoto %s\nregion %s\nviewaspairs\nsort base\nsnapshot %s\n" % (options.input_bam_path, bamname, output_dir, outname, output_dir, interval, region , snapshot)
                #f2.write("new\ngenome hg19\nload %s/%s_post.bam,%s/%s_mean_coverage.tdf\nsnapshotDirectory %s\ngoto %s\nregion %s\nviewaspairs\nsort base\nsnapshot %s\n") % (options.input_bam_path,bamname,output_dir, outname, output_dir, interval, region , snapshot)
         
                #with open("%s/%s" % (output_dir,outname)) as f1:
                #    out_wig=("%s/%s_mean_coverage.wig" % (output_dir, outname))
                #    with open(out_wig, 'w') as out:
                #        out.write('track type=wiggle_0 name="All Samples Mean Coverage"\n')
                #        for line in f1:
                #            if re.match(r'^chr', line):
                #                #match=re.match(r'(\w+):(\d+)\t\d+\t(\d+).*', line)
                #                match=re.match(r'(\w+):(\d+)\t\d+\t([-+]?\d*\.\d+|\d+).*', line)
                #                outname=' '.join([match.group(1),match.group(2),match.group(2),match.group(3)])
                #                out.write(outname)
                #                out.write("\n")
                          
    #shutil.copy("%s/%s" % (output_dir,stdout),options.output)
if __name__ == "__main__": __main__()
