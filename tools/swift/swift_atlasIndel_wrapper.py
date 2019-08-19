#!/usr/bin/env python
# vasa

"""
A wrapper script for running swift on globus-galaxy
"""

## MODULES
import sys, optparse, os, tempfile, subprocess, shutil, glob

## GLOBALS    
swift_bin = '/mnt/galaxyTools/tools/swift/swift-0.94.1/bin/swift'
atlas_bin = '/mnt/galaxyTools/tools/ruby/1.9.3-p488/bin/ruby /opt/galaxy/tools/atlas2/Atlas2_v1.4.3/Atlas-Indel2/Atlas-Indel2.rb'
sites_file = '/opt/galaxy/tools/swift/sites.xml'
tc_file = '/opt/galaxy/tools/swift/tc.data'
swift_file = '/opt/galaxy/tools/swift/atlas-indel_nonArray.swift'
CHUNK_SIZE = 2**20 #1mb

#Parse Command Line     
parser = optparse.OptionParser()
parser.add_option('--bam-dir', dest='bam_dir', type="string", help='Input directory where all BAM file are located' )
parser.add_option('--output-dir', dest='outputdir', type="string", help='Output directory to store VCF files' )
parser.add_option('--fasta', dest='fasta', type="string", help='Fasta file' )

## exposing additional command line arguments
parser.add_option('--type', dest='platform', type="string", help='Platform specific calls.')
parser.add_option('--p_cutoff', dest='p_cutoff', type="string", help='Atlas posterior probability cutoff. See Atlas-Indel2 documentation.')
parser.add_option('--p1b_cutoff', dest='p1b_cutoff', type="string", help='The indel probability (p) cutoff value for 1bp deletions.')
parser.add_option('--bed', dest='bed', type="string", help='Target region file.')
parser.add_option('--orig_qual', dest='orig_qual', type="string", help='Not recommeded for Illumina')
parser.add_option('--norm_qual', dest='norm_qual', type="string", help='Not recommeded for SOLiD')
parser.add_option('--depth', dest='depth', type="string", help='')
parser.add_option('--min_var_reads', dest='min_var_reads', type="string", help='')
parser.add_option('--min_var_ratio', dest='min_var_ratio', type="string", help='')
parser.add_option('--strand_dir_filter', dest='strand_dir_filter', type="string", help='')
parser.add_option('--near_read_end_ratio', dest='near_read_end_ratio', type="string", help='')
parser.add_option('--homo_var_cutoff', dest='homo_var_cutoff', type="string", help='')
parser.add_option('--site_list', dest='site_list', type="string", help='')
parser.add_option( '-d', '--dummy', dest='dummy_out', type="string")
parser.add_option( '-l', '--log', dest='log_out', type="string")

(opts, args) = parser.parse_args()

## construct a command line template and have swift fill in the "parallelizalbe"
##  parameters using regexps.

## parse mandatory commands
swift_params = list()
swift_params.append('-bam_dir=' + opts.bam_dir)
swift_params.append('-outputdir=' + opts.outputdir)
tmp_dir = None

if not os.path.isdir(opts.outputdir):
    os.mkdir(opts.outputdir)

completed_samples = 0
input_samples = {}
input_samples_count = 0
for i in os.listdir(opts.bam_dir):
    if i.endswith("bam"):
        sample_name = os.path.basename(i)
        bam_file = "%s/%s" % (opts.bam_dir, i)
        index_name = "%s/%s.bai" % (opts.bam_dir, i)
        log_name = "%s/%s.log" % (opts.outputdir, sample_name)
        input_samples_count += 1
        input_samples[sample_name] = {'index' : index_name, 'log' : log_name, 'bam' : bam_file}
flag = 0
#while flag==0:
if flag==0:
    # figure out which samples need to be run
    incomplete = []
    for sample in input_samples.keys():
        if os.path.exists(input_samples[sample]['log']):
            completed_samples += 1
        else:
            incomplete.append(input_samples[sample]['bam'])

    if completed_samples == input_samples_count:
        ## all samples are complete
        flag = 1
    else:
        ## samples still need to be analyzed
        ## link bam files to temporary directory for input dataset
        print "Samples Complete: %s" % (completed_samples)
        print "Samples Incomplete: %s" % (input_samples_count-completed_samples)
        print "Submitting Incomplete Samples: %s" % (incomplete)
        tmp_dir = tempfile.mkdtemp()
        print tmp_dir
        for bam in incomplete:
            sample_name = os.path.basename(bam)
            os.symlink(bam, "%s/%s" % (tmp_dir, sample_name))


        ## parse the optional commands
        ## These are parsed explicitly for clarity -- atlas command flags are opaque!
        atlas_params = list()
        if opts.platform: atlas_params.append(" " + opts.platform)
        if opts.p_cutoff: atlas_params.append(" -p " + opts.p_cutoff)
        if opts.p1b_cutoff: atlas_params.append(" -P " + opts.p1b_cutoff)
        if opts.bed: atlas_params.append(" -B " + opts.bed)
        if opts.orig_qual == "Yes": atlas_params.append(" -O ")
        if opts.norm_qual == "Yes": atlas_params.append(" -N ")
        if opts.depth: atlas_params.append(" -t " + opts.depth)
        if opts.min_var_reads: atlas_params.append(" -m " + opts.min_var_reads)
        if opts.min_var_ratio: atlas_params.append(" -v " + opts.min_var_ratio)
        #if opts.strand_dir_filter: atlas_params.append(" -f ")
        if opts.near_read_end_ratio: atlas_params.append(" -n " + opts.near_read_end_ratio)
        if opts.homo_var_cutoff: atlas_params.append(" -h " + opts.homo_var_cutoff)
        if opts.site_list: atlas_params.append(" -a " + opts.site_list)

        ## initial atlas command template
        atlas_cmd = "%s -r %s -b -s -o %s" % (atlas_bin, opts.fasta, ' '.join(atlas_params))

        ## check ruby version
        ruby_cmd = "ruby -v"
        proc = subprocess.Popen( args=ruby_cmd, shell=True )
        returncode = proc.wait()

        ## construct the swift command
        swift_cmd = "%s -sites.file %s -tc.file %s %s " %   (swift_bin, sites_file, tc_file, swift_file)
        cmd = "%s %s %s 2>&1" % (swift_cmd, ' '.join(swift_params), '-atlas_cmd=\"'+atlas_cmd+'\"')
        print cmd
        proc = subprocess.Popen( args=cmd, shell=True, cwd=tmp_dir )
        returncode = proc.wait()

        # delete tmp_dir
        #shutil.rmtree(tmpdir)

swift_log_files = glob.glob("%s/*.log" % tmp_dir)
cmdSummary = "export PYTHONPATH=/mnt/galaxyTools/tools/pymodules/python2.7/lib/python:/opt/galaxy/lib:\$PYTHONPATH;. /mnt/galaxyTools/tools/pymodules/python2.7/env.sh;/opt/galaxy/tools/swift/parse_swift_log.py "
for logF in swift_log_files:
    if "swift.log" in logF:
        continue
    cmdSummary += " -l %s " % logF
cmdSummary += " -o %s" % opts.log_out

return_code = None
stderr = tempfile.NamedTemporaryFile( prefix="TOOL-stderr-", dir=tmp_dir )
stdout = tempfile.NamedTemporaryFile( prefix="TOOL-stdout-", dir=tmp_dir )
if return_code is None or not return_code:
    proc = subprocess.Popen( args=cmdSummary, stdout=stdout, stderr=stderr, shell=True, cwd=tmp_dir )
    return_code = proc.wait()
    if return_code:
        stderr_target = sys.stderr
    else:
        if stdout:
            stderr_target = stdout
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

stdout = open( opts.dummy_out, 'w' )
stdout.write("BAM\t%s\nINDEL\t%s\nREF\t%s\n" % (opts.bam_dir, opts.outputdir, opts.fasta))
stdout.close

