#!/usr/bin/env python
# vasa

"""
A wrapper script for running swift on globus-galaxy
"""

## MODULES
import sys, optparse, os, tempfile, subprocess, shutil, glob

## GLOBALS    
swift_bin = '/mnt/galaxyTools/tools/swift/swift-0.94.1/bin/swift'
atlas_bin = '/mnt/galaxyTools/tools/ruby/1.9.3-p488/bin/ruby /opt/galaxy/tools/atlas2/Atlas2_v1.4.3/Atlas-SNP2/Atlas-SNP2.rb'
sites_file = '/opt/galaxy/tools/swift/sites.xml'
tc_file = '/opt/galaxy/tools/swift/tc.data'
swift_file = '/opt/galaxy/tools/swift/atlas-snp2_v2_nonArray.swift'
CHUNK_SIZE = 2**20 #1mb

#Parse Command Line     
parser = optparse.OptionParser()
parser.add_option('--bam-dir', dest='bam_dir', type="string", help='Input directory where all BAM file are located' )
parser.add_option('--output-dir', dest='outputdir', type="string", help='Output directory to store VCF files' )
parser.add_option('--fasta', dest='fasta', type="string", help='Fasta file' )

## exposing additional command line arguments
parser.add_option('--target-region', dest='target_region', type="string", help='Target region file.')
parser.add_option('--qual-filter', dest='qual_filter', action="store_true", help='Include filtered lines in the output that have a QUAL of at least 1')
parser.add_option('--platform', dest='platform', type="choice", choices=["454_FLX", "454_XLR", "Illumina"], help='Platform specific calls.')
parser.add_option('--post-cutoff', dest='post_cutoff', type="string", help='Atlas posterior probability cutoff. See Atlas-SNP2 documentation.')
parser.add_option('--min-coverage', dest='min_coverage', type="string", help='Minimum coverage for a site to be called.')
parser.add_option('--prior-e', dest='prior_e', type="string", help='Prior(error|c) when variant coverage number is above 2 for 454 and Illumina data (Default is 0.1)')
parser.add_option('--prior-l', dest='prior_l', type="string", help='Prior(error|c) when variant coverage number is 1 or 2 for 454 data (Default is 0.9)')
parser.add_option('--base-sub-filter', dest='base_sub_max', type="string", help='maximum percentage of substitution bases allowed in the alignment (Default is 5.0)')
parser.add_option('--base-indel-filter', dest='base_indel_max', type="string", help='maximum percentage of insertion and deletion bases allowed in the alignment (Default is 5.0)')
parser.add_option('--insert-size-filter', dest='insert_size_max', type="string", help='insertion size for pair-end re-sequencing data (optional)"')
parser.add_option('--max-alignments', dest='alignment_max', type="string", help='maximum number of alignments allowed to be piled up on a site (Default is 1024)')
parser.add_option('--sites-list', dest='sites_list', type="string", help='File containing sites will always be included(optional)')
parser.add_option('--only-eval-sites', dest='eval_sites_flag', action="store_true", help='only evaluate sites in the list (use with -a, optional)')
parser.add_option( '-d', '--dummy', dest='dummy_out', type="string")
parser.add_option( '-l', '--log', dest='log_out', type="string")

(opts, args) = parser.parse_args()
print "BASE_SUB_MAX: %s" % opts.base_sub_max
print "ARGS: %s" % args
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
tmp_dir = None
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

        ## construct a command line template and have swift fill in the "parallelizalbe"
        ##  parameters using regexps.

        ## parse mandatory commands
        swift_params = list()
        swift_params.append('-bam_dir=' + opts.bam_dir)
        swift_params.append('-outputdir=' + opts.outputdir)

        ## parse the optional commands
        ## These are parsed explicitly for clarity -- atlas command flags are opaque!
        atlas_params = list()
        if opts.platform: atlas_params.append(" --" + opts.platform)
        if opts.post_cutoff: atlas_params.append(" -c " + opts.post_cutoff)
        if opts.min_coverage: atlas_params.append(" -y " + opts.min_coverage)
        if opts.prior_e: atlas_params.append(" -e " + opts.prior_e)
        if opts.prior_l: atlas_params.append(" -l " + opts.prior_l)
        if opts.base_sub_max: atlas_params.append(" -m " + opts.base_sub_max)
        if opts.base_indel_max: atlas_params.append(" -g " + opts.base_indel_max)
        if opts.insert_size_max: atlas_params.append(" -p " + opts.insert_size_max)
        if opts.alignment_max: atlas_params.append(" -f " + opts.alignment_max)
        if opts.sites_list: atlas_params.append(" -a " + opts.sites_list)
        if opts.eval_sites_flag: atlas_params.append(" -w ")
        if opts.qual_filter: atlas_params.append(" -F ")

        ## initial atlas command template
        atlas_cmd = "%s -r %s -i -n -o %s; echo completed_job > OUTPUT" % (atlas_bin, opts.fasta, ' '.join(atlas_params))

        ## check ruby version
        ruby_cmd = "ruby -v"
        proc = subprocess.Popen( args=ruby_cmd, shell=True )
        returncode = proc.wait()

        ## construct the swift command
        swift_cmd = "%s -sites.file %s -tc.file %s %s " %   (swift_bin, sites_file, tc_file, swift_file)
        cmd = "%s %s %s 2>&1" % (swift_cmd, ' '.join(swift_params), '-atlas_cmd=\"'+atlas_cmd+'\"')
        print cmd
        proc = subprocess.Popen( args=cmd, shell=True, cwd=tmp_dir)
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
stdout.write("BAM\t%s\nVCF\t%s\nREF\t%s\n" % (opts.bam_dir, opts.outputdir, opts.fasta))
stdout.close

