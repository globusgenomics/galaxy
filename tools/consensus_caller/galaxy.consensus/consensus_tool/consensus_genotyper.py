#from ensemble_walker import concordant_walker
from vcf_ensemble import vcf_ensemble
from variant_ensemble import variant_ensemble
from consensus_writer import consensus_vcf
import argparse as arg
import pysam
import sys, tempfile, os, shutil
import subprocess

CHUNK_SIZE = 2**20 #1mb

def run_cmd ( cmd , wd_tmpdir, descriptor):
    #print cmd
    stderr = tempfile.NamedTemporaryFile( prefix = "cmd_stderr" )
    proc = subprocess.Popen( args=cmd, stderr=stderr, shell=True, cwd=wd_tmpdir )

    exit_code = proc.wait()

    if exit_code:
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


def __main__():

  ## parse command line
  parser = arg.ArgumentParser(description='Find sites and genotypes that aggree among an arbitrary number of VCF files.')
  parser.add_argument('vcfFiles', nargs='+', metavar='VCFS', help='List of VCF files for input.')
  parser.add_argument('--site-threshold', '-s', dest='siteThresh', action='store', type=int, help='Number of inputs which must agree for a site to be included in the output.')
  parser.add_argument('--genotype-threshold', '-g', dest='genoThresh', action='store', type=int, help='Number of inputs which must agree for a genotype to be marked as non-missing.')
  parser.add_argument('--ignore-missing', '-m', dest='ignoreMissing', action='store_true', help='Flag specifying how to handle missing genotypes in the vote. If present, missing genotypes are excluded from the genotype concordance vote unless all genotypes are missing.')
  args = parser.parse_args()

  tmpdir = tempfile.mkdtemp(prefix='consensus_')  ## make tmpdir

  indexedVcfFiles = list()
  for vcf in args.vcfFiles:
    baseVCF = os.path.basename(vcf)
    link_name = "%s/%s" % (tmpdir, baseVCF)
    ##   link vcf in temporary directory
    os.symlink(vcf, link_name)

    ##   create tabix index of vcf in new directory
    # sort
    sorted_vcf = "%s.vcf.sorted" % (link_name)
    sort_cmd = 'vcf-sort %s > %s' % (link_name, sorted_vcf)
    run_cmd(sort_cmd, tmpdir, "sort %s" % (link_name))

    # bgzip
    bgzip_cmd = 'bgzip %s.vcf.sorted' % (link_name)
    run_cmd(bgzip_cmd, tmpdir, "bgzip %s" % (link_name))

    # tabix
    tabix_cmd = 'tabix -f -p vcf %s.vcf.sorted.gz' % (link_name)
    run_cmd(tabix_cmd, tmpdir, "tabix %s" % (link_name))

    ##   append vcf to list stored in indexedVcfFiles list()
    indexedVcfFiles.append('%s.vcf.sorted.gz' % (link_name))
 
  ## create the VCF ensemble
  ensemble = vcf_ensemble(vcfList = indexedVcfFiles, ignoreMissing = args.ignoreMissing)

  ## setup output VCF file. Dummy fields are created for downstream parsing with other tools.
  outVcf = consensus_vcf()
  outVcf.add_format(id="CN", number="1", type="Character", description="Consensus status. \'C\' is concordant, \'D\' is discordant, and \'A\' is ambiguous (i.e. no majority at the given genotype threshold).")
  outVcf.add_format(id="GT", number="1", type="String", description="Genotype")
  outVcf.add_info(id="X1", number="1", type="String", description="Placeholder for INFO parsing")
  outVcf.add_info(id="X2", number="1", type="String", description="Placeholder 2 for INFO parsing")
  outVcf.samples = ensemble.samples
  outVcf.write_header()

  ## iterate over the concordant sites 
  for records, genotypes in ensemble.concordant_variants(siteThresh=args.siteThresh, genoThresh=args.genoThresh):
    outVcf.write_record( records, genotypes )

  #shutil.rmtree(tmpdir)

if __name__ == '__main__':
  __main__()

