'''
Simple tool to append the string 'chr' to all contig IDs in a vcf file
'''

import optparse
import os
import gzip

## pulls the SM tag from the header of a BAM
def extract_sm(bam, direc):
  bamReader = gzip.open(direc + '/' + bam, 'r')
  for line in bamReader:
    if 'SM' in line:
      fields = line.split()
      sm = [f for f in fields if f.startswith('SM')][0].lstrip('SM:')
      break
  bamReader.close()
  return sm  

## builds a sample ID by the same logic used in the swift script
def construct_id(bam):
  iid = bam.strip('_newRG_sorted.bam') + '_newRG.snp.vcf'
  return iid



def main():

  ## parse command line arguments
  parser = optparse.OptionParser()
  parser.add_option('-v', dest = 'oldvcf', help = 'Filepath for input vcf file.')
  parser.add_option('-o', dest = 'newvcf', help = 'Filepath for output vcf file.')
  parser.add_option('-d', dest = 'bamdir', help = 'Directory containing source BAM files')
  (opts, args) = parser.parse_args()

  ## open file connections
  oldvcf = open(opts.oldvcf)
  newvcf = open(opts.newvcf, 'w')

  ## match SM tags to file names
  sampleIds = dict()
  for (path, dirs, files) in os.walk(opts.bamdir):
    for bam in files:
      if bam.endswith('.bam'):
        sm = extract_sm(bam, opts.bamdir)
        oldId = construct_id(bam)
        sampleIds[oldId] = sm

  ## copy over vcf header
  head = oldvcf.next()
  newvcf.write(head)
  while head.startswith('##'):
      head = oldvcf.next()
      if head.startswith('#CHROM'):
        break
      newvcf.write(head)

  ## recode the sample Ids
  line = [ sampleIds[f] if sampleIds.get(f) else f for f in head.split() ]
  newvcf.write('	'.join(line) + '\n')

  ## recode the chromosome IDs in the file
  for line in oldvcf:
      newvcf.write("chr"+line)

  ## close file connections
  oldvcf.close()
  newvcf.close()

if __name__ == '__main__':
  main()

