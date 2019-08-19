#!/bin/bash
#
# First, per part phase of QC pipeline that includes fastq conversion, to
# increase the parallelism.
#
# Generates BAM and index files from a cga part, splits the bam by chromosome,
# and converts the chromosome bams to fastq. The first argument must be the bam
# output file location without the .bam extension, the second argument is the
# chromosome base file (_N.bam is appended to each out file, for N 0 to 23),
# and the rest of the arguments are passed to cgatools map2sam.  For example:
#
#  $ qc2-part.sh output/L01_001 output/L01_001-chrm --reads=reads.tsv \
#                --mappings=maps.tsv --reference=build37.crr
#

if [ $# -lt 3 ]; then
    echo "Usage: $0 bam_base_file chrm_base_file MAP2SAM_ARGS..."
fi

cd $(dirname $0)
bindir=$(pwd)

export PATH=$bindir:/usr/local/bin:/usr/bin:/bin:/usr/local/tools/weblogo:/usr/local/tools/blat:/usr/local/tools/homer/bin

bam_out="$1"
chrm_out="$2"
shift 2

echo "cga to bam $(date)" \
&& cgatools map2sam --add-mate-sequence "$@" | dd bs=1M \
    | samtools view -uS - \
    | samtools sort - "$bam_out" \
&& echo "index $(date)" \
&& samtools index "${bam_out}.bam" \
&& echo "split $(date)" \
&& bam splitChromosome --in "${bam_out}.bam" --out "$chrm_out" \
    --bamIndex "${bam_out}.bam.bai" \
&& echo "sam2fastq $(date)" \
&& for chrm in "$chrm_out"_*.bam; do
      fastq_out=${chrm%.bam}.fastq
      (picard_sam2fastq.sh QUIET=TRUE \
       INPUT="$chrm" FASTQ=/dev/stdout \
       VALIDATION_STRINGENCY=LENIENT 2>/dev/null \
       | dd bs=1M > "$fastq_out" &)
   done

wait
echo "done $(date)"
