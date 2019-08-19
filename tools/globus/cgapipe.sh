#!/bin/bash
#
# Command pipe to generate BAM and FASTQ files from a cgatools part. The
# first two arguments must be the bam and fastq output file locations, and
# the rest of the arguments are passed to cgatools map2sam. For example:
#
#  $ cgapipe.sh output output.fastq --reads=reads.tsv \
#               --mappings=mappings.tsv --reference=ref.crr
#
# Note: .bam will be appended to the bam output filename.
#

if [ $# -lt 3 ]; then
    echo "Usage: $0 bam_outfile fastq_outfile MAP2SAM_ARGS..."
fi

DD_BS="1M"
output_bam="$1"
output_fastq="$2"
shift 2

cd $(dirname $0)
bindir=$(pwd)

export PATH=$bindir:/usr/local/bin:/usr/bin:/bin

cgatools map2sam --add-mate-sequence "$@" \
 | tee \
  >(/nfs/software/bin/picard_sam2fastq.sh \
    QUIET=TRUE INPUT=/dev/stdin FASTQ=/dev/stdout \
    VALIDATION_STRINGENCY=LENIENT 2>/dev/null \
    | dd bs=$DD_BS > "$output_fastq") \
 | samtools view -uS - \
 | samtools sort - "$output_bam"
