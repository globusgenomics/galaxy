#!/bin/bash
#
# Second, per chromosome phase of the QC pipeline, with the fastq conversion
# done in the first per part phase to increase parallelism.
#
# Merges fastq files for a single chromosome, then generates stats and quality
# plots from the data.
#
# Assumes *chrm_N convention for chromosome files in input_dir. Outputs files
# with a chrmN prefix.
#

if [ $# -ne 3 ]; then
    echo "Usage: $0 chrm_number input_dir output_dir"
fi

cd $(dirname $0)
bindir=$(pwd)

export PYTHONPATH=/nfs/software/galaxy/lib

export PATH=$bindir:/usr/local/bin:/usr/bin:/bin:/usr/local/tools/weblogo:/usr/local/tools/blat:/usr/local/tools/homer/bin

chrm="$1"
indir="$2"
outdir="$3"

fastq_out="$outdir/chrm$chrm.fastq"
solexa_out="$outdir/chrm$chrm.solexa.fastq"
stats_out="$outdir/chrm$chrm.stats"
box_out="$outdir/chrm${chrm}_boxplot.png"
dist_out="$outdir/chrm${chrm}_distribution.png"

echo "join chrm $(date)" \
&& cat $indir/*chrm_${chrm}.fastq | dd bs=1M > "$fastq_out" \
&& echo "solexa $(date)" \
&& fastq_groomer.py "$fastq_out" sanger "$solexa_out" \
     solexa ascii summarize_input \
&& echo "stats $(date)" \
&& fastx_quality_stats -i "$solexa_out" -o "$stats_out" \
&& echo "graphs $(date)" \
&& (fastq_quality_boxplot_graph.sh -i "$stats_out" -o "$box_out" \
    -t "chr Fastq Quality Boxplot" &) \
&& (fastx_nucleotide_distribution_line_graph.sh -i "$stats_out" \
    -o "$dist_out" -t "chr Fastq Nucleotide Distribution" &)

wait
echo "done $(date)"
