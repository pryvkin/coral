#!/bin/bash
#  Copyright (c) 2013 University of Pennsylvania
#
#  Permission is hereby granted, free of charge, to any person obtaining a 
#  copy of this software and associated documentation files (the "Software"), 
#  to deal in the Software without restriction, including without limitation 
#  the rights to use, copy, modify, merge, publish, distribute, sublicense, 
#  and/or sell copies of the Software, and to permit persons to whom the 
#  Software is furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included in 
#  all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
#  DEALINGS IN THE SOFTWARE.


# Annotate small RNA loci using a genome annotation gff
# and a list of class priorities (to resolve cases of 
# multiple annotation overlaps)

if [ $# -lt 2 ]; then
    echo "USAGE: $0 bam_file config_file" >&2
    exit 1
fi

if [ -z $TMPDIR ]; then
    export TMPDIR=/tmp
fi

bam=$1
config=$2

source $config

outdir=`dirname $bam`/coral
mkdir -p $outdir

outprefix=`basename ${1%.bam}`

###
echo "Converting BAM to bigWig..."

bam_to_bigwig.sh $bam $outdir/loci $outprefix


###
echo "Segmenting coverage into transcribed loci..."

segment_bigwig_into_loci.sh \
  $outdir/loci.pos.bigwig $outdir/loci.neg.bigwig \
  $seg_threshold $seg_maxgap $seg_minrun $outdir/loci.bed


###
echo "Counting reads at loci..."

count_reads_at_loci.sh $bam $outdir/loci.bed \
  $min_trimmed_read_len $max_trimmed_read_len > $outdir/loci.cov

# put read counts in 'score' column of bed
cut -f7 $outdir/loci.cov | paste $outdir/loci.bed - | \
  awk 'BEGIN{ OFS="\t" }
       { $5=$7; print }' | cut -f1-6 \
  > $outdir/loci_tmp.bed
mv $outdir/loci_tmp.bed $outdir/loci.bed


