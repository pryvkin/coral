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

# compute length vector features for all loci

if [ $# -lt 2 ]; then
    echo "USAGE: $0 bam_file config_file" >&2
    exit 1
fi

bam=$1
config=$2

source $config

outdir=`dirname $bam`/coral
mkdir -p $outdir

###
echo "Computing genomic length vectors..." >&2

compute_genomic_lenvectors $bam $outdir/genomic_lenvec.plus \
  $outdir/genomic_lenvec.minus $min_trimmed_read_len $max_trimmed_read_len

echo "Indexing genomic length vectors..." >&2

index_genomic_lenvectors $outdir/genomic_lenvec.plus
index_genomic_lenvectors $outdir/genomic_lenvec.minus

echo "Computing locus length features..." >&2

compute_locus_lenvectors $outdir/loci.bed $outdir/genomic_lenvec \
  $min_trimmed_read_len > $outdir/feat_lengths.txt




