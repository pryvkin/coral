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

# create a .cov file containing the read count at each locus

if [ $# -lt 4 ]; then
    echo "USAGE: $0 bam_file bed_file min_read_len max_read_len" >&2
    exit 1
fi

bam=$1
bed=$2
minlen=$3
maxlen=$4

# only count reads in length range used to generate features
#samtools view -h $bam | \
#  awk '/^@/ || (length($10) >= '$minlen' && length($10) <= '$maxlen')' | \
#  samtools view -Sb - | \
#  coverageBed -s -counts -abam - -b $bed | sort -k4,4 > ${bed%.*}.cov

coverageBed -s -counts -abam $bam -b $bed | sort -k4,4 > ${bed%.*}.cov
