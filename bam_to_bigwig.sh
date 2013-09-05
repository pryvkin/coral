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


if [ $# -lt 3 ]; then
    echo "USAGE: $0 in_bam out_bigwig_prefix track_name" >&2
    exit 1
fi

BAMIN=$1
BIGWIGPREFIX=$2
TRACKNAME=$3

# generate table of chromosomes and their lengths from the SAM header
samtools view -H $BAMIN | \
  awk 'BEGIN{OFS="\t"}
      /@SQ/ { gsub("SN:", "", $2)
              gsub("LN:", "", $3)
              print $2, $3}' | \
  sort | uniq > $TMPDIR/$$.chromsizes

# make bedgraphs
genomeCoverageBed -split -bg -ibam $BAMIN -g $TMPDIR/$$.chromsizes \
  -strand "+" > $TMPDIR/$$.pos.bedgraph
genomeCoverageBed -split -bg -ibam $BAMIN -g $TMPDIR/$$.chromsizes \
  -strand "-" > $TMPDIR/$$.neg.bedgraph

bedGraphToBigWig $TMPDIR/$$.pos.bedgraph $TMPDIR/$$.chromsizes $TMPDIR/$$.pos.bigwig
bedGraphToBigWig $TMPDIR/$$.neg.bedgraph $TMPDIR/$$.chromsizes $TMPDIR/$$.neg.bigwig

cat - > $TMPDIR/$$.pos.track <<EOF
track type=bigWig name="${TRACKNAME},Pos" description="${TRACKNAME}, Pos strand" bigDataUrl=URLURLURL color=0,0,192 group="${TRACKNAME}" priority=1
EOF

cat - > $TMPDIR/$$.neg.track <<EOF
track type=bigWig name="${TRACKNAME},Neg" description="${TRACKNAME}, Neg strand" bigDataUrl=URLURLURL color=192,0,0 group="${TRACKNAME}" priority=2
EOF

mv $TMPDIR/$$.pos.bigwig ${BIGWIGPREFIX}.pos.bigwig
mv $TMPDIR/$$.pos.track ${BIGWIGPREFIX}.pos.track
mv $TMPDIR/$$.neg.bigwig ${BIGWIGPREFIX}.neg.bigwig
mv $TMPDIR/$$.neg.track ${BIGWIGPREFIX}.neg.track

rm $TMPDIR/$$.*
