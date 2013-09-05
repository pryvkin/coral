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

# segments RNAseq data into "transcriptionally active regions"
#
# NOTE: USES A LOT OF RAM! (>4G) due to bgrSegmenter. watch out

if [ $# -lt 6 ]
then
  echo "USAGE: `basename $0` pos.bigwig neg.bigwig threshold maxgap minrun output.bed"
  exit 1
fi

POSBW=$1
NEGBW=$2
THRESHOLD=$3
MAXGAP=$4
MINRUN=$5
OUTBED=$6

bigWigToBedGraph $POSBW $TMPDIR/$$.plus.bgr&
bigWigToBedGraph $NEGBW $TMPDIR/$$.minus.bgr&
wait

# split bedGraph files up by chromosome name
awk '{print > "'${TMPDIR}/$$'.plus."$1"_split.bgr"}' $TMPDIR/$$.plus.bgr&
awk '{print > "'${TMPDIR}/$$'.minus."$1"_split.bgr"}' $TMPDIR/$$.minus.bgr&
wait

for inbgr in `find $TMPDIR/$$.*_split.bgr | sort`; do
  echo `basename ${inbgr}`...
  good=`awk '{tot += $4; lines++} END{if(tot >= '$THRESHOLD' && lines > 1) {
      print 1 } else { print 0 }}' $inbgr`

  if [ $good -eq 0 ]
  then
    echo "skipped."
    continue
  fi

  inprefix=`echo $inbgr | sed -e 's/.bgr$//'`
  strand="+"
  t=`echo basename $inbgr | grep -c minus`
  if [ `echo basename $inbgr | grep -c minus` -eq 1 ]
  then
      strand="-"
  fi
  bgrSegmenter $inprefix $THRESHOLD $MAXGAP $MINRUN  | \
      awk '{print $0"\tx\t0\t'$strand'" }' > \
      ${inprefix}.bed
done

cat $TMPDIR/$$.*_split.bed | \
  awk '{printf "%s\t%s\t%s\tSL%010d\t%s\t%s\n",$1,$2,$3,NR,$5,$6 }' >\
  $TMPDIR/$$.bgrseg.bed

mv $TMPDIR/$$.bgrseg.bed $OUTBED

#rm $TMPDIR/$$.*
