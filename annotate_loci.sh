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

if [ $# -lt 3 ]
then
    echo "USAGE: $0 loci_bed genome_gff class_pri" >&2
    exit 1
fi

IN=$1
GFF=$2
CLASSPRI=$3

OUT=${1%.bed}.annot

if [ -z $TMPDIR ]; then
    export TMPDIR=/tmp
fi

cp $IN ${TMPDIR}/$$.in
cp $GFF ${TMPDIR}/$$.gff

# add antisense annotations and append them to sense
cat ${TMPDIR}/$$.gff | awk 'BEGIN{OFS="\t"}
  { if($7 == "+")
    { $7 = "-" }
    else { $7 = "+" }
    $3 = "as-"$3; print $0}' > ${TMPDIR}/$$.gff2
cat ${TMPDIR}/$$.gff ${TMPDIR}/$$.gff2 > ${TMPDIR}/$$.gff3
mv ${TMPDIR}/$$.gff3 ${TMPDIR}/$$.gff

# sort inputs
cat ${TMPDIR}/$$.in | sort -k1,1 -k2,2n > ${TMPDIR}/$$.sorted.bed
cat ${TMPDIR}/$$.gff | sort -k1,1 -nk4,4 > ${TMPDIR}/$$.sorted.gff

intersectBed -s -wo -a ${TMPDIR}/$$.sorted.gff \
  -b ${TMPDIR}/$$.sorted.bed > ${TMPDIR}/$$.intersect

cat ${TMPDIR}/$$.intersect | sort -k13,13 -k16,16nr | \
  awk '{print $13"\t"$3"\t"$16"\t"$9}' > ${TMPDIR}/$$.intersect2 

cut -f4 ${TMPDIR}/$$.sorted.bed | sort | \
  awk '{print $0"\tintergenic\t-1\t."}' | \
  sort -k1 -m - ${TMPDIR}/$$.intersect2 > ${TMPDIR}/$$.intersect3

# Run the class prioritization
# and split snoRNAs into 3 subclasses
if [ -e `dirname $GFF`/frnadb_snorna.txt ]; then
  use_frnadb=`dirname $GFF`/frnadb_snorna.txt
fi


annotate_smrna_loci ${TMPDIR}/$$.intersect3 $CLASSPRI | 
  sort -k1 | \
  annotate_snornas.sh $use_frnadb > $OUT

rm ${TMPDIR}/$$.*
