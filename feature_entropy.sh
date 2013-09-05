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

# compute read-position-entropy features for all loci

if [ $# -lt 3 ]; then
  echo "USAGE: $0 inbam config_file chromInfo" >&2
  exit 1
fi

bam=$1
config=$2
chrominfo=$3

source $config

outdir=`dirname $bam`/coral
mkdir -p $outdir

##
# NOTE: requires sorted indexed bam
echo "Computing locus read positions..." >&2
rm ${outdir}/locus_readpos.counts
while read line; do
  reg=$(echo $line | awk '{print $1":"($2-1)"-"$3}')
  strand=$(echo $line | awk '{print $6}')
  id=$(echo $line | awk '{print $4}')
  echo $id"..." >&2
  samtools view $bam $reg | \
    awk '("'$strand'" == "+" && (and($2, 16) == 0)) {
           print "'$id'\t5p\t"($4-1)  
           print "'$id'\t3p\t"($4-1+length($10))
         }
         ("'$strand'" == "-" && (and($2, 16) != 0)) {
           print "'$id'\t5p\t"($4-1+length($10))
           print "'$id'\t3p\t"($4-1)
         } '
done < ${outdir}/loci.bed | \
  sort -k1,2 -k3n | uniq -c | \
    awk '{print $2"\t"$3"\t"$1}' > ${outdir}/locus_readpos.counts

compute_locus_entropy.rb ${outdir}/locus_readpos.counts > \
  ${outdir}/feat_posentropy.txt
   
