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

# compute nucleotide frequencies of reads:
# weighted by expression level,
# additively smoothed (+1),
# and normalized to log(odds ratio) versus equal nucleotide frequencies

if [ $# -lt 2 ]; then
  echo "USAGE: $0 inbam config_file" >&2
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
echo "Computing nucleotide frequencies..." >&2

echo -e "name\tnuc_A\tnuc_C\tnuc_G\tnuc_T" > ${outdir}/feat_nuc.txt

cat ${outdir}/loci.bed | \
while read line; do
  reg=$(echo $line | awk '{print $1":"($2-1)"-"$3}')
  strand=$(echo $line | awk '{print $6}')
  id=$(echo $line | awk '{print $4}')
  echo $id"..." >&2
  samtools view $bam $reg | \
    awk '("'$strand'" == "+" && (and($2, 16) == 0)) ||
         ("'$strand'" == "-" && (and($2, 16) != 0)) {
           print "'$id'\t"$10 } ' | sort | uniq -c
done | \
  awk 'function process_counts(locus) {
         s = 4 + count["A"] + count["C"] + count["G"] + count["T"]
         for(x in count) { count[x] = log( (count[x]/s) / 0.25 ) }
         print locus,count["A"],count["C"],count["G"],count["T"]
       }
       function reset_counts() {
         count["A"] = 1; count["C"] = 1
         count["G"] = 1; count["T"] = 1
       }
       BEGIN{ OFS="\t"; reset_counts() }
       {
         curr_locus = $2
         if (prev_locus && (curr_locus != prev_locus)) {
           process_counts(prev_locus)
           reset_counts()
         }
         gsub(/[^ACGT]/,"",$3)
         split($3, a, "")
         for(x in a) count[a[x]] += $1
         prev_locus = curr_locus
       }
       END { process_counts(prev_locus) }' >> \
  ${outdir}/feat_nuc.txt



   
