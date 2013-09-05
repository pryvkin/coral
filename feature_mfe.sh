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

# compute minimum free energy at a locus
# requires RNAfold and BEDtools

if [ $# -lt 4 ]; then
  echo "USAGE: $0 inbam config_file genome_fas chromInfo" >&2
  exit 1
fi

bam=$1
config=$2
genome_fa=$3
chrominfo=$4

source $config

outdir=`dirname $bam`/coral
mkdir -p $outdir

##
# NOTE: requires sorted indexed bam
echo "Computing MFE..." >&2

echo -e "name\tmfe" > ${outdir}/feat_mfe.txt

# first force bed intervals to fit into chr boundaries
# to placate fastaFromBed
awk 'BEGIN{OFS="\t"; FS="\t"}
     NR==FNR {chrlen[$1] = $2; next}
     $2 < 0 { $2 = 0 }
     $3 > chrlen[$1] { $3 = chrlen[$1] }
     { print }
' $chrominfo ${outdir}/loci.bed | \
while read line; do
  reg=$(echo $line | awk '{print $1":"($2-1)"-"$3}')
  strand=$(echo $line | awk '{print $6}')
  id=$(echo $line | awk '{print $4}')
  echo $id"..." >&2
  echo "$line" | \
    fastaFromBed -s -fi $genome_fa -bed - -tab -fo /dev/stdout | \
    cut -f2
done | RNAfold | \
  awk 'NR % 2 == 0' | \
  sed -e 's/^[^ ]* [(]//; s/[)]$//; s/ //g' | \
  paste <(cut -f4 ${outdir}/loci.bed ) - \
 >> ${outdir}/feat_mfe.txt



   
