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

# combine loci and features into one file

if [ $# -lt 1 ]; then
  echo "USAGE: $0 data_dir" >&2
  exit 1
fi

indir=$1

echo -e "chr\tstart\tend\tname\treads\tstrand" | \
  cat - $indir/loci.bed | \
  paste - <(cut -f2- $indir/feat_antisense.txt) \
          <(cut -f2- $indir/feat_posentropy.txt) \
          <(cut -f2- $indir/feat_nuc.txt) \
          <(cut -f2- $indir/feat_mfe.txt) \
          <(cut -f2- $indir/feat_lengths.txt) \
  > $indir/data_x.txt

if [ -e $indir/loci.annot ]; then
  echo -e "name\tannot" > $indir/data_y.txt
  cut -f1,3 $indir/loci.annot >> $indir/data_y.txt
fi

