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

# split snoRNAs into 3 subclasses based on:
# rmsk IDs  (I use this for human)
# OR
# fRNAdb annotation if given (I use this for mouse)

if [ -z $1 ]; then
    awk 'BEGIN {FS="\t"; OFS="\t"}
     $3 != "snoRNA" {print; next}
     {
       match($6,/Subtype=([^;]+)/,a)
         if (a[1] ~ /^CD/) $3 = "snoRNA_CD"
         if (a[1] ~ /^HAca/) $3 = "snoRNA_HACA"
         if (a[1] ~ /^sca/) $3 = "snoRNA_sca"
       print
     }'
else
    awk 'BEGIN{ FS="\t"; OFS="\t" }
         NR==FNR {annot[$1] = $2; next}
         $3 != "snoRNA" {print; next}
         { 
           match($6, /^ID=(FR.+)/, ma)
           id = ma[1]
           if (!annot[id]) {
             print; next
           } else {
             $3 = annot[id]; print
           }
         }' $1 -
fi


