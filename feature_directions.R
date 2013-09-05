#!/usr/bin/env Rscript
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

# for each (class,feature) compare its mean to the other classes
# in order to color the feature importance heatmap with
# directionality information

args = commandArgs(T)

if (length(args) < 4) {
  cat("USAGE: feature_directions.R x_data y_data params.txt outfile\n")
  q()
}

xfn = args[1]
yfn = args[2]
paramsfn = args[3]
outfn = args[4]

x = read.table(xfn, as.is=T, sep="\t", header=T, row.names=4)
y = read.table(yfn, as.is=T, sep="\t", header=T, row.names=1)[,1]
params = read.table(paramsfn, sep=";", as.is=T)
params = lapply(params[1,], function(s) { unlist(strsplit(s, split="=")) })
names(params) = lapply(params, function (s) { s[1] } )
params = lapply(params, function(s) { s[2] } )

min.reads = as.numeric(params$min_reads)
use.classes = unlist(strsplit(params$use_classes, ","))

# limit loci to >= N reads and belonging to classes of interest
use.loci = x[,4] >= min.reads & y %in% use.classes
x = x[use.loci,]
y = y[use.loci]

x = x[,-(1:5)]

feat.dirs = array(dim=c(length(use.classes), ncol(x) ))
rownames(feat.dirs) = use.classes
colnames(feat.dirs) = colnames(x)
for(cls in use.classes) {
  me = x[y == cls,]
  others = x[y != cls,]
  me.feat.mu = colMeans(me)
  others.feat.mu = colMeans(others)
  feat.dirs[cls,] = sign(me.feat.mu - others.feat.mu)
}

write.table(feat.dirs, col.names=NA, row.names=T, sep="\t", quote=F,
            file=outfn)
