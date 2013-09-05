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

# plot heatmap of feature importance+directionality

library('gplots')

args = commandArgs(T)

if (length(args) < 3) {
  cat("USAGE: plot_feature_heatmap.R feature_importance feature_directions outpdf\n")
  q()
}

feat.imp.infn = args[1]
feat.dir.infn = args[2]
outpdf = args[3]

feat.imp = read.table(feat.imp.infn, header=T, sep="\t", as.is=T, row.names=1)
feat.dir = read.table(feat.dir.infn, header=T, sep="\t", as.is=T, row.names=1)

feat.imp = sweep(feat.imp, 2, colSums(feat.imp), '/')
feat.imp = feat.imp * feat.dir

col=c( rev(rgb(0, seq(0,1,1/127), seq(0,1,1/127))), rgb(seq(0,1,1/127), 0, 0))

pdf(outpdf, points=12)
heatmap.2(t(as.matrix(feat.imp)), trace='none',Rowv=F, Colv='Rowv',dendrogram='none',
          col=col,density.info='none', scale='none',
          margins=c(8, 6))
dev.off()
