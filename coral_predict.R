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

library('optparse')
library('randomForest')
library('caret')

option_list = list(
  make_option(c("-r", "--min-read-count"), type="integer", default=20,
              help="Minimum # of reads at a locus (20)", dest='min.reads'),
  make_option(c("-y", "--labels"), type="character", default='',
              help="File with labels, used to assess performance", dest='yfn')
  )
parser = OptionParser(usage = "%prog [options] data_x modeldir outdir", option_list=option_list)
arguments = parse_args(parser, positional_arguments = TRUE)
opt = arguments$options

if(length(arguments$args) != 3) {
  print_help(parser)
  q(status=1)
} else {
  xfn = arguments$args[1]
  modeldir = arguments$args[2]
  outdir = arguments$args[3]
}

dir.create(outdir)

my.opt = opt
my.outdir = outdir
my.yfn = opt$yfn
my.xfn = xfn
load(sprintf("%s/data.Rdata", modeldir))
outdir = my.outdir
y = NULL
yfn = my.yfn
opt = my.opt
xfn = my.xfn

x = read.table(xfn, sep="\t", as.is=T, header=T, row.names=4)
have.labels = (yfn != '')

### subset data
use.loci = x$reads >= opt$min.reads

x = x[use.loci,]
orig.x = x
x = x[,-(1:5)]

preds = data.frame(lapply(rf.models, function(rf.m) {
  predict(rf.m, x)} ) )

preds = apply(preds, 1, function(r) { names(sort(-table(r)))[1] })
preds = factor(preds, levels=use.classes)

write.table(cbind(orig.x[,1:3], name=rownames(orig.x), orig.x[,4:5], label=preds),
            col.names=T, row.names=F, sep="\t", quote=F,
            file=sprintf("%s/pred.txt", outdir) )

if (have.labels) {
  y = read.table(yfn, sep="\t", as.is=T, header=T, row.names=1)

  y = factor(y[use.loci,1], levels = use.classes)
  use.y = as.character(y) %in% use.classes

  cm = confusionMatrix(preds[use.y], y[use.y])
  
  ### assess performance
  perf.measures = c('Sensitivity', 'Pos Pred Value', 'Specificity', 'Prevalence')
  perf.measure.labels = c('sensitivity', 'ppv', 'specificity', 'prevalence')
  class.perf = array(NA, dim=c(length(use.classes), length(perf.measures)))
  rownames(class.perf) = use.classes
  colnames(class.perf) = perf.measure.labels

  for(cls in use.classes) {
    class.perf[cls,] = cm$byClass[sprintf("Class: %s", cls),
      perf.measures]
  }

  write.table(class.perf, col.names=NA, row.names=T, sep="\t", quote=F,
              file=sprintf("%s/pred_class_performance.txt", outdir))
  
  ### overall accuracy
  write.table(cm$overall,
              col.names=NA, row.names=T, sep="\t",quote=F,
              file=sprintf("%s/pred_overall_accuracy.txt", outdir))
  
  ### output class sizes
  ytab = table(y)
  cls.sizes = rep(0, length(use.classes))
  names(cls.sizes) = use.classes
  for(i in 1:length(ytab))
  cls.sizes[names(ytab)[i]] = ytab[i]
  write.table(cbind(n=cls.sizes), col.names=NA, row.names=T, sep="\t", quote=F,
              file=sprintf("%s/pred_class_sizes.txt", outdir))
  
}

