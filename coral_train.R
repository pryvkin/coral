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
library('varSelRF')
library('randomForest')
library('caret')
library('foreach')
library('doParallel')
library('digest')

option_list = list(
  make_option(c("-r", "--min-read-count"), type="integer", default=20,
              help="Minimum # of reads at a locus (20)", dest='min.reads'),
  make_option(c("-c", "--classes"), type="character",
              default="lincRNA_exon,miRNA,scRNA,snRNA,snoRNA_CD,snoRNA_HACA,transposon",
              help="Comma-delimited list of classes to use", dest='use.classes.str'),
  make_option(c("-p", "--cores"), type="integer", default=1,
              help="Number of processes to run in parallel (1)", dest='ncpu'),
  make_option(c("-v", "--varsel-runs"), type="integer", default=100,
              help="Number of times to run variable selection (100)", dest='n.varsel'),
  make_option(c("-n", "--rf-runs"), type="integer", default=100,
              help="Number of times to run random forest (100)", dest='n.rf'),
  make_option(c("-t", "--trees"), type="integer", default=1000,
              help="Number of trees in forest (1000)", dest='n.trees'),
  make_option(c("-e", "--no-feature-selection"), action="store_true", default=FALSE,
              help="Skip feature importance computation (off)", dest='no.feat.sel'),
  make_option(c("-f", "--feature-selection-only"), action="store_true", default=FALSE,
              help="Skip feature importance computation (off)", dest='only.feat.sel')
  )
parser = OptionParser(usage = "%prog [options] data_x data_y outdir", option_list=option_list)
arguments = parse_args(parser, positional_arguments = TRUE)
opt = arguments$options

if(length(arguments$args) != 3) {
  print_help(parser)
  q(status=1)
} else {
  xfn = arguments$args[1]
  yfn = arguments$args[2]
  outdir = arguments$args[3]
}

use.classes = sort(unlist(strsplit(opt$use.classes.str,",")))

param.str=sprintf("min_reads=%d;use_classes=%s;n.varsel=%d;n.rf=%d;trees=%d",
  opt$min.reads, paste(use.classes,collapse=','),
  opt$n.varsel, opt$n.rf, opt$n.trees)
param.hash = digest(param.str, 'md5')
outdir = sprintf("%s/run_%s", outdir, param.hash)
dir.create(outdir)

cat(sprintf("%s\n", param.str),
    file=sprintf("%s/params.txt", outdir))

x = read.table(xfn, sep="\t", as.is=T, header=T, row.names=4)
y = read.table(yfn, sep="\t", as.is=T, header=T, row.names=1)


### subset data
use.loci = x$reads >= opt$min.reads & y[,1] %in% use.classes

x = x[use.loci,]
y = y[use.loci,]

x = x[,-(1:5)]
y = as.factor(y)

clu = makeCluster(opt$ncpu)
registerDoParallel(clu)

### variable importance
# generate OAA (one-against-all) subsets of labels
if (!opt$no.feat.sel) {
  y.oaa = list()
  varsel.oaa = list()
  for(cls in use.classes) {
    y.oaa[[cls]] = as.character(y)
    y.oaa[[cls]][y.oaa[[cls]] != cls] = 0
    y.oaa[[cls]][y.oaa[[cls]] == cls] = 1
    y.oaa[[cls]] = as.factor(y.oaa[[cls]])
    # run varSelRF in parallel
    varsel.oaa[[cls]] =
      unlist( foreach(i=1:(opt$n.varsel), .packages=c('varSelRF')) %dopar% {
        varSelRF(x, y.oaa[[cls]], mtryFactor=4,
                 ntree = opt$n.trees, vars.drop.frac = 0.35)$selected.vars
      } )
  }
  # summarize varsel data
  varsel.sum = array(0, dim=c(length(use.classes), ncol(x)))
  rownames(varsel.sum) = use.classes
  colnames(varsel.sum) = colnames(x)
  for(cls in use.classes) {
    cls.vartab = table(varsel.oaa[[cls]])
    for(i in 1:length(cls.vartab))
      varsel.sum[cls, names(cls.vartab)[i]] = cls.vartab[i]
  }
  write.table(varsel.sum, col.names=NA, row.names=T, sep="\t", quote=F,
              file=sprintf("%s/feature_importance.txt", outdir))
}

### build RF model
rf.models = foreach(i=1:(opt$n.rf), .packages=c('randomForest')) %dopar% {
  randomForest(y~., data = data.frame(x,y),
               keep.forest=TRUE, proximity=T, ntrees=opt$n.trees)
}
if (opt$only.feat.sel) {
  q()
}

### assess performance by averaging across all runs
perf.measures = c('Sensitivity', 'Pos Pred Value', 'Specificity', 'Prevalence')
perf.measure.labels = c('sensitivity', 'ppv', 'specificity', 'prevalence')
class.perf = array(NA, dim=c(length(use.classes), 2*length(perf.measures)))
rownames(class.perf) = use.classes
colnames(class.perf) = as.vector(sapply(perf.measure.labels, paste, c('avg', 'sd'), sep='_'))

for(cls in use.classes) {
  cls.allperf = data.frame(lapply(rf.models, function(rf.m) {
    confusionMatrix(rf.m$predicted, y)$byClass[sprintf("Class: %s", cls),
                                               perf.measures]
  } ))
  class.perf[cls, paste(perf.measure.labels, "_avg", sep='')] =
             rowMeans(cls.allperf, na.rm=T)
  class.perf[cls, paste(perf.measure.labels, "_sd", sep='')] =
    apply(cls.allperf, 1, sd, na.rm=T)
}
write.table(class.perf, col.names=NA, row.names=T, sep="\t", quote=F,
            file=sprintf("%s/class_performance.txt", outdir))

### overall accuracy
overall = data.frame(lapply(rf.models, function(rf.m) {
  confusionMatrix(rf.m$predicted,y)$overall } ))
overall.avg = rowMeans(overall)
overall.sd = apply(overall, 1, sd)
write.table(cbind(avg=overall.avg, sd=overall.sd),
            col.names=NA, row.names=T, sep="\t",quote=F,
            file=sprintf("%s/overall_accuracy.txt", outdir))

### output class sizes
ytab = table(y)
cls.sizes = rep(0, length(use.classes))
names(cls.sizes) = use.classes
for(i in 1:length(ytab))
  cls.sizes[names(ytab)[i]] = ytab[i]
write.table(cbind(n=cls.sizes), col.names=NA, row.names=T, sep="\t", quote=F,
            file=sprintf("%s/class_sizes.txt", outdir))

# save workspace
save(list=ls(all=T), file=sprintf("%s/data.Rdata", outdir))

