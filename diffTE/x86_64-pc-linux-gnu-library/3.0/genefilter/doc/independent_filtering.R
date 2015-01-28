### R code from vignette source 'independent_filtering.Rnw'

###################################################
### code chunk number 1: options
###################################################
options(digits=3, width=100)
library("pasilla") # make sure this is installed, since we need it in the next section


###################################################
### code chunk number 2: libraries
###################################################
library("pasilla")
data("pasillaGenes")


###################################################
### code chunk number 3: DESeq1
###################################################
library("DESeq")


###################################################
### code chunk number 4: DESeq2
###################################################
cds  = estimateSizeFactors( pasillaGenes )
cds  = estimateDispersions( cds )
fit1 = fitNbinomGLMs( cds, count ~ type + condition )
fit0 = fitNbinomGLMs( cds, count ~ type  )


###################################################
### code chunk number 5: DESeq3
###################################################
res = data.frame(
filterstat = rowMeans(counts(cds)),
pvalue    = nbinomGLMTest( fit1, fit0 ),
row.names = featureNames(cds) )


###################################################
### code chunk number 6: headres
###################################################
dim(res)
head(res)


###################################################
### code chunk number 7: pass
###################################################
theta = 0.4
pass = with(res, filterstat > quantile(filterstat, theta))


###################################################
### code chunk number 8: figscatterindepfilt
###################################################
with(res,
  plot(rank(filterstat)/length(filterstat), -log10(pvalue), pch=16, cex=0.45))


###################################################
### code chunk number 9: figecdffilt
###################################################
trsf = function(n) log10(n+1)
plot(ecdf(trsf(res$filterstat)), xlab=body(trsf), main="")


###################################################
### code chunk number 10: badfilter1
###################################################
badfilter = as.numeric(gsub("[+]*FBgn", "", rownames(res)))


###################################################
### code chunk number 11: badfilter2
###################################################
stopifnot(!any(is.na(badfilter)))


###################################################
### code chunk number 12: figbadfilter
###################################################
plot(rank(badfilter)/length(badfilter), -log10(res$pvalue), pch=16, cex=0.45)


###################################################
### code chunk number 13: genefilter
###################################################
library("genefilter")


###################################################
### code chunk number 14: pBH1
###################################################
theta = seq(from=0, to=0.5, by=0.1)
pBH = filtered_p(filter=res$filterstat, test=res$pvalue, theta=theta, method="BH")


###################################################
### code chunk number 15: pBH2
###################################################
head(pBH)


###################################################
### code chunk number 16: figrejection
###################################################
rejection_plot(pBH, at="sample",
               xlim=c(0, 0.5), ylim=c(0, 2000),
               xlab="FDR cutoff (Benjamini & Hochberg adjusted p-value)", main="")


###################################################
### code chunk number 17: filtered_R1
###################################################
theta = seq(from=0, to=0.8, by=0.02)
rejBH = filtered_R(alpha=0.1, filter=res$filterstat, test=res$pvalue, theta=theta, method="BH")


###################################################
### code chunk number 18: fignumreject
###################################################
plot(theta, rejBH, type="l",
     xlab=expression(theta), ylab="number of rejections")


###################################################
### code chunk number 19: differentstats
###################################################
filterChoices = data.frame(
  `mean`   = res$filterstat,
  `geneID` = badfilter,
  `min`    = rowMin(counts(cds)),
  `max`    = rowMax(counts(cds)),
  `sd`     = rowSds(counts(cds))
)
rejChoices = sapply(filterChoices, function(f)
  filtered_R(alpha=0.1, filter=f, test=res$pvalue, theta=theta, method="BH"))


###################################################
### code chunk number 20: colours
###################################################
library("RColorBrewer")
myColours = brewer.pal(ncol(filterChoices), "Set1")


###################################################
### code chunk number 21: figdifferentstats
###################################################
matplot(theta, rejChoices, type="l", lty=1, col=myColours, lwd=2,
        xlab=expression(theta), ylab="number of rejections")
legend("bottomleft", legend=colnames(filterChoices), fill=myColours)


###################################################
### code chunk number 22: histindepfilt
###################################################
h1 = hist(res$pvalue[!pass], breaks=50, plot=FALSE)
h2 = hist(res$pvalue[pass], breaks=50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")


###################################################
### code chunk number 23: fighistindepfilt
###################################################
barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))


###################################################
### code chunk number 24: sortP
###################################################
resFilt = res[pass,]
orderInPlot = order(resFilt$pvalue)
showInPlot = (resFilt$pvalue[orderInPlot] <= 0.06)
alpha = 0.1


###################################################
### code chunk number 25: sortedP
###################################################
plot(seq(along=which(showInPlot)), resFilt$pvalue[orderInPlot][showInPlot],
     pch=".", xlab = expression(rank(p[i])), ylab=expression(p[i]))
abline(a=0, b=alpha/length(resFilt$pvalue), col="red3", lwd=2)


###################################################
### code chunk number 26: doBH
###################################################
whichBH = which(resFilt$pvalue[orderInPlot] <= alpha*seq(along=resFilt$pvalue)/length(resFilt$pvalue))
## Test some assertions:
## - whichBH is a contiguous set of integers from 1 to length(whichBH)
## - the genes selected by this graphical method coincide with those
##   from p.adjust (i.e. padjFilt)
stopifnot(length(whichBH)>0,
          identical(whichBH, seq(along=whichBH)),
          resFilt$FDR[orderInPlot][ whichBH] <= alpha,
          resFilt$FDR[orderInPlot][-whichBH]  > alpha)


###################################################
### code chunk number 27: SchwSpjot
###################################################
j  = round(length(resFilt$pvalue)*c(1, .66))
px = (1-resFilt$pvalue[orderInPlot[j]])
py = ((length(resFilt$pvalue)-1):0)[j]
slope = diff(py)/diff(px)


###################################################
### code chunk number 28: SchwederSpjotvoll
###################################################
plot(1-resFilt$pvalue[orderInPlot],
     (length(resFilt$pvalue)-1):0, pch=".", xaxs="i", yaxs="i",
     xlab=expression(1-p[i]), ylab=expression(N(p[i])))
abline(a=0, slope, col="red3", lwd=2)
abline(h=slope)
text(x=0, y=slope, labels=paste(round(slope)), adj=c(-0.1, 1.3))


###################################################
### code chunk number 29: sessionInfo
###################################################
si = as.character( toLatex( sessionInfo() ) )
cat( si[ -grep( "Locale", si ) ], sep = "\n" )


