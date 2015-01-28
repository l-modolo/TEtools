### R code from vignette source 'DESeq2.Rnw'

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: loadDESeq2
###################################################
library("DESeq2")


###################################################
### code chunk number 3: options
###################################################
options(digits=3, width=80, prompt=" ", continue=" ")


###################################################
### code chunk number 4: quick (eval = FALSE)
###################################################
## dds <- DESeqDataSet(se = se, design = ~ condition)
## dds <- DESeq(dds)
## res <- results(dds)


###################################################
### code chunk number 5: loadSumExp
###################################################
library("parathyroidSE")
data("parathyroidGenesSE")
se <- parathyroidGenesSE
colnames(se) <- colData(se)$run


###################################################
### code chunk number 6: sumExpInput
###################################################
library("DESeq2")
ddsPara <- DESeqDataSet(se = se, design = ~ patient + treatment)
colData(ddsPara)$treatment <- factor(colData(ddsPara)$treatment,
                                     levels=c("Control","DPN","OHT"))
ddsPara


###################################################
### code chunk number 7: loadPasilla
###################################################
library("Biobase")
library("pasilla")
data("pasillaGenes")
countData <- counts(pasillaGenes)
colData <- pData(pasillaGenes)[,c("condition","type")]


###################################################
### code chunk number 8: matrixInput
###################################################
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition)
colData(dds)$condition <- factor(colData(dds)$condition,
                                 levels=c("untreated","treated"))
dds


###################################################
### code chunk number 9: htseqInput
###################################################
library("pasilla")
directory <- system.file("extdata", package="pasilla", mustWork=TRUE)
sampleFiles <- grep("treated",list.files(directory),value=TRUE)
sampleCondition <- sub("(.*treated).*","\\1",sampleFiles)
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)
colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,
                                      levels=c("untreated","treated"))
ddsHTSeq


###################################################
### code chunk number 10: relevel (eval = FALSE)
###################################################
## colData(dds)$condition <- relevel(colData(dds)$condition, "control")


###################################################
### code chunk number 11: deseq
###################################################
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj),]
head(res)


###################################################
### code chunk number 12: deseqNoPrior
###################################################
# this chunk is used to demonstrate the shrinkage
# of log fold changes 
# note: betaPrior can also be used as an argument to DESeq()
ddsNoPrior <- nbinomWaldTest(dds,betaPrior=FALSE)


###################################################
### code chunk number 13: MA
###################################################
plotMA(dds,ylim=c(-2,2),main="DESeq2")


###################################################
### code chunk number 14: MANoPrior
###################################################
plotMA(ddsNoPrior, ylim=c(-2,2),main=expression(unshrunken~log[2]~fold~changes))


###################################################
### code chunk number 15: metadata
###################################################
mcols(res, use.names=TRUE)


###################################################
### code chunk number 16: export (eval = FALSE)
###################################################
## write.csv(as.data.frame(res), 
##           file="condition_treated_results.csv")


###################################################
### code chunk number 17: multifactor
###################################################
colData(dds)


###################################################
### code chunk number 18: replaceDesign
###################################################
design(dds) <- formula(~ type + condition)
dds <- DESeq(dds)


###################################################
### code chunk number 19: multiResults
###################################################
res <- results(dds)
head(res)


###################################################
### code chunk number 20: multiTypeResults
###################################################
resultsNames(dds)
resType <- results(dds, "type_single.read_vs_paired.end")
head(resType)
mcols(resType)


###################################################
### code chunk number 21: defvsd
###################################################
rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)


###################################################
### code chunk number 22: vsd1
###################################################
px     <- counts(dds)[,1] / sizeFactors(dds)[1]
ord    <- order(px)
ord    <- ord[px[ord] < 150]
ord    <- ord[seq(1, length(ord), length=50)]
last   <- ord[length(ord)]
vstcol <- c("blue", "black")
matplot(px[ord],
        cbind(assay(vsd)[, 1], log2(px))[ord, ],
        type="l", lty=1, col=vstcol, xlab="n", ylab="f(n)")
legend("bottomright",
       legend = c(
        expression("variance stabilizing transformation"),
        expression(log[2](n/s[1]))),
       fill=vstcol)


###################################################
### code chunk number 23: vsd2
###################################################
library("vsn")
par(mfrow=c(1,3))
notAllZero <- (rowSums(counts(dds))>0)
meanSdPlot(log2(counts(dds,normalized=TRUE)[notAllZero,] + 1),
           ylim = c(0,2.5))
meanSdPlot(assay(rld[notAllZero,]), ylim = c(0,2.5))
meanSdPlot(assay(vsd[notAllZero,]), ylim = c(0,2.5))


###################################################
### code chunk number 24: heatmap
###################################################
library("RColorBrewer")
library("gplots")
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)


###################################################
### code chunk number 25: figHeatmap2a
###################################################
heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10,6))


###################################################
### code chunk number 26: figHeatmap2b
###################################################
heatmap.2(assay(rld)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))


###################################################
### code chunk number 27: figHeatmap2c
###################################################
heatmap.2(assay(vsd)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))


###################################################
### code chunk number 28: sampleClust
###################################################
distsRL <- dist(t(assay(rld)))


###################################################
### code chunk number 29: figHeatmapSamples
###################################################
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),
                                       paste(condition, type, sep=" : "))
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))


###################################################
### code chunk number 30: figPCA
###################################################
print(plotPCA(rld, intgroup=c("condition", "type")))


###################################################
### code chunk number 31: WaldTest (eval = FALSE)
###################################################
## dds <- estimateSizeFactors(dds)
## dds <- estimateDispersions(dds)
## dds <- nbinomWaldTest(dds)


###################################################
### code chunk number 32: paraSetup
###################################################
ddsCtrst <- ddsPara[, colData(ddsPara)$time == "48h"]
as.data.frame(colData(ddsCtrst)[,c("patient","treatment")])
design(ddsCtrst) <- ~ patient + treatment


###################################################
### code chunk number 33: paraDE
###################################################
ddsCtrst <- DESeq(ddsCtrst)
resultsNames(ddsCtrst)
resPara <- results(ddsCtrst,"treatment_OHT_vs_Control")
head(resPara,2)
mcols(resPara)


###################################################
### code chunk number 34: contrastsChar
###################################################
resCtrst <- results(ddsCtrst, contrast=c("treatment","OHT","DPN"))
head(resCtrst,2)
mcols(resCtrst)


###################################################
### code chunk number 35: contrastsNum
###################################################
resCtrst <- results(ddsCtrst, contrast=c(0,0,0,0,-1,1))
head(resCtrst,2)
mcols(resCtrst)


###################################################
### code chunk number 36: replaceCooksOutlier
###################################################
ddsClean <- replaceOutliersWithTrimmedMean(dds)


###################################################
### code chunk number 37: refitaAfterReplace
###################################################
ddsClean <- DESeq(ddsClean)
tab <- table(initial = results(dds)$padj < .1,
             cleaned = results(ddsClean)$padj < .1)
addmargins(tab)


###################################################
### code chunk number 38: likelihoodRatioTest
###################################################
ddsLRT <- nbinomLRT(dds, reduced = ~ type)
resLRT <- results(ddsLRT)
head(resLRT,2)
mcols(resLRT)
tab <- table(Wald=res$padj < .1, LRT=resLRT$padj < .1)
addmargins(tab)


###################################################
### code chunk number 39: dispFit
###################################################
plotDispEsts(dds)


###################################################
### code chunk number 40: dispFitLocal (eval = FALSE)
###################################################
## ddsLocal <- estimateDispersions(dds, fitType="local")


###################################################
### code chunk number 41: dispFitMean (eval = FALSE)
###################################################
## ddsMean <- estimateDispersions(dds, fitType="mean")


###################################################
### code chunk number 42: dispFitCustom (eval = FALSE)
###################################################
## ddsMed <- estimateDispersionsGeneEst(dds)
## useForMedian <- mcols(ddsMed)$dispGeneEst > 1e-7
## medianDisp <- median(mcols(ddsMed)$dispGeneEst[useForMedian],na.rm=TRUE)
## mcols(ddsMed)$dispFit <- medianDisp
## ddsMed <- estimateDispersionsMAP(ddsMed)


###################################################
### code chunk number 43: filtByMean
###################################################
attr(res,"filterThreshold")
plot(attr(res,"filterNumRej"),type="b",
     ylab="number of rejections")


###################################################
### code chunk number 44: moreFilt
###################################################
resNoFilt <- results(dds, independentFiltering=FALSE)
table(filtering=(res$padj < .1), noFiltering=(resNoFilt$padj < .1))
library(genefilter)
rv <- rowVars(counts(dds,normalized=TRUE))
resFiltByVar <- results(dds, filter=rv)
table(rowMean=(res$padj < .1), rowVar=(resFiltByVar$padj < .1))


###################################################
### code chunk number 45: mcols
###################################################
mcols(dds,use.names=TRUE)[1:4,1:4]
mcols(mcols(dds), use.names=TRUE)[1:4,]


###################################################
### code chunk number 46: normFactors (eval = FALSE)
###################################################
## normFactors <- normFactors / mean(normFactors)
## normalizationFactors(dds) <- normFactors


###################################################
### code chunk number 47: offsetTransform (eval = FALSE)
###################################################
## cqnOffset <- cqnObject$glm.offset
## cqnNormFactors <- exp(cqnOffset)
## EDASeqNormFactors <- exp(-1 * EDASeqOffset)


###################################################
### code chunk number 48: cooksPlot
###################################################
W <- mcols(dds)$WaldStatistic_condition_treated_vs_untreated
maxCooks <- apply(assays(dds)[["cooks"]],1,max)
idx <- !is.na(W)
plot(rank(W[idx]), maxCooks[idx], xlab="rank of Wald statistic", 
     ylab="maximum Cook's distance per gene",
     ylim=c(0,5), cex=.4, col=rgb(0,0,0,.3))
m <- ncol(dds)
p <- 3
abline(h=qf(.99, p, m - p))


###################################################
### code chunk number 49: indFilt
###################################################
plot(res$baseMean+1, -log10(res$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3))


###################################################
### code chunk number 50: histindepfilt
###################################################
use <- res$baseMean > attr(res,"filterThreshold")
table(use)
h1 <- hist(res$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(res$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")


###################################################
### code chunk number 51: fighistindepfilt
###################################################
barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))


###################################################
### code chunk number 52: sortP
###################################################
resFilt <- res[use & !is.na(res$pvalue),]
orderInPlot <- order(resFilt$pvalue)
showInPlot <- (resFilt$pvalue[orderInPlot] <= 0.08)
alpha <- 0.1


###################################################
### code chunk number 53: sortedP
###################################################
plot(seq(along=which(showInPlot)), resFilt$pvalue[orderInPlot][showInPlot],
     pch=".", xlab = expression(rank(p[i])), ylab=expression(p[i]))
abline(a=0, b=alpha/length(resFilt$pvalue), col="red3", lwd=2)


###################################################
### code chunk number 54: doBH
###################################################
whichBH <- which(resFilt$pvalue[orderInPlot] <= alpha*seq(along=resFilt$pvalue)/length(resFilt$pvalue))
## Test some assertions:
## - whichBH is a contiguous set of integers from 1 to length(whichBH)
## - the genes selected by this graphical method coincide with those
##   from p.adjust (i.e. padjFilt)
stopifnot(length(whichBH)>0,
          identical(whichBH, seq(along=whichBH)),
          resFilt$padj[orderInPlot][ whichBH] <= alpha,
          resFilt$padj[orderInPlot][-whichBH]  > alpha)


###################################################
### code chunk number 55: SchwSpjot
###################################################
j  <- round(length(resFilt$pvalue)*c(1, .66))
px <- (1-resFilt$pvalue[orderInPlot[j]])
py <- ((length(resFilt$pvalue)-1):0)[j]
slope <- diff(py)/diff(px)


###################################################
### code chunk number 56: SchwederSpjotvoll
###################################################
plot(1-resFilt$pvalue[orderInPlot],
     (length(resFilt$pvalue)-1):0, pch=".",
     xlab=expression(1-p[i]), ylab=expression(N(p[i])))
abline(a=0, slope, col="red3", lwd=2)


###################################################
### code chunk number 57: sessInfo
###################################################
toLatex(sessionInfo())


###################################################
### code chunk number 58: resetOptions
###################################################
options(prompt="> ", continue="+ ")


