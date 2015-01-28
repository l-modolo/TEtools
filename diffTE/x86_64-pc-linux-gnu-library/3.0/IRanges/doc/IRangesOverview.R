### R code from vignette source 'IRangesOverview.Rnw'

###################################################
### code chunk number 1: options
###################################################
options(width=72)


###################################################
### code chunk number 2: biocLite (eval = FALSE)
###################################################
## source("http://bioconductor.org/biocLite.R")
## biocLite("IRanges")


###################################################
### code chunk number 3: initialize
###################################################
library(IRanges)


###################################################
### code chunk number 4: initialize
###################################################
set.seed(0)
lambda <- c(rep(0.001, 4500), seq(0.001, 10, length = 500), 
            seq(10, 0.001, length = 500))
xVector <- rpois(1e7, lambda)
yVector <- rpois(1e7, lambda[c(251:length(lambda), 1:250)])


###################################################
### code chunk number 5: basic-ops
###################################################
length(xVector)
xVector[1]
zVector <- c(xVector, yVector)


###################################################
### code chunk number 6: seq-extraction
###################################################
xSnippet <- window(xVector, start = 4751, end = 4760)
xSnippet
head(xSnippet)
tail(xSnippet)
rev(xSnippet)
rep(xSnippet, 2)
window(xSnippet, delta = 2)
subset(xSnippet, xSnippet >= 5L)


###################################################
### code chunk number 7: seq-combine
###################################################
c(xSnippet, rev(xSnippet))
append(xSnippet, xSnippet, after = 3)


###################################################
### code chunk number 8: aggregate
###################################################
xSnippet
aggregate(xSnippet, start = 1:8, width = 3, FUN = median)


###################################################
### code chunk number 9: shiftApply-cor
###################################################
cor(xVector, yVector)
shifts <- seq(235, 265, by=3)
corrs  <- shiftApply(shifts, yVector, xVector, FUN = cor)


###################################################
### code chunk number 10: figshiftcorrs
###################################################
plot(shifts, corrs)


###################################################
### code chunk number 11: Rle-construction
###################################################
xRle <- Rle(xVector)
yRle <- Rle(yVector)
xRle
yRle


###################################################
### code chunk number 12: Rle-vector-compare
###################################################
as.vector(object.size(xRle) / object.size(xVector))
identical(as.vector(xRle), xVector)


###################################################
### code chunk number 13: Rle-accessors
###################################################
head(runValue(xRle))
head(runLength(xRle))


###################################################
### code chunk number 14: Rle-ops
###################################################
xRle > 0
xRle + yRle
xRle > 0 | yRle > 0


###################################################
### code chunk number 15: Rle-summary
###################################################
range(xRle)
sum(xRle > 0 | yRle > 0)


###################################################
### code chunk number 16: Rle-math
###################################################
log1p(xRle)


###################################################
### code chunk number 17: Rle-cor
###################################################
cor(xRle, yRle)
shiftApply(249:251, yRle, xRle, FUN = function(x, y) var(x, y) / (sd(x) * sd(y)))


###################################################
### code chunk number 18: list-intro
###################################################
getClassDef("RleList")


###################################################
### code chunk number 19: list-construct
###################################################
args(IntegerList)
cIntList1 <- IntegerList(x = xVector, y = yVector)
cIntList1
sIntList2 <- IntegerList(x = xVector, y = yVector, compress = FALSE)
sIntList2
## sparse integer list
xExploded <- lapply(xVector[1:5000], function(x) seq_len(x))
cIntList2 <- IntegerList(xExploded)
sIntList2 <- IntegerList(xExploded, compress = FALSE)
object.size(cIntList2)
object.size(sIntList2)


###################################################
### code chunk number 20: list-length
###################################################
length(cIntList2)
Rle(elementLengths(cIntList2))


###################################################
### code chunk number 21: list-lapply
###################################################
system.time(sapply(xExploded, mean))
system.time(sapply(sIntList2, mean))
system.time(sapply(cIntList2, mean))
identical(sapply(xExploded, mean), sapply(sIntList2, mean))
identical(sapply(xExploded, mean), sapply(cIntList2, mean))


###################################################
### code chunk number 22: list-groupgenerics
###################################################
xRleList <- RleList(xRle, 2L * rev(xRle))
yRleList <- RleList(yRle, 2L * rev(yRle))
xRleList > 0
xRleList + yRleList
sum(xRleList > 0 | yRleList > 0)


###################################################
### code chunk number 23: list-endoapply
###################################################
safe.max <- function(x) { if(length(x)) max(x) else integer(0) }
endoapply(sIntList2, safe.max)
endoapply(cIntList2, safe.max)
endoapply(sIntList2, safe.max)[[1]]


###################################################
### code chunk number 24: iranges-constructor
###################################################
ir1 <- IRanges(start = 1:10, width = 10:1)
ir2 <- IRanges(start = 1:10, end = 11)
ir3 <- IRanges(end = 11, width = 10:1)
identical(ir1, ir2) & identical(ir2, ir3)
ir <- IRanges(c(1, 8, 14, 15, 19, 34, 40),
  width = c(12, 6, 6, 15, 6, 2, 7))


###################################################
### code chunk number 25: iranges-start
###################################################
start(ir)


###################################################
### code chunk number 26: iranges-end
###################################################
end(ir)


###################################################
### code chunk number 27: iranges-width
###################################################
width(ir)


###################################################
### code chunk number 28: iranges-subset-numeric
###################################################
ir[1:4]


###################################################
### code chunk number 29: iranges-subset-logical
###################################################
ir[start(ir) <= 15]


###################################################
### code chunk number 30: ranges-extraction
###################################################
ir[[1]]


###################################################
### code chunk number 31: plotRanges
###################################################
plotRanges <- function(x, xlim = x, main = deparse(substitute(x)),
                       col = "black", sep = 0.5, ...) 
{
  height <- 1
  if (is(xlim, "Ranges"))
    xlim <- c(min(start(xlim)), max(end(xlim)))
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  plot.new()
  plot.window(xlim, c(0, max(bins)*(height + sep)))
  ybottom <- bins * (sep + height) - height
  rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col = col, ...)
  title(main)
  axis(1)
}


###################################################
### code chunk number 32: ir-plotRanges
###################################################
plotRanges(ir)


###################################################
### code chunk number 33: ranges-reduce
###################################################
reduce(ir)
plotRanges(reduce(ir))


###################################################
### code chunk number 34: rangeslist-contructor
###################################################
rl <- RangesList(ir, rev(ir))


###################################################
### code chunk number 35: rangeslist-start
###################################################
start(rl)


###################################################
### code chunk number 36: bracket-ranges
###################################################
irextract <- IRanges(start = c(4501, 4901) , width = 100)
xRle[irextract]


###################################################
### code chunk number 37: overlap-ranges
###################################################
ol <- findOverlaps(ir, reduce(ir))
as.matrix(ol)


###################################################
### code chunk number 38: ranges-coverage
###################################################
cov <- coverage(ir)
plotRanges(ir)
cov <- as.vector(cov)
mat <- cbind(seq_along(cov)-0.5, cov)
d <- diff(cov) != 0
mat <- rbind(cbind(mat[d,1]+1, mat[d,2]), mat)
mat <- mat[order(mat[,1]),]
lines(mat, col="red", lwd=4)
axis(2)


###################################################
### code chunk number 39: ranges-shift
###################################################
shift(ir, 10)


###################################################
### code chunk number 40: ranges-plus
###################################################
ir + seq_len(length(ir))


###################################################
### code chunk number 41: ranges-asterisk
###################################################
ir * -2 # half the width


###################################################
### code chunk number 42: ranges-narrow
###################################################
narrow(ir, start=1:5, width=2)


###################################################
### code chunk number 43: ranges-threebands
###################################################
threebands(ir, start=1:5, width=2)


###################################################
### code chunk number 44: ranges-restrict
###################################################
restrict(ir, start=2, end=3)


###################################################
### code chunk number 45: ranges-disjoin
###################################################
disjoin(ir)
plotRanges(disjoin(ir))


###################################################
### code chunk number 46: ranges-disjointBins
###################################################
disjointBins(ir)


###################################################
### code chunk number 47: ranges-reflect
###################################################
reflect(ir, IRanges(start(ir), width=width(ir)*2))


###################################################
### code chunk number 48: ranges-flank
###################################################
flank(ir, width = seq_len(length(ir)))


###################################################
### code chunk number 49: ranges-gaps
###################################################
gaps(ir, start=1, end=50)
plotRanges(gaps(ir, start=1, end=50), c(1,50))


###################################################
### code chunk number 50: ranges-pgap
###################################################



###################################################
### code chunk number 51: ranges-union
###################################################



###################################################
### code chunk number 52: ranges-punion
###################################################



###################################################
### code chunk number 53: ranges-intersect
###################################################



###################################################
### code chunk number 54: ranges-pintersect
###################################################



###################################################
### code chunk number 55: ranges-setdiff
###################################################



###################################################
### code chunk number 56: ranges-psetdiff
###################################################



###################################################
### code chunk number 57: Views-constructors
###################################################
xViews <- Views(xRle, xRle >= 1)
xViews <- slice(xRle, 1)
xViewsList <- slice(xRleList, 1)


###################################################
### code chunk number 58: views-looping
###################################################
head(viewSums(xViews))
viewSums(xViewsList)
head(viewMaxs(xViews))
viewMaxs(xViewsList)


###################################################
### code chunk number 59: RangedData-construct
###################################################
values <- rnorm(length(ir))
rd <- RangedData(ir, name = letters[seq_len(length(ir))], values)
rd


###################################################
### code chunk number 60: RangedData-construct-space
###################################################
rd <- RangedData(ir, name = letters[seq_len(length(ir))], values, 
                 space = rep(c("chr1", "chr2"), c(3, length(ir) - 3)))
rd


###################################################
### code chunk number 61: RangedData-ranges
###################################################
ranges(rd)


###################################################
### code chunk number 62: RangedData-values
###################################################
values(rd)


###################################################
### code chunk number 63: RangedData-subspace
###################################################
rd["chr1"]


###################################################
### code chunk number 64: RangedData-subspace-2
###################################################
all(identical(rd["chr1"], rd[1]),
    identical(rd[1], rd[c(TRUE, FALSE)]))


###################################################
### code chunk number 65: RangedData-names
###################################################
names(rd)


###################################################
### code chunk number 66: RangedData-length
###################################################
length(rd)


###################################################
### code chunk number 67: RangedData-lapply
###################################################
lapply(rd, names)


###################################################
### code chunk number 68: RangedData-extract
###################################################
rd[[2]]


###################################################
### code chunk number 69: RangedData-dollar
###################################################
rd$values


###################################################
### code chunk number 70: RangedData-subset-2d
###################################################
rd[1:3, "name"]


###################################################
### code chunk number 71: sessionInfo
###################################################
toLatex(sessionInfo())


