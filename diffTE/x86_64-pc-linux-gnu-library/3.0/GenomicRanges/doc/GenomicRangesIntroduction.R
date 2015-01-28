### R code from vignette source 'GenomicRangesIntroduction.Rnw'

###################################################
### code chunk number 1: options
###################################################
options(width=72)


###################################################
### code chunk number 2: biocLite (eval = FALSE)
###################################################
## source("http://bioconductor.org/biocLite.R")
## biocLite("GenomicRanges")


###################################################
### code chunk number 3: initialize
###################################################
library(GenomicRanges)


###################################################
### code chunk number 4: example-GRanges
###################################################
gr <-
  GRanges(seqnames =
          Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
          ranges =
          IRanges(1:10, end = 7:16, names = head(letters, 10)),
          strand =
          Rle(strand(c("-", "+", "*", "+", "-")),
              c(1, 2, 2, 3, 2)),
          score = 1:10,
          GC = seq(1, 0, length=10))
gr


###################################################
### code chunk number 5: GRanges-location-accessors
###################################################
seqnames(gr)
ranges(gr)
strand(gr)


###################################################
### code chunk number 6: metadataAccess
###################################################
mcols(gr)
mcols(gr)$score


###################################################
### code chunk number 7: setSeqLengths
###################################################
seqlengths(gr) <- c(249250621,243199373,198022430)


###################################################
### code chunk number 8: setSeqLengths2
###################################################
seqlengths(gr)


###################################################
### code chunk number 9: names
###################################################
names(gr)
length(gr)


###################################################
### code chunk number 10: splitAppendGRanges
###################################################
sp <- split(gr, rep(1:2, each=5))
sp


###################################################
### code chunk number 11: combine
###################################################
c(sp[[1]], sp[[2]])


###################################################
### code chunk number 12: subset1
###################################################
gr[2:3]


###################################################
### code chunk number 13: subset2
###################################################
gr[2:3, "GC"]


###################################################
### code chunk number 14: assign1
###################################################
singles <- split(gr, names(gr))
grMod <- gr
grMod[2] <- singles[[1]]
head(grMod, n=3)


###################################################
### code chunk number 15: assign2
###################################################
grMod[2,1] <- singles[[3]][,1]
head(grMod, n=3)


###################################################
### code chunk number 16: other
###################################################
rep(singles[[2]], times = 3)
rev(gr)
head(gr,n=2)
tail(gr,n=2)
window(gr, start=2,end=4)
gr[IRanges(start=c(2,7), end=c(3,9))]


###################################################
### code chunk number 17: IRangesStuff
###################################################
g <- gr[1:3]
g <- append(g, singles[[10]])
start(g)
end(g)
width(g)
range(g)


###################################################
### code chunk number 18: flank
###################################################
flank(g, 10)


###################################################
### code chunk number 19: flank2
###################################################
flank(g, 10, start=FALSE)


###################################################
### code chunk number 20: shiftAndResize
###################################################
shift(g, 5)
resize(g, 30)


###################################################
### code chunk number 21: reduce
###################################################
reduce(g)


###################################################
### code chunk number 22: gaps
###################################################
gaps(g)


###################################################
### code chunk number 23: disjoin
###################################################
disjoin(g)


###################################################
### code chunk number 24: coverage
###################################################
coverage(g)


###################################################
### code chunk number 25: intervals1
###################################################
g2 <- head(gr, n=2)
union(g, g2)
intersect(g, g2)
setdiff(g, g2)


###################################################
### code chunk number 26: intervals2
###################################################
g3 <- g[1:2]
ranges(g3[1]) <- IRanges(start=5, end=12)
punion(g2, g3)
pintersect(g2, g3)
psetdiff(g2, g3)


###################################################
### code chunk number 27: manPage (eval = FALSE)
###################################################
## ?GRanges


###################################################
### code chunk number 28: example-GRangesList
###################################################
gr1 <-
  GRanges(seqnames = "chr2", ranges = IRanges(3, 6),
          strand = "+", score = 5L, GC = 0.45)
gr2 <-
  GRanges(seqnames = c("chr1", "chr1"),
          ranges = IRanges(c(7,13), width = 3),
          strand = c("+", "-"), score = 3:4, GC = c(0.3, 0.5))
grl <- GRangesList("txA" = gr1, "txB" = gr2)
grl


###################################################
### code chunk number 29: basicGRLAccessors
###################################################
seqnames(grl)
ranges(grl)
strand(grl)


###################################################
### code chunk number 30: exceptions
###################################################
length(grl)
names(grl)
seqlengths(grl)


###################################################
### code chunk number 31: elementLengths
###################################################
elementLengths(grl)


###################################################
### code chunk number 32: isEmpty
###################################################
isEmpty(grl)


###################################################
### code chunk number 33: mcolsGRL
###################################################
mcols(grl) <- c("Transcript A","Transcript B")
mcols(grl)


###################################################
### code chunk number 34: unlistGRL
###################################################
ul <- unlist(grl)


###################################################
### code chunk number 35: intOpsGRL
###################################################
start(grl)
end(grl)
width(grl)


###################################################
### code chunk number 36: coverageGRL
###################################################
shift(grl, 20)
coverage(grl)


###################################################
### code chunk number 37: subsetGRL (eval = FALSE)
###################################################
## grl[1]
## grl[[1]]
## grl["txA"]
## grl$txB


###################################################
### code chunk number 38: subsetGRL2
###################################################
grl[1, "score"]
grl["txB", "GC"]


###################################################
### code chunk number 39: otherSubsetGRL
###################################################
rep(grl[[1]], times = 3)
rev(grl)
head(grl, n=1)
tail(grl, n=1)
window(grl, start=1, end=1)
grl[IRanges(start=2, end=2)]


###################################################
### code chunk number 40: lapply
###################################################
lapply(grl, length)
sapply(grl, length)


###################################################
### code chunk number 41: mapply
###################################################
grl2 <- shift(grl, 10)
names(grl2) <- c("shiftTxA", "shiftTxB")

mapply(c, grl, grl2)
Map(c, grl, grl2)


###################################################
### code chunk number 42: endoapply
###################################################
endoapply(grl,rev)


###################################################
### code chunk number 43: mendoapply
###################################################
mendoapply(c,grl,grl2)


###################################################
### code chunk number 44: ReduceGRL
###################################################
Reduce(c,grl)


###################################################
### code chunk number 45: manPage2 (eval = FALSE)
###################################################
## ?GRangesList


###################################################
### code chunk number 46: findOverlaps
###################################################
mtch <- findOverlaps(gr, grl)
as.matrix(mtch)


###################################################
### code chunk number 47: countOL
###################################################
countOverlaps(gr, grl)


###################################################
### code chunk number 48: subsetByOverlaps
###################################################
subsetByOverlaps(gr,grl)


###################################################
### code chunk number 49: select-first
###################################################
findOverlaps(gr, grl, select="first")
findOverlaps(grl, gr, select="first")


###################################################
### code chunk number 50: readGAlignments
###################################################
library(Rsamtools)
aln1_file <- system.file("extdata", "ex1.bam", package="Rsamtools")
aln1 <- readGAlignments(aln1_file)
aln1
length(aln1)


###################################################
### code chunk number 51: accessors
###################################################
head(seqnames(aln1))
seqlevels(aln1)
head(strand(aln1))
head(cigar(aln1))
head(qwidth(aln1))
head(start(aln1))
head(end(aln1))
head(width(aln1))
head(ngap(aln1))


