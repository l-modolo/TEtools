### R code from vignette source 'howtogenefilter.Rnw'

###################################################
### code chunk number 1: howtogenefilter.Rnw:41-47
###################################################
library("Biobase")
library("genefilter")
data(sample.ExpressionSet)
varLabels(sample.ExpressionSet)
table(sample.ExpressionSet$sex)
table(sample.ExpressionSet$type)


###################################################
### code chunk number 2: howtogenefilter.Rnw:70-74
###################################################
f1 <- kOverA(5, 200)
ffun <- filterfun(f1)
wh1 <- genefilter(exprs(sample.ExpressionSet), ffun)
sum(wh1)


###################################################
### code chunk number 3: howtogenefilter.Rnw:85-88
###################################################
f2 <- ttest(sample.ExpressionSet$type, p=0.1)
wh2 <- genefilter(exprs(sample.ExpressionSet), filterfun(f2))
sum(wh2)


###################################################
### code chunk number 4: howtogenefilter.Rnw:100-103
###################################################
ffun_combined <- filterfun(f1, f2)
wh3 <- genefilter(exprs(sample.ExpressionSet), ffun_combined)
sum(wh3)


###################################################
### code chunk number 5: aggregate
###################################################

 knnCV <- function(EXPR, selectfun, cov, Agg, pselect = 0.01, Scale=FALSE) {
   nc <- ncol(EXPR)
   outvals <- rep(NA, nc)
   for(i in 1:nc) {
      v1 <- EXPR[,i]
      expr <- EXPR[,-i]
      glist <- selectfun(expr, cov[-i], p=pselect)
      expr <- expr[glist,]
      if( Scale ) {
        expr <- scale(expr)
        v1 <- as.vector(scale(v1[glist]))
      }
      else
         v1 <- v1[glist]
      out <- paste("iter ",i, " num genes= ", sum(glist), sep="")
      print(out)
      Aggregate(row.names(expr), Agg)
      if( length(v1) == 1)
         outvals[i] <- knn(expr, v1, cov[-i], k=5)
      else
          outvals[i] <- knn(t(expr), v1, cov[-i], k=5)
    }
    return(outvals)
  }


###################################################
### code chunk number 6: aggregate
###################################################
 gfun <- function(expr, cov, p=0.05) {
    f2 <- ttest(cov, p=p)
    ffun <- filterfun(f2)
    which <- genefilter(expr, ffun)
  }



###################################################
### code chunk number 7: aggregate
###################################################
  library("class")

  ##scale the genes
  ##genescale is a slightly more flexible "scale"
  ##work on a subset -- for speed only
  geneData <- genescale(exprs(sample.ExpressionSet)[1:75,], 1)

  Agg <- new("aggregator")

  testcase <- knnCV(geneData, gfun, sample.ExpressionSet$type, 
         Agg, pselect=0.05)


###################################################
### code chunk number 8: aggregate
###################################################
sort(sapply(aggenv(Agg), c), decreasing=TRUE)


###################################################
### code chunk number 9: howtogenefilter.Rnw:207-208
###################################################
toLatex(sessionInfo())


