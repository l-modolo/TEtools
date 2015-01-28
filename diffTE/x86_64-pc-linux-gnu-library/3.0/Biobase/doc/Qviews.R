### R code from vignette source 'Qviews.Rnw'

###################################################
### code chunk number 1: doa
###################################################
if (!("Biobase" %in% search())) library(Biobase)
if (!("ALL" %in% search())) library(ALL)
if (!("ALL" %in% objects())) data(ALL)


###################################################
### code chunk number 2: getaEF (eval = FALSE)
###################################################
## library(Biobase)
## library(ALL)
## data(ALL)
## ALL


###################################################
### code chunk number 3: mkprov
###################################################
dataSource = function(dsn) {
 if (!is(dsn, "character")) dsn = try(deparse(substitute(dsn)))
 if (inherits(dsn, "try-error")) stop("can't parse dsn arg")
 dd = data()$results
 if (is.na(match(dsn, dd[,"Item"]))) return(NULL)
 paste("package:", dd[ dd[,"Item"] == dsn, "Package" ], sep="")
}
 



###################################################
### code chunk number 4: newmeth
###################################################
setGeneric("peek", function(x,maxattr=10)standardGeneric("peek"))
setMethod("peek", c("eSet", "numeric"), function(x,maxattr=10) {
 ds = dataSource(deparse(substitute(x)))
 if (!is.null(ds)) ds = paste(" [from ", ds, "]", sep="")
  else ds = ""
 cat(deparse(substitute(x)), ds, ":\n", sep="")
 cat("Platform annotation: ", annotation(x),"\n")
 cat("primary assay results are:\n")
 print(dim(x))
 cat("sample attributes are:\n")
   vn = rownames(varMetadata(x))
   ld = substr(varMetadata(x)$labelDescription,1,50)
   dd = data.frame("labelDescription[truncated]"=ld)
   rownames(dd) = vn
 if ((ndd <- nrow(dd)) <= maxattr) show(dd)
 else {
    cat("first", maxattr, "of", ndd, "attributes:\n")
    show(dd[1:maxattr,,drop=FALSE])
    }
 cat("----------\n")
 cat("use varTable to see values/freqs of all sample attributes\n")
 cat("----------\n")
})
setMethod("peek", c("eSet", "missing"), function(x,maxattr=10) {
 ds = dataSource(deparse(substitute(x)))
 if (!is.null(ds)) ds = paste(" [from ", ds, "]", sep="")
  else ds = ""
 cat(deparse(substitute(x)), ds, ":\n", sep="")
 cat("Platform annotation: ", annotation(x),"\n")
 cat("primary assay results are:\n")
 print(dim(x))
 cat("sample attributes are:\n")
   vn = rownames(varMetadata(x))
   ld = substr(varMetadata(x)$labelDescription,1,50)
   dd = data.frame("labelDescription[truncated]"=ld)
   rownames(dd) = vn
 if ((ndd <- nrow(dd)) <= maxattr) show(dd)
 else {
    cat("first", maxattr, "of", ndd, "attributes:\n")
    show(dd[1:maxattr,,drop=FALSE])
    }
 cat("----------\n")
 cat("use varTable to see values/freqs of all sample attributes\n")
 cat("----------\n")
})
setGeneric("varTable", function(x, full=FALSE, max=Inf) standardGeneric("varTable"))
setMethod("varTable", c("eSet", "missing", "ANY"), function(x, full=FALSE, max=Inf) varTable(x, FALSE, max))
setMethod("varTable", c("eSet", "logical", "ANY"), function(x, full=FALSE, max=Inf) {
   ans = lapply( names(pData(x)), function(z)table(x[[z]]) )
   tans = lapply(ans, names)
   kp = 1:min(max,length(tans))
   if (!full) ans = sapply(tans, selectSome, 3)[kp]
   else ans = tans[kp]
   names(ans) = names(pData(x))[kp]
   ans
})
setGeneric("varNames", function(x) standardGeneric("varNames"))
setMethod("varNames", "eSet", function(x) names(pData(x)))


###################################################
### code chunk number 5: lka
###################################################
peek(ALL)


###################################################
### code chunk number 6: lkv
###################################################
varNames(ALL)


###################################################
### code chunk number 7: lkvn
###################################################
varTable(ALL, max=4)


###################################################
### code chunk number 8: lkvn
###################################################
varTable(ALL, full=TRUE, max=4)


