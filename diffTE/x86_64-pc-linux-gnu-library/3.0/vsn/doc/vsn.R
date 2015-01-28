### R code from vignette source 'vsn.Rnw'

###################################################
### code chunk number 1: setup
###################################################
options(width=41, signif=3, digits=3)
set.seed(0xdada)

## To create bitmap versions of plots with many dots, circumventing
##   Sweave's fig=TRUE mechanism...
##   (pdfs are too large)
openBitmap = function(nm, rows=1, cols=1) {
  png(paste("vsn-", nm, ".png", sep=""), 
       width=600*cols, height=700*rows, pointsize=14)
  par(mfrow=c(rows, cols), cex=2)
}


###################################################
### code chunk number 2: overv1
###################################################
require("vsn")
data("kidney")
xnorm = justvsn(kidney)


###################################################
### code chunk number 3: overv2
###################################################
M = exprs(xnorm)[,1] - exprs(xnorm)[,2]


###################################################
### code chunk number 4: overv3
###################################################
fit = vsn2(kidney)
ynorm = predict(fit, kidney)


###################################################
### code chunk number 5: overv4
###################################################
stopifnot(
  identical(exprs(xnorm), exprs(ynorm)), 
  identical(exprs(xnorm), exprs(fit)))


###################################################
### code chunk number 6: nkid-scp1
###################################################
openBitmap("nkid-scp", cols=2)


###################################################
### code chunk number 7: nkid-scp2
###################################################
select = (0==rowSums(exprs(kidney)<=0))
plot(log2(exprs(kidney)[select, ]), 
  main = "a) raw", pch = ".", asp=1)
plot(exprs(xnorm), main = "b) vsn", 
  pch = ".", asp=1,
  col=ifelse(select, "black", "orange"))


###################################################
### code chunk number 8: nkid-scp3
###################################################
dev.off()


###################################################
### code chunk number 9: nkid-sdmean1
###################################################
openBitmap("nkid-sdmean", cols=2)


###################################################
### code chunk number 10: nkid-sdmean2
###################################################
meanSdPlot(xnorm, ranks=TRUE)
meanSdPlot(xnorm, ranks=FALSE)


###################################################
### code chunk number 11: nkid-sdmean3
###################################################
dev.off()


###################################################
### code chunk number 12: nkid-histM
###################################################
hist(M, breaks=33, col="#d95f0e")


###################################################
### code chunk number 13: lymphoma
###################################################
data("lymphoma")
dim(lymphoma)


###################################################
### code chunk number 14: pDatalym
###################################################
pData(lymphoma)


###################################################
### code chunk number 15: lymjustvsn
###################################################
lym = justvsn(lymphoma)


###################################################
### code chunk number 16: lym-sdmean1
###################################################
openBitmap("lym-sdmean")


###################################################
### code chunk number 17: lym-sdmean2
###################################################
meanSdPlot(lym)


###################################################
### code chunk number 18: lym-sdmean3
###################################################
dev.off()


###################################################
### code chunk number 19: lym-M
###################################################
iref = seq(1, 15, by=2)
ismp = seq(2, 16, by=2)
M= exprs(lym)[,ismp]-exprs(lym)[,iref] 
A=(exprs(lym)[,ismp]+exprs(lym)[,iref])/2
colnames(M) = lymphoma$sample[ismp]
colnames(A) = colnames(M)

j = "DLCL-0032"
smoothScatter(A[,j], M[,j], main=j, 
        xlab="A", ylab="M", pch=".")
abline(h=0, col="red")


###################################################
### code chunk number 20: affy1
###################################################
require("affydata")
data("Dilution")
d_vsn = vsnrma(Dilution)


###################################################
### code chunk number 21: affy2
###################################################
d_rma = rma(Dilution)


###################################################
### code chunk number 22: figaffy1
###################################################
openBitmap("affy", cols=3, rows=1)


###################################################
### code chunk number 23: figaffy2
###################################################
par(pch=".")
ax = c(2, 16)

plot(exprs(d_vsn)[,c(1,3)], 
  main = "vsn: array 1 vs 3",
  asp=1, xlim=ax, ylim=ax)

plot(exprs(d_rma)[,c(1,3)], 
  main = "rma: array 1 vs 3",
  asp=1, xlim=ax, ylim=ax)

plot(exprs(d_rma)[,1], 
     exprs(d_vsn)[,1], 
     xlab="rma", ylab="vsn", 
     asp=1, xlim=ax, ylim=ax,
     main = "array 1")
abline(a=0, b=1, col="#ff0000d0")


###################################################
### code chunk number 24: figaffy3
###################################################
dev.off()


###################################################
### code chunk number 25: limma
###################################################
require("limma")
wg = which(lymphoma$dye=="Cy3")
wr = which(lymphoma$dye=="Cy5")

lymRG = new("RGList", list(
  R=exprs(lymphoma)[, wr],
  G=exprs(lymphoma)[, wg]))

lymNCS = justvsn(lymRG)


###################################################
### code chunk number 26: addmeta
###################################################
vmd = data.frame(
  labelDescription=I(c("array ID", 
      "sample in G", "sample in R")),
  channel=c("_ALL", "G", "R"),
  row.names=c("arrayID", "sampG", "sampR"))

arrayID = lymphoma$name[wr]
stopifnot(identical(arrayID, 
          lymphoma$name[wg])) 

## remove sample number suffix 
sampleType = factor(sub("-.*", "", 
                   lymphoma$sample)) 

v = data.frame(
  arrayID = arrayID,
  sampG   = sampleType[wg],
  sampR   = sampleType[wr])
v
              
adf = new("AnnotatedDataFrame", 
  data=v,
  varMetadata=vmd)

phenoData(lymNCS) = adf 


###################################################
### code chunk number 27: width1
###################################################
oo=options(width=75L)


###################################################
### code chunk number 28: lymNCS
###################################################
lymNCS


###################################################
### code chunk number 29: width2
###################################################
options(oo)


###################################################
### code chunk number 30: lymM
###################################################
lymM = (assayData(lymNCS)$R - 
        assayData(lymNCS)$G)


###################################################
### code chunk number 31: design
###################################################
design = model.matrix( ~ lymNCS$sampR)
lf = lmFit(lymM, design[, 2, drop=FALSE])
lf = eBayes(lf)


###################################################
### code chunk number 32: figlimma
###################################################
par(mfrow=c(1,2))
hist(lf$p.value, 100, col="orange") 
pdat=t(lymM[order(lf$p.value)[1:5],])
matplot(pdat, 
  lty=1, type="b", lwd=2, 
  col=hsv(seq(0,1,length=5), 0.7, 0.8), 
  ylab="M", xlab="arrays")


###################################################
### code chunk number 33: makebg
###################################################
rndbg=function(x, off, fac)
 array(off+fac*runif(prod(dim(x))),
       dim=dim(x))

lymRGwbg = lymRG
lymRGwbg$Rb = rndbg(lymRG, 100, 30)
lymRGwbg$Gb = rndbg(lymRG,  50, 20)


###################################################
### code chunk number 34: justvsnwbg
###################################################
lymESwbg = justvsn(lymRGwbg[, 1:3], 
   backgroundsubtract=TRUE)


###################################################
### code chunk number 35: pinId
###################################################
ngr = ngc = 4L
nsr = nsc = 24L
arrayGeometry = data.frame(
  spotcol = rep(1:nsc, 
    times = nsr*ngr*ngc),
  spotrow = rep(1:nsr, 
    each = nsc, times=ngr*ngc),
  pin = rep(1:(ngr*ngc), 
    each = nsr*nsc))


###################################################
### code chunk number 36: strata
###################################################
EconStr = justvsn(lymRG[,1], 
     strata=arrayGeometry$pin)


###################################################
### code chunk number 37: nostrata
###################################################
EsenzaStr = justvsn(lymRG[,1])


###################################################
### code chunk number 38: figstrata1
###################################################
openBitmap("figstrata")


###################################################
### code chunk number 39: figstrata2
###################################################
j = 1L
plot(assayData(EsenzaStr)$R[,j],
     assayData(EconStr)$R[,j],
     pch = ".", asp = 1, 
     col = hsv(seq(0, 1, length=ngr*ngc), 
       0.8, 0.6)[arrayGeometry$pin], 
     xlab = "without strata", 
     ylab = "print-tip strata",
     main = sampleNames(lymNCS)$R[j])


###################################################
### code chunk number 40: figstrata3
###################################################
dev.off()


###################################################
### code chunk number 41: miss1
###################################################
lym2 = lymphoma
nfeat = prod(dim(lym2))
wh = sample(nfeat, nfeat/10)
exprs(lym2)[wh] = NA
table(is.na(exprs(lym2)))


###################################################
### code chunk number 42: miss2
###################################################
fit1 = vsn2(lymphoma, lts.quantile=1)
fit2 = vsn2(lym2, lts.quantile=1)


###################################################
### code chunk number 43: figmiss
###################################################
par(mfrow=c(1,2))
for(j in 1:2){
  p1 = coef(fit1)[,,j]
  p2 = coef(fit2)[,,j]
  d  = max(abs(p1-p2))
  stopifnot(d < c(0.05, 0.03)[j])
  plot(p1, p2, pch = 16, asp = 1,
    main = paste(letters[j], 
      ": max diff=", signif(d,2), sep = ""),
    xlab = "no missing data",
    ylab = "10% of data missing")
  abline(a = 0, b = 1, col = "blue")
}


###################################################
### code chunk number 44: spikein
###################################################
spikeins = 100:200
spfit = vsn2(kidney[spikeins,], 
             lts.quantile=1)
nkid = predict(spfit, newdata=kidney)


###################################################
### code chunk number 45: ref1
###################################################
ref = vsn2(lymphoma[, ismp[1:7]])


###################################################
### code chunk number 46: ref2
###################################################
f8 = vsn2(lymphoma[, ismp[8]], 
          reference = ref)


###################################################
### code chunk number 47: ref3
###################################################
fall = vsn2(lymphoma[, ismp])


###################################################
### code chunk number 48: ref4
###################################################
coefficients(f8)[,1,]
coefficients(fall)[,8,]


###################################################
### code chunk number 49: figref1
###################################################
openBitmap("figref")


###################################################
### code chunk number 50: figref2
###################################################
plot(exprs(f8), exprs(fall)[,8], 
     pch=".", asp=1)
abline(a=0, b=1, col="red")


###################################################
### code chunk number 51: figref3
###################################################
dev.off()


###################################################
### code chunk number 52: hiddenchecks
###################################################
stopifnot(length(ismp)==8L)
maxdiff = max(abs(exprs(f8) - exprs(fall)[,8])) 
if(maxdiff>0.3)
  stop(sprintf("maxdiff is %g", maxdiff))


###################################################
### code chunk number 53: nkid-calib1
###################################################
coef(fit)[1,,]


###################################################
### code chunk number 54: nkid-calib2
###################################################
bkid = kidney
exprs(bkid)[,1]=0.25*(500+exprs(bkid)[,1])


###################################################
### code chunk number 55: nkid-calib3
###################################################
bfit = vsn2(bkid)


###################################################
### code chunk number 56: nkid-calib4
###################################################
oldwarn = options(warn=-1)
openBitmap("nkid-calib", cols=2)


###################################################
### code chunk number 57: nkid-calib5
###################################################
plot(exprs(bkid), main="raw", 
     pch=".", log="xy")
plot(exprs(bfit), main="vsn", 
     pch=".")
coef(bfit)[1,,]


###################################################
### code chunk number 58: nkid-calib6
###################################################
dev.off()
options(oldwarn)
rm(list="oldwarn")


###################################################
### code chunk number 59: vsnQ
###################################################
lym_q = normalizeQuantiles(exprs(lymphoma))
lym_qvsn = vsn2(lym_q, calib="none")


###################################################
### code chunk number 60: figqvsn1
###################################################
openBitmap("qvsn", cols=2, rows=1)


###################################################
### code chunk number 61: figqvsn
###################################################
plot(exprs(lym_qvsn)[, 1:2], pch=".", 
   main="lym_qvsn")
plot(exprs(lym)[,1], exprs(lym_qvsn)[, 1], 
  main="lym_qvsn vs lym", pch=".", 
  ylab="lym_qvsn[,1]", xlab="lym[,1]")


###################################################
### code chunk number 62: figqvsn3
###################################################
dev.off()


###################################################
### code chunk number 63: calcshrink
###################################################
log2.na = function(x){
  w = which(x>0)
  res = rep(as.numeric(NA), length(x))
  res[w] = log2(x[w])
  return(res)
}
  
fc = 0.5                    ## true fold change
x2 = seq(0.5, 15, by=0.5)   ## 'true values' in sample 1
x1 = x2/fc                  ## 'true values' in sample 2
m = s = numeric(length(x1))
n  = 10000
sa = 1
b  = 1
sb = 0.1
sdeta = 0.1
for(i in 1:length(x1)){
  z1 = sa*rnorm(n)+x1[i]*b*exp(sb*rnorm(n))
  z2 = sa*rnorm(n)+x2[i]*b*exp(sb*rnorm(n))
  if(i%%2==1) {
    q = log2.na(z1/z2)
    m[i] = mean(q, na.rm=TRUE)
    s[i] = sd(q, na.rm=TRUE)
  } else {
    h = (asinh(z1/(sa*b/sb))-asinh(z2/(sa*b/sb)))/log(2)
    m[i] = mean(h)
    s[i] = sd(h)
  }
}

colq = c("black", "blue")
lty  = 1
pch  = c(20,18)
cex  = 1.4
lwd  = 2


###################################################
### code chunk number 64: figshrink
###################################################
mai=par("mai"); mai[3]=0; par(mai)
plot(x2, m, ylim=c(-2, 3.5), type="n", xlab=expression(x[2]), ylab="")
abline(h=-log2(fc), col="red", lwd=lwd, lty=1)
abline(h=0, col="black", lwd=1, lty=2)
points(x2, m, pch=pch, col=colq, cex=cex)
segments(x2, m-2*s, x2, m+2*s, lwd=lwd, col=colq, lty=lty)
legend(8.75, -0.1, c("q","h"), col=colq, pch=pch, lty=lty, lwd=lwd)


###################################################
### code chunk number 65: figgraph
###################################################
par(mai = c(0.8, 0.6, 0.01, 0.01))
x = seq(-70.5, 170.5, by=1)
cols = c("black", "blue", "darkgrey")
xoff = cc = 50
ymat = cbind(log2.na(x), log2( (x+sqrt(x^2+cc^2))/2 ), log2.na(x+xoff))
ylim = range(ymat, na.rm=TRUE)
matplot(x, ymat, 
        lty=c(1,1,2), col=cols, lwd=3, type="l", ylim=ylim, ylab="")
abline(v = 0, col = "#80808080", lty = 2)
legend(x = par("usr")[2], y = par("usr")[3], 
       legend = c(expression(y = log[2](x)), expression(y = glog[2](x,c)), expression(y = log[2](x+x[off]))), 
       fill=cols, density = c(NA, NA, 50), xjust=1.1, yjust=-0.1)


###################################################
### code chunk number 66: lymbox
###################################################
colours = hsv(seq(0,1,length=nsr),0.6,1)
j = "CLL-13"
boxplot(A[, j] ~ arrayGeometry$spotrow, 
   col=colours, main=j, 
   ylab="A", xlab="spotrow")


###################################################
### code chunk number 67: lymquscp1
###################################################
openBitmap("lymquscp")


###################################################
### code chunk number 68: lymquscp2
###################################################
plot(A[,j], M[,j], pch=16, cex=0.3,
col=ifelse(arrayGeometry$spotrow%in%(22:23),
           "orange", "black"))
abline(h=0, col="blue")


###################################################
### code chunk number 69: lymquscp3
###################################################
dev.off()


###################################################
### code chunk number 70: sessionInfo
###################################################
toLatex(sessionInfo())


