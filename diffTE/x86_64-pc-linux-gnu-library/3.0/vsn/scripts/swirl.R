library(marrayNorm)
library(vsn)
data(swirl)

if(!exists("nyloess")) {
  ## stratified vsn:
  nystrat = vsn(swirl, strata=as.integer(maPrintTip(swirl)))

  normparams = nystrat@description@preprocessing$vsnParams
  Mvsnstrat = nystrat@exprs[,1:4] - nystrat@exprs[,5:8]

  ## unstratified vsn
  ny = vsn(dat)
  Mvsn = ny@exprs[,1:4] - ny@exprs[,5:8]

  ## print-tip loess
  nyloess = maNorm(swirl, norm="p")
}

nrpt <- max(maPrintTip(swirl))

## boxplots
x11()
par(mfrow=c(4,2))
for (j in 1:4) {
  boxplot(Mvsnstrat[,j] ~maPrintTip(swirl), col=1:nrpt, main = paste("Array", j, "strat"))
  boxplot(Mvsn[,j] ~maPrintTip(swirl), col=1:nrpt, main = paste("Array", j))
}

dev.copy(pdf, file="swirl-boxplots.pdf", width=7, height=12)
dev.off()

x11()
par(mfrow=c(4,2))
## compare stratified and unstratified vsn, and stratified vsn and print-tip loess:
for (j in 1:4) {
  plot(Mvsn[,j], Mvsnstrat[,j], xlab="vsn unstrat", ylab="vsn strat", col=maPrintTip(swirl), pch=".")
  plot(nyloess@maM[,j], Mvsnstrat[,j], xlab="loess strat", ylab="vsn strat", col=maPrintTip(swirl), pch=".")
}
dev.copy(pdf, file="swirl-scpvsn.pdf", width=7, height=12)
dev.off()

#look at the parameters for the 16 print-tips across arrays:
x11()
par(mfrow=c(1,2))
boxplot(as.data.frame(matrix(normparams[1:128], ncol=16)), main="offsets")
boxplot(as.data.frame(matrix(normparams[129:256], ncol=16)), main="factors")
dev.copy(pdf, file="swirl-params.pdf", width=7, height=4)
dev.off()



