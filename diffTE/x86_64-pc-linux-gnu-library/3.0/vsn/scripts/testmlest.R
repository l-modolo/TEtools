##------------------------------------------------------
## Test whether the likelihood estimator recovers the
## correct parameters, and plot the likelihood landscape
##--------------------------------------------------

library("vsn")
options(error=recover)
graphics.off()
## set.seed(0xbeeb)

## Generate data
n = 500*(2^(1:9))
dat = sagmbSimulateData(n=n[length(n)], d=1, de=0, nrstrata=1, miss=0, log2scale=TRUE)
ref = new("vsn", mu=dat$mu, sigsq=dat$sigsq, coefficients=dat$coefficients, hoffset=rep(NA_real_, dim(dat$coefficients)[1]))
vin = new("vsnInput", x=dat$y, pstart=array(as.numeric(NA), dim=c(1,1,2)), lts.quantile = 1)

## series of fits with ascending amounts of data
fitpar = sapply(n, function(k) {
  subs = 1:k
  fit = vsn2(dat$y[subs,], reference=ref[subs,], lts.quantile=1, verbose=FALSE)
  coef(fit)
})

graphics.off()

cat("----------------------------------------------------------------------------------\n",
    " Test 1:\n",
    " Do the estimated parameters converge to the true ones when the data is increased?\n",
    "----------------------------------------------------------------------------------\n\n", sep="")
x11(height=6, width=5)
par(mfrow=c(nrow(fitpar), 1))
for(i in 1:nrow(fitpar)) {
  target = if(i<=2) dat$coefficients[1,1,i] else 0
  plot(sqrt(n), fitpar[i,], pch=20, xlab=expression(sqrt(n)), main="Test 1", 
       ylab=letters[i], ylim=range(c(fitpar[i,], target)))
  abline(h=target, col="blue")
}

## pedestrian calculation of logLik
myLLfun = function(object, p, mu, sigsq) {
  res = matrix(as.numeric(NA), nrow=1, ncol=ncol(p))
  stopifnot(nrow(p)==2)
  for(j in 1:ncol(p))  {
    fb = scalingFactorTransformation(p[2,j])
    Y = p[1,j]+(object@x*fb)
    h = asinh(Y)
    scale = nrow(object)/2*log(2*pi*sigsq)
    residuals = sum((h-mu)^2/(2*sigsq))
    jacobi = -nrow(object)*log(fb) + 0.5*sum(log(1+Y^2))
    res[, j] = -(scale+residuals+jacobi)
  }
  return(res)
}

cat("----------------------------------------------------------------------------------\n",
    " Test 2:\n",
    " How does the likelihood landscape look around the true parameters?\n",
    "----------------------------------------------------------------------------------\n\n", sep="")
dat = sagmbSimulateData(n=16000, d=1, de=0, nrstrata=1, miss=0, log2scale=TRUE)
ref = new("vsn", mu=dat$mu, sigsq=dat$sigsq, coefficients=dat$coefficients, hoffset=rep(NA_real_, dim(dat$coefficients)[1]))
vin = new("vsnInput", x=dat$y, pstart=array(as.numeric(NA), dim=c(1,1,2)), lts.quantile = 1)

ex = c(.2, .2)
x11(xpos=400)
lp1 = plotVsnLogLik(vin, dat$coefficients, mu=ref@mu, sigsq=ref@sigsq, expand=ex, main="Test 2")

cat("----------------------------------------------------------------------------------\n",
    " Test 3:\n",
    " Does the likelihood landscape around the true parameters look the same when\n",
    " computed using 'myLLfun'?\n",
    "----------------------------------------------------------------------------------\n\n", sep="")
x11(xpos=800)
lp2 = plotVsnLogLik(vin, dat$coefficients, mu=ref@mu, sigsq=ref@sigsq, expand=ex, fun=myLLfun, main="Test 3")
if(max(abs(lp1$logLik-lp2$logLik))>1e-7) warning("lp1 and lp2 are different!")

wm = which.max(lp1$logLik)
cat(sprintf("True parameters   %11g %11g\n", dat$coefficients[1], dat$coefficients[2]))
cat(sprintf("Maximum of logLik %11g %11g\n\n", lp1[wm,1], lp1[wm,2]))


cat("----------------------------------------------------------------------------------\n",
    " Test 4:\n",
    " What happen if noise is only multiplicative (no profiling)?\n",
    "----------------------------------------------------------------------------------\n\n", sep="")
n = 32000L
mu = runif(n)*2+4
sigsq = 0.01
y  = exp(mu+matrix(rnorm(2L*length(mu), sd=sqrt(sigsq)), ncol=2L))

stopifnot(identical(scalingFactorTransformation(1), exp(1)))

ref = new("vsn", mu=mu, sigsq=sigsq,
          coefficients = array(c(10, -log(2)), dim=c(1,1,2)),
          hoffset = rep(NA_real_, 1))

v1 = vsn2(x=y, ref=ref, lts.quantile=1)
hy = predict(v1, newdata=y)

x11(width=7, height=4); par(mfrow=c(1,2))
plot(log2(y), hy, pch=".", main="Test 4")
abline(a=0, b=1, col="orange")

cat("----------------------------------------------------------------------------------\n",
    " Test 5:\n",
    " What happen if noise is only multiplicative (with profiling)?\n",
    "----------------------------------------------------------------------------------\n\n", sep="")
v2 = vsn2(x=y, lts.quantile=1)

hy = predict(v2, newdata=y)
plot(log2(y), hy, pch=".", main="Test 5")
abline(a=0, b=1, col="orange")
