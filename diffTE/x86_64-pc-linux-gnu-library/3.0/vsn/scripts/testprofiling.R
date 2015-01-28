##
## Test concordance of maxima of profile likelihood and normal likelihood
##

library("vsn")

n = 1000L
d = 2L
nrstr = 1L
nrpar = 2L*d*nrstr

dat = sagmbSimulateData(n=n, d=d, de=0, nrstrata=nrstr, miss=0, log2scale=TRUE)
fit = vsn2(dat$y, lts.quantile=1)

stopifnot(all(diff(dat$strata)>=0)) ## >1 strata not yet implemented here
v = new("vsnInput", x=dat$y, pstart=array(as.numeric(NA), dim=c(nrstr, d, 2)),
  strata=factor(dat$strata), ordered=TRUE)

nplot = 31
par(mfcol=c(2, nrpar))

for(i in seq_len(nrpar))  {
  pars = matrix(coef(fit), nrow=nrpar, ncol=nplot)
  pars[i,] = coef(fit)[i] +  seq(-1, 1, length=nplot)
  logarg  = ""
  if(i<=d*nrstr) {
    xlab = substitute(a[k], list(k=i))
  } else {
    xlab = substitute(b[k], list(k=i-nrpar/2))
  }
  
  for(what in 1:2){
    ll = switch(what,
      logLik(v, pars),                             ## without reference
      logLik(v, pars, mu=fit@mu, sigsq=fit@sigsq)) ## with reference

    plot(pars[i,], ll[1, ], type="b",
         xlab = xlab, log=logarg, main=round(ll[1, (ncol(ll)+1)/2], 1),
         ylab = expression(-log(L)))
    abline(v=coefficients(fit)[i], col="red")
  } ## for what
}

cat("PLL=", logLik(v, matrix(coef(fit), ncol=1)),
   "\n LL=", logLik(v, matrix(coef(fit), ncol=1),  mu=fit@mu, sigsq=fit@sigsq), "\n")




