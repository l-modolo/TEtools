## A quick demo to overview the package
##         -- Laurent

if(dev.cur() <= 1) get(getOption("device"))()

opar <- par(ask= (interactive() &&
                  (.Device %in% c("X11", "GTK", "windows", "Macintosh"))),
            bg="cornsilk", mfrow=c(1,1))



## load the data
library(affydata)
data(cdf.example)
data(Dilution)

## display the image of the data in the CEL file

cel <- Dilution[[1]]

image(cel,transfo=log)

image(cel,transfo=log)

## find the locations for probes corresponding to a given ID

l.pm <- locate.name("AFFX-BioC-5_at", cdf.example, type="pm")
plotLocation(l.pm, cdf.example, col="red")
l.mm <- locate.name("AFFX-BioC-5_at", cdf.example, type="mm")
plotLocation(l.mm, cdf.example, col="blue")

arrows(20,60, min(l.pm[,1]), min(l.pm[,2]), col="white")
legend(20,60,c("perfect match","mismatch"),c("red","blue"),bg="white")

rm(l.pm, l.mm)



## The '5', 'M' or '3' in the names means the set of probes relates to
## the 5-prime, middle or 3-prime sectors  respectively.
## The at or st ending means the sequence relates to the complementary
## sequence of the gene or not respectively.

namesspot <- c("AFFX-BioB-5_at","AFFX-BioB-M_at", "AFFX-BioB-3_at")

p <- probeset(Dilution, genenames=namesspot)

par(mfrow=c(3,3))

for (pps in p)
  barplot(pps)


## normalize
plot.new()
par(mfrow=c(2,2))

nat <- pmormm(cdf.example)

cel2 <- Dilution[[2]]
plot(intensity(cel), intensity(cel2), xlab="CEL file 1", ylab="CEL file 2",main="raw values",sub="all probes plotted",type="n")
points(intensity(cel)[nat], intensity(cel2)[nat], col="red")
points(intensity(cel)[!nat], intensity(cel2)[!nat], col="blue")
points(intensity(cel)[is.na(nat)], intensity(cel2)[is.na(nat)], pch="+")
legend(25000, 15000, c("PM","MM","Unknown","identity line"), c("red","blue","black","grey"), bg="white")
abline(0, 1, type="l", col="gray")
rm(nat)

abatch.n <- normalize(Dilution, method="constant", refindex=2)
plot(intensity(abatch.n[[1]]), intensity(abatch.n[[2]]), xlab="CEL file 1", ylab="CEL file 2",main="normalized by constant",sub="all probes plotted")
abline(0, 1, type="l", col="gray")


abatch.n <- normalize(Dilution, method="invariantset")
i.set <- history(abatch.n[[1]])$invariantset

plot(intensity(cel), intensity(cel2), xlab="CEL file 1", ylab="CEL file 2",main="raw values",sub="all probes plotted")

abline(0, 1, type="l", col="gray")
legend(25000,15000,c("invariant set","identity line","spline through the invariant set"),c("orange","grey","red"),bg="white")

plot(intensity(abatch.n[[1]]), intensity(abatch.n[[2]]), xlab="CEL file 1", ylab="CEL file 2",main="normalized by invariant set",sub="all probes plotted")
points(intensity(abatch.n[[1]])[i.set], intensity(abatch.n[[2]])[i.set], col="orange",pch=16)
abline(0, 1, type="l", col="gray")
legend(20000,10000,c("invariant set","identity line"),c("orange","grey"),bg="white")

##
## Normalization and its effects can be observed in a couple of commands
##

#par(mfrow=c(2,2))
#data(Dilution)
#boxplot(Dilution)
#boxplot(normalize(Dilution))


p <- probeset(Dilution, "1001_at")[[1]]

par(mfcol=c(5,2))
mymethods <- express.summary.stat.methods
nmet <- length(mymethods)
nc <- ncol(pm(p))

layout(matrix(c(1:nc, rep(nc+1, nc)), nc, 2), width = c(1, 1))
##p@pm <- log(p@pm)
##p@mm <- log(p@mm)

barplot(p)

results <- matrix(0, nc, nmet)
rownames(results) <- paste("sample", 1:nc)
colnames(results) <- mymethods

for (i in 1:nmet) {
  ev <- express.summary.stat(p, summary=mymethods[i], pmcorrect="pmonly")
  if (mymethods[[i]] != "medianpolish")
    results[, i] <- log(ev$exprs)
  else
    results[, i] <- ev$exprs
}
##matplot(results, matrix(1:3, nr=nc, nc=nmet), type="l", lty=1:3, col=rainbow(nmet),
##        xlab="expression value", lab="sample")
dotchart(results, labels=paste("sample", 1:nc))
##legend()

par(opar)
rm(opar)



