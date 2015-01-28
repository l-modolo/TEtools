### R code from vignette source 'independent_filtering_plots.Rnw'

###################################################
### code chunk number 1: setup
###################################################
options( width = 80 )


###################################################
### code chunk number 2: libraries
###################################################
library("genefilter")
library("ALL")
data("ALL")


###################################################
### code chunk number 3: sample_data
###################################################
bcell <- grep("^B", as.character(ALL$BT))
moltyp <- which(as.character(ALL$mol.biol) %in% 
                c("NEG", "BCR/ABL"))
ALL_bcrneg <- ALL[, intersect(bcell, moltyp)]
ALL_bcrneg$mol.biol <- factor(ALL_bcrneg$mol.biol)
n1 <- n2 <- 3
set.seed(1969)
use <- unlist(tapply(1:ncol(ALL_bcrneg), 
                     ALL_bcrneg$mol.biol, sample, n1))
subsample <- ALL_bcrneg[,use]


###################################################
### code chunk number 4: stats
###################################################
S <- rowSds( exprs( subsample ) )
temp <- rowttests( subsample, subsample$mol.biol )
d <- temp$dm
p <- temp$p.value
t <- temp$statistic


###################################################
### code chunk number 5: filter_volcano
###################################################
S_cutoff <- quantile(S, .50)
filter_volcano(d, p, S, n1, n2, alpha=.01, S_cutoff)


###################################################
### code chunk number 6: kappa
###################################################
t <- seq(0, 5, length=100)
plot(t, kappa_t(t, n1, n2) * S_cutoff, 
     xlab="|T|", ylab="Fold change bound", type="l")


###################################################
### code chunk number 7: table
###################################################
table(ALL_bcrneg$mol.biol)


###################################################
### code chunk number 8: filtered_p
###################################################
S2 <- rowVars(exprs(ALL_bcrneg))
p2 <- rowttests(ALL_bcrneg, "mol.biol")$p.value
theta <- seq(0, .5, .1)
p_bh <- filtered_p(S2, p2, theta, method="BH")


###################################################
### code chunk number 9: p_bh
###################################################
head(p_bh)


###################################################
### code chunk number 10: rejection_plot
###################################################
rejection_plot(p_bh, at="sample",
               xlim=c(0,.3), ylim=c(0,1000),
               main="Benjamini & Hochberg adjustment")


###################################################
### code chunk number 11: filtered_R
###################################################
theta <- seq(0, .80, .01)
R_BH <- filtered_R(alpha=.10, S2, p2, theta, method="BH")


###################################################
### code chunk number 12: R_BH
###################################################
head(R_BH)


###################################################
### code chunk number 13: filtered_R_plot
###################################################
plot(theta, R_BH, type="l",
     xlab=expression(theta), ylab="Rejections",
     main="BH cutoff = .10"
     )


###################################################
### code chunk number 14: sessionInfo
###################################################
si <- as.character( toLatex( sessionInfo() ) )
cat( si[ -grep( "Locale", si ) ], sep = "\n" )


